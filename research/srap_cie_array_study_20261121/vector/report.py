"""Offline summary, object-clustered paired intervals and fixed orthogonal slices."""
import numpy as np
import json,csv,time
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
GEOM={'plane_grid':'矩形面阵','hemi':'半球阵'}
FAM=['偏心壳','开缝框','分离多体','交叉棒薄层']
def read(p):return json.loads(p.read_text())
def fmt(a):return ' / '.join(f'{x:.3f}' for x in a)
def stats(rows,geom,bud,view):
    group=[r for r in rows if (geom=='all' or r['array']==geom) and r['budget']==str(bud) and r['view']==view]
    gs=['plane_grid','hemi'] if geom=='all' else [geom];d={}
    seq=dict()
    for metric in ['rel_L2','half_IoU','inside_L2','background_L2','data_rel','peak_chi','foreground_mean_ratio']:
        mat=np.array([[[float(next(r[metric] for r in group if r['array']==g and r['object']==str(i) and r['method']==m)) for m in ['Dual','SRAP']] for g in gs] for i in range(12)])
        delta=mat[:,:,0]-mat[:,:,1] if metric!='half_IoU' else mat[:,:,1]-mat[:,:,0]
        obj=delta.mean(axis=1);rng=np.random.default_rng(120613);bs=rng.integers(0,12,size=(10000,12))
        d[metric]=dict(medians=np.median(mat,axis=(0,1)).tolist(),median_paired_gain=float(np.median(obj)),CI=np.percentile(np.median(obj[bs],axis=1),[2.5,97.5]).tolist(),W=int(sum(obj>1e-7)),T=int(sum(abs(obj)<=1e-7)),L=int(sum(obj< -1e-7)))
        seq[metric]=mat
    return d,seq
def render():
    tick=time.perf_counter()
    with (ROOT/'objects.csv').open() as f:rows=list(csv.DictReader(f))
    cfg=read(ROOT/'FREEZE.json');gen=read(ROOT/'data/generation.json');arc=read(ROOT/'development_ring/transition.json');oldgen=read(ROOT/'development_ring/data/generation.json');acq=read(ROOT/'acquisition.json');sc=read(ROOT/'scoring_cost.json'); finer=read(ROOT/'data/shell_finer.json')
    receipts=[read(p) for p in sorted((ROOT/'predictions').glob('*/SEALED.json'))];assert len(receipts)==48 and len(rows)==192
    summ={f'{g}_{v}_{b}':stats(rows,g,b,v)[0] for g in ['all','plane_grid','hemi'] for v in ['raw','light'] for b in [60,180]}
    interaction={}
    for v in ['raw','light']:
        for b in [60,180]:
            mat=stats(rows,'all',b,v)[1]['rel_L2'];gain=mat[:,:,0]-mat[:,:,1];diff=gain[:,0]-gain[:,1]
            resamples=np.random.default_rng(120613).integers(0,12,(10000,12))
            interaction[f'{v}_{b}']=dict(median=float(np.median(diff)),CI=np.percentile(np.median(diff[resamples],axis=1),[2.5,97.5]).tolist(),object_differences=diff.tolist())
    summ['array_interaction']=interaction
    (ROOT/'summary.json').write_text(json.dumps(summ,indent=2))
    figdir=ROOT/'figures';figdir.mkdir(exist_ok=True)
    for geom in GEOM:
        for stem,indices in [('slices',[0,3,6,9]),('atlas_01',range(0,4)),('atlas_02',range(4,8)),('atlas_03',range(8,12))]:
            fig,axes=plt.subplots(4,9,figsize=(14.3,6.7),layout='constrained')
            for ir,idx in enumerate(indices):
                truth=np.load(ROOT/'sealed_truth'/f'{idx:02d}.npz')['chi']
                a=[truth]+[np.load(ROOT/'predictions'/f'{geom}_{idx:02d}_{m}'/'materials.npz')['light_180'] for m in ['Dual','SRAP']]
                for im,image in enumerate(a):
                    x=image.reshape(12,12,12)
                    for j,(cut,name) in enumerate(zip([x[6,:,:],x[:,6,:],x[:,:,6]],['YZ','XZ','XY'])):
                        ax=axes[ir,3*im+j];col=ax.imshow(cut.T,vmin=0,vmax=4,cmap='magma',origin='lower',interpolation='nearest',extent=[-.3,.3,-.3,.3]);ax.set_xticks([]);ax.set_yticks([])
                        if ir==0:ax.set_title(['Truth','Dual','SRAP'][im]+' '+name,fontsize=8)
                        if j==0 and im==0:ax.set_ylabel(['Shell','Frame','Bodies','Bars/layer'][idx//3]+f' #{idx:02d}',fontsize=9)
            fig.colorbar(col,ax=axes,shrink=.58,label='Contrast');fig.suptitle(f'{"Planar grid" if geom=="plane_grid" else "Hemispherical"} array · 180 s with physical correction · central voxel slices',fontsize=12)
            fig.savefig(figdir/f'{stem}_{geom}.png',dpi=175 if stem=='slices' else 148);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(9,5.8),layout='constrained')
    for gi,g in enumerate(GEOM):
        for vi,(v,b) in enumerate([('raw',180),('light',180)]):
            metric=stats(rows,g,b,v)[1]['rel_L2'][:,0,:];ax=axes[gi,vi]
            for i,(x,y) in enumerate(metric):ax.plot([0,1],[x,y],color=['#37739f','#b26939','#568d71','#7a5e94'][i//3],alpha=.72,lw=1.4,marker='o',markersize=3)
            ax.set_xticks([0,1],['Dual','SRAP']);ax.set_xlim(-.18,1.18);ax.set_ylabel('Relative material L2');ax.grid(axis='y',alpha=.15);ax.set_title(['Planar grid','Hemisphere'][gi]+' · '+['raw','physical correction'][vi],fontsize=10)
    fig.savefig(figdir/'paired_l2.png',dpi=180);plt.close(fig)
    a=summ['all_raw_180']['rel_L2'];b=summ['all_light_180']['rel_L2'];ss=[r for r in receipts if r['method']=='SRAP'];dd=[r for r in receipts if r['method']=='Dual']
    win=lambda v:f"{v['W']}/{v['T']}/{v['L']}"
    gains=lambda v:f"{v['median_paired_gain']:+.4f} [{v['CI'][0]:+.4f}, {v['CI'][1]:+.4f}]"
    raw_j='支持' if a['CI'][0]>0 else ('不利于' if a['CI'][1]<0 else '没有辨明')
    corr_j='支持' if b['CI'][0]>0 else ('不利于' if b['CI'][1]<0 else '没有辨明')
    txt=['# 三维双极化有限阵列确认实验','',f"12个独立材料分别由两个阵列测量。180秒预算下，SRAP 相对逐源 Dual 的材料相对 L2 配对收益为 {gains(a)}，共同物理校正后为 {gains(b)}。正值表示 SRAP 误差较低。两阵列先在同一物体内取平均，再对12个物体计数，原始胜/平/负 {win(a)}，校正后 {win(b)}。原始区间{raw_j} SRAP 的材料恢复优势，校正后的区间{corr_j}该优势。",'','## 阵列与预算','', '各行以12个物体为单位。误差和 IoU 分栏均列 Dual / SRAP 中位数。收益列是配对材料 L2 差的中位数及物体簇自助抽样95%区间（1万次）。','','|阵列|预算|读出|材料相对 L2|体积 IoU|配对收益及95%区间|胜/平/负|','|---|---:|---|---:|---:|---:|---:|']
    for g in GEOM:
        for bud in [60,180]:
            for v in ['raw','light']:
                u=summ[f'{g}_{v}_{bud}'];l=u['rel_L2'];txt.append(f"|{GEOM[g]}|{bud}s|{'原始' if v=='raw' else '物理校正'}|{fmt(l['medians'])}|{fmt(u['half_IoU']['medians'])}|{gains(l)}|{win(l)}|")
    interaction_text='，'.join(f"{'原始' if v=='raw' else '校正'} {interaction[v+'_180']['median']:+.4f} [{interaction[v+'_180']['CI'][0]:+.4f}, {interaction[v+'_180']['CI'][1]:+.4f}]" for v in ['raw','light'])
    txt+=['',f"180秒下同物体的面阵收益减半球收益为 {interaction_text}。正差表示路由在面阵的相对改善更大。阵列同时改变视角与场分量采样，几何效应限于这两种具体布置。"]
    txt+=['','![逐物体配对相对误差](figures/paired_l2.png)','','## 空间与场','', '体素支撑阈值为真值峰值的一半，两臂共用。内部和背景误差的分母均为完整材料真值的 L2 范数。数据残差为预测复场与实测复场的 Frobenius 相对误差。下表均为180秒物理校正读出。','','|阵列|内部 L2，Dual/SRAP|背景 L2|目标平均幅值比|峰值对比度|数据残差|','|---|---:|---:|---:|---:|---:|']
    for g in GEOM:
        u=summ[f'{g}_light_180'];txt.append('|'+GEOM[g]+'|'+ '|'.join(fmt(u[k]['medians']) for k in ['inside_L2','background_L2','foreground_mean_ratio','peak_chi','data_rel'])+'|')
    txt+=['','四类物体各取预定的第一个编号。每组三列为 YZ、XZ、XY 中央体素切片，切面坐标 +0.025λ，显示范围统一为 χ∈[0,4]。','','![矩形面阵固定切片](figures/slices_plane_grid.png)','','![半球阵固定切片](figures/slices_hemi.png)','','## 物理与算法边界','', '![矩形面阵与半球阵真实端口位置](figures/acquisition.png)','', '立方域宽0.6λ。8个物理发射位置各有两个正交单位偶极，16个物理接收位置各读两个正交分量，共512个复测量。两种布局最近端口到域中心均为0.9λ，矩形面阵外侧端口距中心略大；完整坐标与至立方体边界距离见 acquisition.json。矩形面阵含 2×4 发射、4×4 接收，左右错开，单位偶极沿全局 x、y。半球阵按等面积纬度和黄金角方位分布。方位布局与两个横向极化分别定义在代码中。各阵列使用同一物体和配对标准噪声方向，各物体间独立。噪声幅值为每个干净数据矩阵范数的10%。','', '正演为三维自由空间完整 dyadic Green、体元体积和 Clausius–Mossotti 加辐射反作用自项。24³细网以每体素64点体积分数生成观测，反演用12³，评分比较细网在同一粗体素内的均值。Mie 半径0.12λ，球 εr=3，比较平面波激励下181角度的绝对复矢量远场。相对误差随12³至24³由7.946%降至6.152%。偏心壳不同网格的完整复端口差见下表。进一步独立使用36³作参考，24³生成数据与36³参考的差见下列资源记录。','', '成熟逐源 Dual 从统一 χ=0.25 起步。接收 DPCS 截断为最大奇异值的0.1，三级傅里叶 mf=2,3,4，电流 β=6,4,3，材料 β=6,3,1。逐源 PR+ 用二次精确步长，跨三个分量和全部源拟合一个实标量 χ∈[0,4]。20步后检验源系数相对 RMS<0.001，每级上限500，沿用三维 R0（无额外材料先验）。两臂每16步或阶段末都用原方程重算场并审计，其间的材料不参与选择。','', 'SRAP 同期在 X=J/6、Q 空间上依次做材料投影、仿射图与数据球投影，再将两个增量投到接收弱且传播强的全谱子空间，重新拟合合法投影傅里叶系数并重启 PR。两组矢量均用相同路由。数据球半径为 0.1‖Y‖/√1.01，单乘子约束所有端口。独立于目标的接收谱阈值为0.1，传播能量阈值为0.01 λmax(G*G)。每臂现场重建；完整5184维稠密 Hermitian 分解仅用于路由。图投影的 H=I+36G*G 使用零填充 FFT-CG，32维接收 Schur 系统求共享乘子；物理前向和伴随为 FFT-BiCGSTAB。底层复数为 complex64，接收基、合法回注的小矩阵使用 complex128。','', '原始预算读取同一冷轨迹的60秒/180秒前缀。light 从48秒/144秒边界前公开原物理目标最小值出发，分别在独立12秒/36秒窗口完成最多20个 L-BFGS-B 接受步，按完成时刻选择。Phi 为逐源相对复数据平方误差均值。跨过截止的审计和线搜索试探无法入选。两种预算共享相同前缀，程序与CUDA上下文启动耗时在臂外，采集算子、观测初始化、正式路由构建、投影和正反演计入臂内。评分在48个预测包均验证并封存后单独进行。','','## 实际资源','']
    shape_note="在180秒校正切片中，偏心壳00在两阵列的两臂均呈实心团块；其粗网真值有一个闭合空腔，四个恢复体均为零。开缝框03在半球XY面可辨开口，在矩形面阵中的响应幅值较弱。三个多体对象的粗网真值各有三个分离支撑，两个阵列、两臂的整幅恢复体均无体素达到固定半峰阈值。交叉棒09在半球阵找到了主要响应位置，细杆端部和薄层仍有缺失。面阵有8/12个Dual、7/12个SRAP的体积IoU为零。其物理场残差均接近10%，数据拟合可以与轮廓丢失并存。"
    atlas_text=['',shape_note,'', '### 全部12个物体的切片图册','', '每页四个物体；每行从左到右依次为真值、Dual、SRAP，各给三个方向的中央体素切片。完整支撑的三维交并比、六邻接连通分量和闭合空腔计数见逐物体表。','']
    for g in GEOM:
        links=' · '.join(f'[{4*j:02d}–{4*j+3:02d}](figures/atlas_{j+1:02d}_{g}.png)' for j in range(3))
        atlas_text.append(f"- {GEOM[g]} {links}")
    at=txt.index('## 物理与算法边界');txt[at:at]=atlas_text+['']
    for g in GEOM:
        tx=acq[g]['transmit'];rx=acq[g]['receive']
        txt.append(f"- {GEOM[g]}端口距中心范围（发射/接收）{fmt(tx['center_range'])}λ 和 {fmt(rx['center_range'])}λ，距实际立方体边界范围 {fmt(tx['cube_boundary_range'])}λ 和 {fmt(rx['cube_boundary_range'])}λ。两布局最近端口距包围球均为0.3804λ。")
    paid_ring=[read(p) for p in (ROOT/'development_ring/predictions').glob('*/SEALED.json')]
    txt.append(f"- 开发圆环支路有 {len(paid_ring)} 个已封存臂，累计墙钟 {sum(r['total_s'] for r in paid_ring):.2f} 秒，没有打开其材料真值评分。V1 生成圆环/半球24个观测与Mie合计 {oldgen['wall']:.2f} 秒，包含在V1内的保留半球工作只计一次。V1控制器停止前主计算进程耗时 {max(p['wall_since_start'] for p in arc['proc']):.2f} 秒，含已完成半球和圆环、启动与中断在途的半球；此项与分臂账目重叠。V1平面圆环单独路由预检额外耗时3.09秒。")
    for rec in finer['rows']:
        txt.append(f"- {GEOM[rec['array']]}偏心壳 12³→24³ / 24³→36³ 的复端口差分别为 {100*rec['relative_12_24']:.3f}% / {100*rec['relative_24_36']:.3f}%。三网校核单独开销为两阵合计 {finer['wall']:.2f} 秒。")
    for title,subs in [('Dual',dd),('SRAP',ss)]:
        raw=np.array([r['raw_s'] for r in subs]);tot=sum(r['total_s'] for r in subs);cp=sum(r['cpu'] for r in subs);stops=[sum(r['raw_s']>v for r in subs) for v in [60,180]]
        txt.append(f"- {title} 共24臂，原始轨迹中位 {np.median(raw):.2f} 秒，范围 {min(raw):.2f}–{max(raw):.2f} 秒，跨60/180秒的轨迹分别为 {stops[0]}/{stops[1]}。含两个独立校正的实付墙钟 {tot:.1f} 秒，进程 CPU {cp:.1f} 秒。")
    for arm,subset in [('Dual',dd),('SRAP',ss)]:
        counts={}
        for r in subset:
            for counter in [r['raw_counts']]+[k['counts'] for k in r['correctors']]:
                for key,value in counter.items():counts[key]=counts.get(key,0)+value
        txt.append(f"- {arm} 总计 {counts.get('PR_steps',0)} 个逐源PR更新，物理正演/伴随 RHS {counts.get('forward_solve_rhs',0)} / {counts.get('adjoint_solve_rhs',0)}，G及G*的单源作用数 {counts.get('G_rhs',0)} / {counts.get('Gadj_rhs',0)}，图CG迭代 {counts.get('H_cg_iterations',0)}。稠密路由谱与主特征值的稠密乘法另计，不归入FFT作用数。")
    for g in GEOM:
        subset=[r for r in ss if r['array']==g];r=subset[0];st=[rr['route']['setup_s'] for rr in subset]
        txt.append(f"- {GEOM[g]}接收强秩 {r['rank']}，路由秩 {r['route']['route_rank']}/5184，完整路由及图预计算费用范围 {min(st):.2f}–{max(st):.2f} 秒/臂。")
    proj=[p for r in ss for p in r['projection_log']];allcorr=[c for r in receipts for c in r['correctors']]
    txt+= [f"- 设备为单张 RTX 4060 Laptop 8 GB。PyTorch统计的CUDA张量分配峰值为 {max(r['peak_gpu_bytes'] for r in receipts)/2**30:.2f} GiB（V2各臂重置；继承的半球V1记录为进程累计高水位，不含驱动分配）。投影调用 {len(proj)} 次，图CG最大相对残差 {max(p['cg_residual'] for p in proj):.2e}，投影后最大数据残差/噪声球半径 {max(p['data_over_delta'] for p in proj):.5f}。96条校正分支中 {sum(c['time_exceeded'] for c in allcorr)} 条触及墙钟截止，接受步范围 {min(c['num_accepted'] for c in allcorr)}–{max(c['num_accepted'] for c in allcorr)}。",f"- 新面阵12组数据生成的增量墙钟 {gen['wall']:.2f} 秒，CPU {gen['cpu']:.2f} 秒。正式评分墙钟 {sc['wall']:.2f} 秒，CPU {sc['cpu']:.2f} 秒。各审计的FFT作用、物理正/伴随 RHS、BiCGSTAB 与图CG迭代及校正尝试均保存在逐臂封存记录。",'','## 逐物体表','', '完整192行含两种预算、两类读出、峰值、目标/背景误差、连通分量、封闭空腔、实际计费与数据残差。 [下载 CSV](objects.csv)','', '|物体|阵列|180秒 L2 原始 Dual/SRAP|180秒 L2 校正|180秒 IoU 校正|180秒数据残差校正|','|---|---|---:|---:|---:|---:|']
    for i in range(12):
        for g in GEOM:
            query=lambda v,k: [float(next(r[k] for r in rows if r['object']==str(i) and r['array']==g and r['budget']=='180' and r['view']==v and r['method']==m)) for m in ['Dual','SRAP']]
            txt.append(f"|{i:02d} {FAM[i//3]}|{GEOM[g]}|{fmt(query('raw','rel_L2'))}|{fmt(query('light','rel_L2'))}|{fmt(query('light','half_IoU'))}|{fmt(query('light','data_rel'))}|")
    txt+=['','## 来源与复现','', '成熟逐源更新可对照既有项目的 research_max_struct_r9/src/vector_original_cuda.py 第298–342行。三级双β、起点、系数停止和 R0 对照 research_max_struct_r11/configs/vector_confirmation_gpu_recovery_attempt1_execution.json 第63–100行。早期配置的对比度上界为8，本实验按题设改为4。',  '在安装 PyTorch CUDA、NumPy 与 SciPy 的显卡环境执行 generate.py、invert.py、score.py。report.py 另需 Matplotlib。三个独立进程分别生成、反演及评分。所有输入和选定材料已保存，可从评分或作图阶段单独开始。生成脚本内有私有物体设计，反演脚本只导入 physics.py 和 projection.py。']
    (ROOT/'REPORT_ZH.md').write_text('\n'.join(txt)+'\n')
    (ROOT/'report_cost.json').write_text(json.dumps(dict(wall=time.perf_counter()-tick)))
if __name__=='__main__':render()
