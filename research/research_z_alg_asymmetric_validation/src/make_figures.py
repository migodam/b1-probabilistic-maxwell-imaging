"""Source-only campaign figures. No Book, forward calls, or network operations."""
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from PIL import Image, ImageOps, ImageDraw

P = Path(__file__).resolve().parents[1]
OUT = P / 'figures'
COLORS = {'A':'#2857a5','B':'#d57622','F':'#73818c','bad':'#b7373f'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                     'axes.spines.right':False,'savefig.dpi':180,'pdf.fonttype':42,
                     'axes.titlesize':11,'figure.titlesize':14})
SOURCES, ARRAYS, MANIFEST = [], {}, []


def read(relative):
    SOURCES.append(relative)
    return json.loads((P/relative).read_text())


def save(fig, name, caption, rows, scope):
    fig.savefig(OUT/(name+'.png'), bbox_inches='tight', facecolor='white')
    fig.savefig(OUT/(name+'.pdf'), bbox_inches='tight', facecolor='white')
    plt.close(fig)
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with (OUT/'source_data'/(name+'.csv')).open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=columns);writer.writeheader()
        writer.writerows({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in row.items()} for row in rows)
    MANIFEST.append(dict(name=name,caption_zh=caption,scope=scope,source_csv='source_data/'+name+'.csv'))


def main():
    if not (P/'CLOSURE.json').exists():raise RuntimeError('Plot only after campaign closure')
    if OUT.exists() and any(OUT.iterdir()):raise FileExistsError('Preserve existing figure outputs')
    OUT.mkdir(exist_ok=True);(OUT/'source_data').mkdir()
    contract=read('EXPERIMENT_CONTRACT.json')
    e1={(p,m):read(f'results/E1/P{p}/{m}/SUMMARY.json') for p in range(3) for m in ['A','B']}
    e2=read('results/E2/SUMMARY.json');cost=read('results/COST_CONTROL/SUMMARY.json')
    local=read('results/E4_local/SUMMARY.json');target=read('results/E4_local/CONDITIONAL_TARGET.json')
    legacy=read('audit/E4_legacy_replay/e4_two_cell_posterior.json')
    certs=read('audit/INTERVAL_CHART_CERTIFICATES.json');basis=read('results/E3/SUMMARY.json')

    # 1: exact frozen geometry and all three deterministic reference materials.
    fig,ax=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
    c=contract['E1'];centers=np.array(c['centers']);src=np.array(c['sources']);rec=np.array(c['receivers'])
    for i,(x,y) in enumerate(centers):
        ax[0].add_patch(Circle((x,y),contract['disk_radius'],color='#c6d6ed',ec=COLORS['A']))
        ax[0].text(x,y+.075,f'cell {i}',ha='center')
    ax[0].scatter(*src.T,marker='*',s=130,color='#aa4a37',label='Fixed sources')
    ax[0].scatter(*rec.T,marker='v',s=75,color='#28735a',label='Fixed receiver')
    for i,(x,y) in enumerate(src):ax[0].text(x+.035,y,f'S{i}')
    for i,(x,y) in enumerate(rec):ax[0].text(x+.035,y,f'R{i}')
    ax[0].set(xlim=(-.8,.8),ylim=(-.9,1.03),xlabel='x / wavelength',ylabel='y / wavelength',title='Asymmetric fixed-label Green VIE')
    ax[0].set_aspect('equal');ax[0].legend(loc='lower center',fontsize=8)
    refs=np.array(c['public_reference_materials']);ARRAYS['reference_materials']=refs
    im=ax[1].imshow(refs,vmin=.5,vmax=4,cmap='viridis',aspect='auto')
    for i in range(3):
        for j in range(6):ax[1].text(j,i,f'{refs[i,j]:.1f}',ha='center',va='center',color='white' if refs[i,j]<2.3 else '#222')
    ax[1].set(xticks=range(6),yticks=range(3),yticklabels=['P0','P1','P2'],xlabel='Cell label',title='Three public deterministic material contexts')
    fig.colorbar(im,ax=ax[1],label='Real material contrast')
    fig.suptitle('6 real material variables; 2 sources x 1 complex receiver = 4 real data')
    rows=[dict(kind='cell',label=i,x=x,y=y,radius=contract['disk_radius']) for i,(x,y) in enumerate(centers)]
    rows += [dict(kind=kind,label=i,x=x,y=y) for kind,points in [('source',src),('receiver',rec)] for i,(x,y) in enumerate(points)]
    rows += [dict(kind='reference_material',parent=i,cell=j,value=refs[i,j]) for i in range(3) for j in range(6)]
    save(fig,'01_geometry','非对称固定标签小系统与三个公开参考材料。6实材料、4实观测；不是Full66。',rows,'New six-cell contract; three deterministic references; no population inference')

    # 2: include all material-valued failures, not only successful endpoints.
    fig,axes=plt.subplots(1,3,figsize=(13,4),sharey=True,layout='constrained');rows=[]
    for p,ax in enumerate(axes):
        for method in ['A','B']:
            for qualified in [True,False]:
                rs=[r for r in e1[p,method]['rows'] if r.get('material_RMS') is not None and r.get('pair_relative') is not None and r['qualified']==qualified]
                ax.scatter([r['material_RMS'] for r in rs],[max(r['pair_relative'],1e-18) for r in rs],
                           s=18 if qualified else 45,marker='o' if qualified else 'x',alpha=.6,
                           color=COLORS[method] if qualified else COLORS['bad'],
                           label=(('Direct' if method=='A' else 'Current')+' qualified') if qualified else ('Failed' if method=='A' else None))
            rows += [{**r,'parent':p,'method':method} for r in e1[p,method]['rows']]
        ax.axhline(1e-10,color='#555',ls='--',lw=1,label='Response gate' if p==0 else None)
        ax.set(yscale='log',xlabel='Material RMS from reference',title=f'P{p}: sampled grid + radial endpoints')
    axes[0].set_ylabel('Independent full-response relative error');axes[0].legend(fontsize=8)
    fig.suptitle('Distinct materials, fixed response: finite sampled endpoints; failures retained')
    save(fig,'02_material_response','三参考上下文的材料位移与独立全前向响应误差。叉号保留有材料端点的失败；无合法端点失败见覆盖图。',rows,'No continuous guarantee between samples; no global fiber coverage')

    # 3: exact matching target keys, including failed free targets.
    fig,axes=plt.subplots(1,3,figsize=(12,4.2),layout='constrained');rows=[]
    for p,ax in enumerate(axes):
        aa={r['key']:r for r in e1[p,'A']['rows']};bb={r['key']:r for r in e1[p,'B']['rows']}
        for key,a in aa.items():
            b=bb[key];same=a['qualified']==b['qualified'];qualified=a['qualified'] and b['qualified']
            ax.scatter(*a['free_values'],s=19 if qualified else 60,marker='o' if qualified else 'x',
                       color='#3d867b' if qualified else COLORS['bad'] if same else '#b066ba')
            rows.append(dict(parent=p,key=key,kind=a['kind'],free0=a['free_values'][0],free1=a['free_values'][1],A_qualified=a['qualified'],B_qualified=b['qualified']))
        qa=e1[p,'A']['qualified'];attempt=e1[p,'A']['attempts']
        ax.set(xlabel='First selected free material coordinate',ylabel='Second selected free coordinate',title=f'P{p}: A = B, {qa}/{attempt} qualified\n49/49 fixed-grid targets each')
    fig.suptitle('Paired coverage is identical: green = both qualified; red = both failed')
    save(fig,'03_paired_coverage','同一目标逐项配对：两方法覆盖一致，失败未删除。固定grid均49/49；radial边界终止不代表全域边界认证。',rows,'Sampled target coverage only; no additional branches or global completeness')

    # 4: both windows, no substitution of solver time for total time.
    fig,ax=plt.subplots(1,3,figsize=(14,4.8),layout='constrained');rows=[]
    historical_audit=read('audit/results_independent/AUDIT.json')
    main_times={m:np.array([next(r[m+'_job_wall'] for r in historical_audit['E1_contexts'] if r['parent']==p) for p in range(3)]) for m in ['A','B']}
    x=np.arange(3);w=.35
    for m,shift in [('A',-w/2),('B',w/2)]:ax[0].bar(x+shift,main_times[m],w,label='Direct Full-J' if m=='A' else 'Current',color=COLORS[m])
    main_saving=1-main_times['B'].sum()/main_times['A'].sum()
    ax[0].set(xticks=x,xticklabels=['P0','P1','P2'],ylabel='Complete controlled job wall (s)',title=f'Original exploratory window\nPooled saving {100*main_saving:.1f}%')
    ax[0].legend(fontsize=8)
    methods=['DYNAMIC_FULLJ','FROZEN_J','CURRENT_FULL_KERNEL'];labels=['Dynamic\nFull-J','Frozen J','Current\nfull kernel']
    pooled={}
    for method in methods:
        rr=[r for r in cost['rows'] if r['method']==method]
        pooled[method]={k:sum(r[k] for r in rr) for k in ['solver_seconds','validation_seconds','persistence_seconds','total_seconds','RHS']}
        # Setup was stored outside each row's total_seconds. Include it in both
        # complete total and solver subtotal, matching the frozen grouped report.
        setup=sum(r['setup_seconds'] for r in rr)
        pooled[method]['solver_seconds']+=setup
        pooled[method]['total_seconds']+=setup
        rows += [dict(window='optimized_equal_trace',**r) for r in rr]
    bottoms=np.zeros(3)
    for field,color,label in [('solver_seconds','#4b74b5','Solver + setup'),('validation_seconds','#75ab9d','Independent validation'),('persistence_seconds','#c5b9a8','Persistence')]:
        values=[pooled[m][field] for m in methods];ax[1].bar(range(3),values,bottom=bottoms,color=color,label=label);bottoms+=values
    total_saving=1-pooled[methods[2]]['total_seconds']/pooled[methods[0]]['total_seconds']
    sv_saving=1-sum(pooled[methods[2]][k] for k in ['solver_seconds','validation_seconds'])/sum(pooled[methods[0]][k] for k in ['solver_seconds','validation_seconds'])
    ax[1].set(xticks=range(3),xticklabels=labels,ylabel='Pooled controlled wall (s)',title=f'Optimized equal-trace window\nTotal saving {100*total_saving:.1f}%: FAIL 20% gate')
    ax[1].legend(fontsize=8)
    for method,color in zip(methods,[COLORS['A'],COLORS['F'],COLORS['B']]):
        gg=[r for r in cost['grouped'] if r['method']==method];ax[2].plot([r['parent'] for r in gg],[100*(1-r['total_seconds']/next(z['total_seconds'] for z in cost['grouped'] if z['parent']==r['parent'] and z['method']=='DYNAMIC_FULLJ')) for r in gg],'-o',color=color,label=method.replace('_',' '))
    ax[2].axhline(20,color=COLORS['bad'],ls='--',label='20% total-cost gate');ax[2].axhline(0,color='#888',lw=.7)
    ax[2].set(xticks=range(3),xticklabels=['P0','P1','P2'],ylabel='Total saving vs Dynamic Full-J (%)',title=f'Optimized context-level totals\nSolver + validation alone: {100*sv_saving:.1f}% pooled')
    ax[2].legend(fontsize=7)
    for p in range(3):
        for m in ['A','B']:rows.append(dict(window='original_exploratory',parent=p,method=m,controlled_job_wall=main_times[m][p]))
    ARRAYS['original_window_times_A_B']=np.stack([main_times['A'],main_times['B']])
    fig.suptitle('Same coverage; no net-total-cost GO. Three deterministic contexts, no population CI.')
    save(fig,'04_cost_windows','两个费用窗口分别保留。优化equal-trace总费用仅约13.6%改善，未过20%门；solver加validation约28.6%不替代总费用。FrozenJ同时展示。',rows,'Timing repeats are not independent parents; original and optimized windows not merged')

    # 5: starts + retraced cycles; no second root and no uniqueness theorem.
    fig,ax=plt.subplots(1,3,figsize=(13,4.1),layout='constrained');rows=[]
    for r in e2['starts']:
        ax[0].scatter(r['start'],r['iterations'],marker='o' if r['qualified'] else 'x',s=35,color='#3d867b' if r['qualified'] else COLORS['bad'])
        rows.append(dict(type='multistart',**r))
    nqual=sum(r['qualified'] for r in e2['starts'])
    ax[0].set(xlabel='Start index (reference + 32 starts)',ylabel='Corrector iterations',title=f'{nqual}/{len(e2["starts"])} qualified; 1 root found\nRed crosses retain failed starts')
    for cycle in e2['cycles']:
        nodes=cycle['nodes'];t=[0]+[r['t'] for r in nodes];rms=[0]+[r['material_RMS'] for r in nodes]
        ax[1].plot(t,rms,'-o',ms=2.8,label=f'Cycle {cycle["cycle"]}, return {cycle["return_RMS"]:.1e}')
        ax[2].semilogy([r['t'] for r in nodes],[r['pair_relative'] for r in nodes],'-o',ms=2.8,label=f'Cycle {cycle["cycle"]}')
        rows += [dict(type='data_path',cycle=cycle['cycle'],**r) for r in nodes]
    ax[1].set(xlabel='Auxiliary path time',ylabel='Material RMS from starting root',title='Two changing-data out-and-back cycles');ax[1].legend(fontsize=8)
    ax[2].set(xlabel='Auxiliary path time',ylabel='Full response error to each moving target',title='Paths retrace to the original root');ax[2].legend(fontsize=8)
    fig.suptitle('E2: no isolated second branch found; NOT_FOUND is not a uniqueness proof')
    save(fig,'05_branch_search','33起点中31个合格、仅找到1根，2条数据往返路径回到原根。循环途中目标变化，非固定数据纤维，不声称monodromy。',rows,'No bounded-domain uniqueness or branch-completeness certificate')

    # 6: different probability objects on different axes.
    fig,ax=plt.subplots(1,3,figsize=(13,4.3),layout='constrained');rows=[]
    z=legacy['results']['32'];vals=[z['original_box']['event_x1_gt_x2'],z['ordered_quotient']['event_x1_gt_x2'],z['ordered_quotient']['uniform_branch_wrong_control']['event_x1_gt_x2']]
    ax[0].bar(range(3),vals,color=[COLORS['A'],COLORS['B'],COLORS['bad']]);ax[0].set(xticks=range(3),xticklabels=['Original\nmaterial','Correct\nquotient','Wrong uniform\nbranches'],ylim=(.45,.72),ylabel='P(x1 > x2 | noisy data)',title='Legacy symmetric two-cell posterior')
    for i,v in enumerate(vals):ax[0].text(i,v+.006,f'{v:.8f}',ha='center',fontsize=8)
    rows += [dict(object='legacy_noisy_posterior',method=m,event_probability=v) for m,v in zip(['original','quotient','wrong_uniform'],vals)]
    q=local['orders']['8'];ms=['A','B','wrong_uniform_free','wrong_equal_nodes'];vs=[q[m]['event_probability'] for m in ms]
    ax[1].bar(range(4),vs,color=[COLORS['A'],COLORS['B'],COLORS['bad'],'#ab6173']);ax[1].set(xticks=range(4),xticklabels=['Direct\ncoarea','Current\ncoarea','Uniform\nfree','Equal\nnodes'],ylim=(.499,.5055),ylabel='P(first free material > reference | local tile)',title='New six-cell LOCAL conditional kernel')
    for i,v in enumerate(vs):ax[1].text(i,v+.00015,f'{v:.8f}',ha='center',fontsize=7.5)
    for order,d in local['orders'].items():
        for m in ms:rows.append(dict(object='new_local_conditional',order=int(order),method=m,**d[m]))
    orders=[int(k) for k in local['orders']]
    for metric in ['mean','event_probability']:
        discrepancies=[np.max(np.abs(np.array(local['orders'][str(o)]['A'][metric])-np.array(local['orders'][str(o)]['B'][metric]))) for o in orders]
        ax[2].semilogy(orders,discrepancies,'-o',label='Mean max' if metric=='mean' else 'Event abs')
    ax[2].set(xticks=orders,xlabel='Gauss-Legendre order',ylabel='Direct/current absolute discrepancy',title='Numerical agreement; no global mass bound');ax[2].legend(fontsize=8)
    fig.suptitle('Correct weights matter. Legacy posterior and new tile-conditional law are distinct.')
    save(fig,'06_probability','旧二单元有限噪声后验的正确事件概率约.67516447，错误均匀分支为.5；新六单元只计算局部tile条件核，正确约.50395231，错误权重为.5。两者不可互换。',rows,'New local tile mass is an unnormalized coarea density integral, not a probability or noisy evidence; global missing mass unbounded')

    # 7: full basis roundtrip versus support-changing projections, from saved arrays.
    SOURCES.append('results/E3/COMPACT_ARRAYS.npz')
    fig,ax=plt.subplots(1,2,figsize=(11,4.4),layout='constrained');rows=[]
    with np.load(P/'results/E3/COMPACT_ARRAYS.npz',allow_pickle=False) as saved:
        for i,name in enumerate(basis['transforms']):
            full=np.sqrt(np.mean(saved[name+'_roundtrip_error']**2,axis=1));trunc=np.sqrt(np.mean(saved[name+'_truncated3_error']**2,axis=1))
            # Native full roundtrip is exactly zero; zeros stay visible on symlog.
            ax[0].scatter(np.full(len(full),i),full,s=8,alpha=.16,color=COLORS['A'])
            ax[1].boxplot([trunc],positions=[i],widths=.5,showfliers=False,patch_artist=True,boxprops=dict(facecolor='#d7ba9a'))
            ax[1].scatter(np.full(len(trunc),i),trunc,s=5,alpha=.10,color=COLORS['B'])
            ARRAYS[name+'_full_roundtrip_RMS']=full;ARRAYS[name+'_trunc3_RMS']=trunc
            rows += [dict(transform=name,source_occurrence=j,full_roundtrip_RMS=float(a),truncated3_material_RMS=float(b)) for j,(a,b) in enumerate(zip(full,trunc))]
    for a in ax:a.set(xticks=range(3),xticklabels=['Native','DFT real\nconjugacy','Haar6'])
    ax[0].set(yscale='symlog',ylabel='Complete basis material RMS error',title='Complete invertible maps preserve material');ax[0].set_yscale('symlog',linthresh=1e-16)
    ax[1].set(ylabel='Rank-three projected material RMS error',title='Three retained coefficients change support')
    fig.suptitle('728 saved endpoint occurrences; same original prior pushed through complete bases')
    save(fig,'07_representation','完整正交基只换表示，原prior经推前保持；三系数投影改变材料与支持。728条封存端点含重复，不是独立样本或后验抽样。',rows,'No forward responses computed for truncation; no fiber/gauge or Full66 extension claim')

    # 8: keep failed broad certificate attempts, numeric tile versus narrow proof.
    fig,ax=plt.subplots(1,3,figsize=(13,4.4),layout='constrained');rows=[]
    for p,a in enumerate(ax):
        cert=certs[p];attempts=cert['attempts'];widths=[r['halfwidth'] for r in attempts]
        a.loglog(widths,[max(r['mapping_inclusion_ratios']) for r in attempts],'-o',label='Max inclusion ratio')
        a.loglog(widths,[max(r['contraction_ratios']) for r in attempts],'-s',label='Max contraction ratio')
        a.axhline(1,color=COLORS['bad'],ls='--',label='Strict < 1 needed')
        final=attempts[-1];a.axvline(final['halfwidth'],color='#378377',ls=':',lw=2)
        a.axvline(target['halfwidth'],color='#9b803e',ls=':',lw=2)
        a.set(xlabel='Free-coordinate box halfwidth',ylabel='Interval bound ratio',title=f'P{p}: certified halfwidth {final["halfwidth"]:g}')
        rows += [dict(parent=p,**r) for r in attempts]
    ax[0].legend(fontsize=7);fig.suptitle('Frozen binary64 coefficient certificates are narrower than the 0.05 numerical tile')
    save(fig,'08_certificate_scope','保留所有较宽区间失败尝试：仅冻结binary64系数的窄局部graph获证；.05数值积分tile没有整体区间证书。真实Hankel系数包围、连续PDE与全局质量均未认证。',rows,'Exact frozen-coefficient local graph only; no global missing-mass bound; failed bounds are not nonexistence proofs')

    np.savez_compressed(OUT/'source_data'/'COMPACT_PLOT_ARRAYS.npz',**ARRAYS)
    manifest=dict(status='RENDERED_REVIEW_REQUIRED',campaign_closed=True,new_RHS=0,Book_used=False,
                  source_paths=sorted(set(SOURCES)),figure_count=len(MANIFEST),figures=MANIFEST,
                  scientific_scope='Finite asymmetric local 2D graphs; same coverage; net total-cost NO_GO; no second isolated branch found; no Full66 extension',
                  statistics='Only three deterministic reference contexts; no population CI')
    (OUT/'MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    (OUT/'CAPTIONS_ZH.md').write_text('\n\n'.join(f'### {r["name"]}\n\n{r["caption_zh"]}\n\nScope: {r["scope"]}' for r in MANIFEST)+'\n')
    thumbs=[]
    for r in MANIFEST:
        im=Image.open(OUT/(r['name']+'.png')).convert('RGB');im.thumbnail((1050,430))
        tile=Image.new('RGB',(1100,480),'white');tile.paste(im,((1100-im.width)//2,28));ImageDraw.Draw(tile).text((15,8),r['name'],fill='black');thumbs.append(tile)
    sheet=Image.new('RGB',(2200,1920),'#e4e6e8')
    for i,tile in enumerate(thumbs):sheet.paste(tile,((i%2)*1100,(i//2)*480))
    sheet.save(OUT/'CONTACT_SHEET.png')
    print(json.dumps(dict(figures=len(MANIFEST),png=8,pdf=8,output=str(OUT),new_RHS=0)))


if __name__=='__main__':main()
