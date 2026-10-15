#!/usr/bin/env python3
"""Read existing arrays only: no physics backend, solve, training or search."""
import csv
import json
import time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'research/prior_constrained_20261015'
FIG = OUT/'figures'
TROOT = ROOT.parent/'b1-research-t/research/research_t_20261015'
A30 = ROOT.parents[1]/'A30/a30-maxwell-material-manifold'
A26 = ROOT.parents[1]/'A26/a26-om-feedback-prior-completion'
V = (1+.2j)/np.sqrt(1.04)
SAVED = json.loads((OUT/'GALLERY_DATA.json').read_text()) if (OUT/'GALLERY_DATA.json').exists() else None


def load_chi(path):
    if not path.exists() and SAVED:
        for r in SAVED['E0']:
            if str(path).endswith(r['path']):
                a=r['arrays']; return np.asarray(a['epsilon_real'])-1+1j*np.asarray(a['epsilon_imag'])
        for r in SAVED['E1']:
            for name,end in [('anchor_array','A'),('candidate_array','B')]:
                if str(path).endswith(r['pair'][name]):
                    a=r['arrays']; return np.asarray(a['epsilon_real_'+end])-1+1j*np.asarray(a['epsilon_imag_'+end])
    with np.load(path, allow_pickle=False) as z:
        return np.asarray(z['chi']).reshape(24,24)


def main():
    started = time.perf_counter()
    FIG.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.titlesize':11,'figure.facecolor':'white','savefig.facecolor':'white'})
    config_path=A26/'configs/a25/campaign.json'
    cfg=json.loads((config_path if config_path.exists() else OUT/'source_snapshots/campaign.json').read_text())
    geometry_path=A30/'inputs/a29_release/inputs/a25/geometry.json'
    points=np.asarray(json.loads((geometry_path if geometry_path.exists() else OUT/'source_snapshots/geometry.json').read_text())['points'])
    def bg(name):
        b=np.full(576,complex(*cfg['background_weak']),complex)
        if name=='ring':
            r=np.linalg.norm(points,axis=1)
            b[(r>=cfg['ring_radius_inner'])&(r<=cfg['ring_radius_outer'])]=complex(*cfg['background_ring'])
        return b.reshape(24,24)
    audit=list(csv.DictReader((OUT/'delegated/cached_audit/E0_AUDIT.csv').open()))
    decision=json.loads((OUT/'E1_DECISION.json').read_text())
    pairs=json.loads((OUT/'WITNESS_PAIRS.json').read_text())['pairs']
    replay=list(csv.DictReader((OUT/'PAIR_NOISE_REPLAY.csv').open()))
    search=list(csv.DictReader((OUT/'E1_SEARCH.csv').open()))
    data={'scope':'cached material-only visualization; no truth used for selection/search',
          'new_physics':{'states':0,'rhs':0,'factorizations':0,'training':0},
          'formula':{'epsilon':'1+chi','rho':'Re(conj(v)*(chi-background))','rho_orth':'Im(conj(v)*(chi-background))','v':[float(V.real),float(V.imag)]},
          'grid':{'shape':[24,24],'display':'array cell index; origin lower; no interpolated pixels','points':points.tolist()},
          'sources':{'protocol':'PROTOCOL_FROZEN.json','E0':'E0_DECISION.md','E1':'E1_DECISION.json','cached':'delegated/cached_audit/E0_AUDIT.csv','pairs':'WITNESS_PAIRS.json','replay':'PAIR_NOISE_REPLAY.csv'},
          'E0':[],'E1':[]}
    paths=[TROOT/'qualification/p250065_ring_a1_snr20/original_n3.npz',TROOT/'qualification/p250065_ring_a1_snr20/qualified_minchange_base.npz',TROOT/'mechanism_v1/p250065_ring_a1_snr20/edit_00/new_valid_candidate.npz']
    ids=['N3','T_QUAL:qualified_minchange_base.npz','mechanism_v1:edit_00:new_valid_candidate.npz']
    names=['Original N3\nHistorical strong Top-1','minchange_base\nHistorical noise-ball feasible','T v1 edit 00 (compensated)\nHistorical noise-ball feasible']
    chis=[load_chi(p) for p in paths]
    layers=[[1+x.real for x in chis],[x.imag for x in chis],[(np.conj(V)*(x-bg('ring'))).imag for x in chis]]
    limits=[(min(x.min() for x in layers[0]),max(x.max() for x in layers[0])),(min(x.min() for x in layers[1]),max(x.max() for x in layers[1]))]
    om=max(abs(x).max() for x in layers[2]);limits.append((-om,om))
    fig,axs=plt.subplots(3,3,figsize=(12.7,10.3),layout='constrained')
    fig.suptitle('Same historical observation: p250065_ring_a1_snr20\nNoise-ball feasibility does not establish original-prior support',fontsize=16)
    for j,(p,ident) in enumerate(zip(paths,ids)):
        row=next(r for r in audit if r['case_id']=='p250065_ring_a1_snr20' and r['candidate_id']==ident)
        for i,label in enumerate([r'Re $\epsilon_r$',r'Im $\epsilon_r$',r'Phase departure $\rho_\perp$']):
            ax=axs[i,j];im=ax.imshow(layers[i][j],origin='lower',cmap='RdBu_r' if i==2 else 'viridis',vmin=limits[i][0],vmax=limits[i][1],interpolation='nearest')
            if j==0:ax.set_ylabel(label+'\ny cell')
            if i==2:ax.set_xlabel('x cell')
            ax.set_xticks([0,12,23]);ax.set_yticks([0,12,23])
            if i==0:ax.set_title(names[j]+f"\nell={float(row['ell']):.3f}; rW={float(row['whitened_residual']):.3f}")
        data['E0'].append({'name':names[j],'path':str(p.relative_to(ROOT.parent)),'audit':row,'arrays':{'epsilon_real':layers[0][j].tolist(),'epsilon_imag':layers[1][j].tolist(),'rho_orth':layers[2][j].tolist()}})
    for i in range(3):fig.colorbar(axs[i,0].images[0],ax=axs[i,:],shrink=.78,label=['relative permittivity','relative permittivity','signed phase departure'][i])
    fig.supxlabel('All three fail necessary original-prior rules. Shared sigma=0.000664225; tau=23.785. N3 fails noise ball.\nShared row color scales; ell is a point likelihood cost, not posterior mass.',fontsize=10,color='#8b1d1d')
    for ext in ['png','pdf']:fig.savefig(FIG/f'e0_same_observation.{ext}',dpi=170)
    plt.close(fig)
    data['E0_color_limits']=limits
    pair_chis=[(load_chi(ROOT/p['anchor_array']),load_chi(ROOT/p['candidate_array'])) for p in pairs]
    shared_real=(min((1+x.real).min() for ab in pair_chis for x in ab),max((1+x.real).max() for ab in pair_chis for x in ab))
    shared_imag=(min(x.imag.min() for ab in pair_chis for x in ab),max(x.imag.max() for ab in pair_chis for x in ab))
    for p,(a,b) in zip(pairs,pair_chis):
        eligible=[r for r in search if r['family']==p['family'] and r['structure_qualified']=='True' and r['prior_status']=='P0_EXACT']
        assert p['pair_id']==min(eligible,key=lambda r:float(r['d_pred']))['pair_id']
        counts={e:sum(r['competitive_primary']=='True' for r in replay if r['pair_id']==p['pair_id'] and r['truth_endpoint']==e) for e in ['a','b']}
        totals={e:sum(r['pair_id']==p['pair_id'] and r['truth_endpoint']==e for r in replay) for e in ['a','b']}
        vals=[[1+a.real,1+b.real,b.real-a.real],[a.imag,b.imag,b.imag-a.imag],[(np.conj(V)*(a-bg(p['background']))).real/1.2,(np.conj(V)*(b-bg(p['background']))).real/1.2,(np.conj(V)*(b-a)).real/1.2]]
        fig,axs=plt.subplots(3,3,figsize=(12.1,10.4),layout='constrained')
        fig.suptitle(f"{p['family']} / {p['background']} — nearest structurally qualified P0_EXACT pair\nd_pred={p['d_pred']:.3f} (required <=2); same fixed sigma={p['sigma_fixed']:.7g}\n32-noise competition: A truth {counts['a']}/{totals['a']}, B truth {counts['b']}/{totals['b']} — NO COMPETING WITNESS",fontsize=14,color='#8b1d1d')
        for i,label in enumerate([r'Re $\epsilon_r$',r'Im $\epsilon_r$',r'Increment $\rho/1.2$']):
            lim=[shared_real,shared_imag,(0,1)][i]
            dm=max(abs(vals[i][2]).max(),1e-15)
            for j in range(3):
                ax=axs[i,j]; im=ax.imshow(vals[i][j],origin='lower',interpolation='nearest',cmap='RdBu_r' if j==2 else 'viridis',vmin=-dm if j==2 else lim[0],vmax=dm if j==2 else lim[1])
                ax.set_xticks([0,12,23]);ax.set_yticks([0,12,23])
                if j==0:ax.set_ylabel(label+'\ny cell')
                if i==2:ax.set_xlabel('x cell')
                if i==0:ax.set_title(['A: fixed anchor','B: closest qualified endpoint','B - A (signed)'][j])
            fig.colorbar(axs[i,0].images[0],ax=axs[i,:2],shrink=.77)
            fig.colorbar(axs[i,2].images[0],ax=axs[i,2],shrink=.77)
        fig.supxlabel(f"Shape RMS={p['shape_multiscale_rms']:.3f}; occupancy={p['occupancy_multithreshold_mean']:.3f}; boundary={p['boundary_multithreshold_mean']:.3f}; persistent levels={p['persistent_levels']}/5\nFinite bounded search: no global uniqueness claim; no posterior-mass estimate. Shared A/B scales across families.",fontsize=10)
        stem='e1_'+p['family']
        for ext in ['png','pdf']:fig.savefig(FIG/f'{stem}.{ext}',dpi=170)
        plt.close(fig)
        data['E1'].append({'pair':p,'filename':stem,'competition_counts':counts,'competition_totals':totals,'arrays':{'epsilon_real_A':vals[0][0].tolist(),'epsilon_real_B':vals[0][1].tolist(),'epsilon_real_difference':vals[0][2].tolist(),'epsilon_imag_A':vals[1][0].tolist(),'epsilon_imag_B':vals[1][1].tolist(),'epsilon_imag_difference':vals[1][2].tolist(),'rho_normalized_A':vals[2][0].tolist(),'rho_normalized_B':vals[2][1].tolist(),'rho_normalized_difference':vals[2][2].tolist()}})
    data['E1_shared_color_limits']={'epsilon_real':shared_real,'epsilon_imag':shared_imag,'rho_normalized':[0,1]}
    data['decision']=decision
    data['postprocess_wall_seconds']=time.perf_counter()-started
    (OUT/'GALLERY_DATA.json').write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')
    cards=''.join(f'<section><h2>{p["family"]}：原先验内的最近结构合格材料对</h2><p>d_pred={p["d_pred"]:.3f}，远高于冻结门槛2；两端各32次固定噪声重放，竞争支持均为0/32。它是有限搜索中的负结果，不是多解成功。</p><a href="figures/e1_{p["family"]}.pdf">PDF</a><img src="figures/e1_{p["family"]}.png" alt="{p["family"]} pair"></section>' for p in pairs)
    (OUT/'gallery.html').write_text(f'''<!doctype html><html lang="zh"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>原先验支持：科学图库</title><style>body{{font:17px/1.7 system-ui,sans-serif;max-width:1200px;margin:30px auto;padding:0 20px;color:#172234;background:#f5f7fa}}h1{{line-height:1.3}}section{{background:white;padding:24px;margin:24px 0;border-radius:10px}}img{{width:100%;height:auto}}.negative{{border-left:6px solid #a32727;background:#fff1f1;padding:18px}}a{{color:#1b5f9e}}</style><h1>原先验支持与竞争材料见证：科学图库</h1><div class="negative"><b>NO_GO_FOR_GENERATOR_COMPARISON</b><br>E1未发现合格竞争见证：4族中结构合格的最近材料对，预测间距均超过门槛2，两端噪声重放竞争均0/32。E2未激活。有限搜索没有证明全局唯一性，歧义仍未解决。</div><p>本图库只读取现存材料数组、已缓存精确预测审计和已完成E1指标。新增物理状态/RHS/分解/训练均为0；不读取truth指导选择或搜索。点似然不是posterior质量。</p><section><h2>同一历史观测：强Top-1与数据可行的补偿材料</h2><p>p250065_ring：展示原N3、minchange_base及T v1 edit00补偿端点。N3保留历史强Top-1身份，但本噪声球检验失败。两种T材料在历史噪声球内可行，却被原固定增量相位/范围必要条件排除；N3同样不属于原先验的必要支持规则。色标逐行固定，εᵣ=1+χ；ρ⊥=Im[conj(v)(χ−χ背景)]，v=(1+0.2i)/√1.04。</p><p>本图只画材料与相位偏离；ell数值直接复用cached_audit，未重新计算物理。绝对空气孔洞不能代替原增量先验的结构定义。</p><a href="figures/e0_same_observation.pdf">PDF</a><img src="figures/e0_same_observation.png" alt="same observation material audit"></section>{cards}<section><h2>可复核数据与边界</h2><p><a href="GALLERY_DATA.json">全部绘图源数据JSON</a> · <a href="PROTOCOL_FROZEN.json">冻结协议</a> · <a href="E1_DECISION.json">E1裁决</a> · <a href="WITNESS_PAIRS.json">候选材料对</a></p><p>图中A/B材料色标跨family固定；差分色标为各材料对对称色标。网格为24×24原数组，显示cell索引，不插值。每个family选取E1_SEARCH中P0_EXACT且结构合格的最小d_pred，未依据truth或误差挑选。噪声支持频率为同时满足两端噪声球及|Δell|≤log10的频率，不是posterior质量。</p><p>实际绘图后处理耗时 {data['postprocess_wall_seconds']:.3f} 秒；物理调用0。E1本身预算与失败账本请见主报告，本图库的0新增物理不覆盖已经完成的E1运行。</p></section></html>''')
    print(json.dumps({'postprocess_wall_seconds':data['postprocess_wall_seconds'],'new_physics_states':0,'new_rhs':0,'figures':5,'source_json':str(OUT/'GALLERY_DATA.json')}))

if __name__=='__main__':main()
