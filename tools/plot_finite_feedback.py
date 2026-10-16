#!/usr/bin/env python3
"""Plot paid saved arrays only. No solver imports, solves, or ledger writes."""
from pathlib import Path
import json
import html
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / 'research/finite_feedback_fiber_diagnostic_20261016'
OUT = RUN / 'figures'
OUT.mkdir(exist_ok=True)
M = json.loads((RUN / 'METRICS.json').read_text())
LEDGER = [json.loads(s) for s in (RUN / 'RESOURCE_LEDGER.jsonl').read_text().splitlines()]
plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': 'white', 'savefig.facecolor': 'white', 'pdf.fonttype': 42})
FIGS = []

def save(fig, name, scope):
    for ext in ['png', 'pdf']:
        fig.savefig(OUT / f'{name}.{ext}', dpi=190, bbox_inches='tight')
    plt.close(fig)
    FIGS.append((name, scope))

def arr(name, key, endpoints=False):
    path = RUN / 'results' / ('endpoints' if endpoints else '') / (name + '.npz')
    with np.load(path, allow_pickle=False) as z:
        return z[key].copy()

# Oracle spectrum and explicit cutoff dependence; rectangular null is retained.
fig, ax = plt.subplots(1, 3, figsize=(16, 4.4), layout='constrained')
for a, color in [('physics_R_base', '#176b9a'), ('N3', '#d95c35')]:
    d = next(d for d in M['decomposition'] if d['anchor'] == a)
    s = arr(a + '_DECOMPOSITION', 'singular_values')
    c = arr(a + '_DECOMPOSITION', 'GT_coefficients')
    h = arr(a + '_DECOMPOSITION', 'GT_mass_displacement')
    total = float(np.sum(h*h))
    ax[0].semilogy(np.arange(1, len(s)+1), s/s[0], label=a, color=color)
    # Weak energy at each singular-value threshold includes unrepresented nullspace.
    weak = np.clip(1 - np.cumsum(c*c)/total, 0, 1)
    ax[1].plot(np.arange(1, len(s)+1), weak, color=color, label=a)
    th = d['spectral_thresholds']
    ax[2].semilogx([x['relative_cutoff'] for x in th],
                  [x['weak_including_rectangular_null_fraction'] for x in th], 'o-', color=color, label=a)
    for x in th:
        ax[2].annotate(f"r={x['strong_rank']}", (x['relative_cutoff'], x['weak_including_rectangular_null_fraction']),
                       xytext=(4, 6), textcoords='offset points', fontsize=9)
ax[0].set(xlabel='Singular index', ylabel='Singular value / largest', title='Whitened mass Jacobian spectrum')
ax[1].set(xlabel='Retained leading singular directions', ylabel='Remaining GT displacement energy fraction', title='GT energy outside retained span', ylim=(0,1.04))
ax[2].set(xlabel='Relative singular-value cutoff', ylabel='Weak GT energy fraction', title='Cutoff dependence (nullspace included)', ylim=(.65,1.04))
for x in ax: x.grid(alpha=.2); x.legend(fontsize=9)
fig.suptitle('GT geometry diagnostic: two anchors | oracle displacement, off-prior straight paths', fontsize=14)
save(fig, '01_gt_spectrum', 'GT 仅用于几何诊断；弱能量包含矩形零空间，不代表已找到等数据解。')

fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True, layout='constrained')
for a, ax in zip(['physics_R_base', 'N3'], axes):
    rows = [r for r in M['feedback'] if r['endpoint'].startswith(a+'_oracle_t')]
    t = [float(r['endpoint'].split('_t')[-1]) for r in rows]
    for order in [1,2,3,6,20]:
        ys = [r['relative_prediction_errors'].get(f'order_{order}', np.nan) for r in rows]
        ax.semilogy(t, ys, 'o-' if order != 20 else 'D', label=f'Order {order}')
    ax.semilogy(t, [r['rational_relative_error'] for r in rows], 'k--s', label='Rational identity')
    ax.set(title=a, xlabel='t on anchor-to-GT path', xticks=t)
    ax.grid(alpha=.2); ax.legend(fontsize=9, ncol=2)
axes[0].set_ylabel('Relative error in finite data increment')
fig.suptitle('Saved finite feedback predictions | order 20 recorded at t=1 only', fontsize=14)
save(fig, '02_gt_feedback', '两条使用 GT 的离先验路径；高阶与有理表达的误差核对，不是逆问题发现。')

methods = [('frozen_weak', 'Frozen direction', '#176b9a'),
           ('second_order_curvature', 'Second-order candidate', '#d95c35'),
           ('updated_J_predictor', 'Updated-J predictor', '#a087be'),
           ('updated_J_corrector', 'Updated-J corrector', '#368b70')]
known = [p for p in M['paths'] if 'rows' in p]
fig, axes = plt.subplots(2, 2, figsize=(13, 9.7), layout='constrained')
for col, p in enumerate(known):
    for method, label, color in methods:
        rows = [r for r in p['rows'] if r['method']==method]
        for row, key in enumerate(['data_displacement', 'whitened_residual']):
            ax = axes[row,col]
            ax.plot([r['material_RMS'] for r in rows], [r[key] for r in rows], marker='o',
                    color=color, label=label, linestyle='--' if method.endswith('predictor') else '-')
            for r in rows:
                ax.annotate(f"{r['t']:g}", (r['material_RMS'],r[key]), xytext=(3,4), textcoords='offset points', fontsize=8, color=color)
    axes[0,col].set_title(p['anchor']+' | known TRAIN calibration')
    axes[1,col].axhline(23.785449, color='black', linestyle=':', label='Unchanged tau = 23.785449')
    for ax in axes[:,col]: ax.grid(alpha=.2); ax.set_xlabel('Material RMS displacement from known TRAIN base')
    axes[0,col].legend(fontsize=8)
    axes[1,col].legend(fontsize=8)
axes[0,0].set_ylabel('Whitened prediction displacement from base')
axes[1,0].set_ylabel('Whitened residual to original observation')
fig.suptitle('Known TRAIN paths | unequal paid work; point labels are t, not equal cost', fontsize=14)
fig.supxlabel('Paid per anchor: base 1 state / 8 RHS + initial J 32 RHS; frozen 4 states / 256 RHS incl. feedback audits;\nsecond-order 3 states / 24 RHS + curvature 8 RHS; updated-J 6 states / 48 RHS + two refreshed J / 64 RHS.\nPredictors and correctors are sequential; shared setup is not allocated as an equal per-point cost.', fontsize=10)
save(fig, '03_known_train_paths', '已知 TRAIN 校准，不是逆问题发现或多模态样本；原始阈值不变，路径成本不相等。')

for p in known:
    a = p['anchor']
    ids = [a+'_base', a+'_weak_t0.75', a+'_curvature_t0.75', a+'_pc_3_correct']
    fields = [arr(x, 'chi', True).reshape(24,24) for x in ids]
    titles = ['Known TRAIN base', 'Frozen direction t=0.75', 'Second-order t=0.75', 'Updated-J step 3 corrector']
    fig, axes = plt.subplots(4,4,figsize=(14,12), layout='constrained')
    for row in range(4):
        values = ([1+f.real for f in fields] if row==0 else [f.imag for f in fields] if row==1
                  else [f.real-fields[0].real for f in fields] if row==2 else [f.imag-fields[0].imag for f in fields])
        if row<2: lo,hi=min(v.min() for v in values),max(v.max() for v in values); cmap='viridis'
        else: hi=max(max(abs(v.min()),abs(v.max())) for v in values); lo=-hi; cmap='RdBu_r'
        for col,v in enumerate(values):
            im=axes[row,col].imshow(v, origin='lower', cmap=cmap, vmin=lo, vmax=hi, interpolation='nearest', extent=[-.5,.5,-.5,.5])
            axes[row,col].set_xticks([]); axes[row,col].set_yticks([])
            if row==0: axes[row,col].set_title(titles[col],fontsize=10)
        axes[row,0].set_ylabel(['Re(1 + chi)','Im(chi)','Delta Re(chi)','Delta Im(chi)'][row])
        fig.colorbar(im, ax=list(axes[row,:]), shrink=.82, label='Shared color limits within this row')
    fig.suptitle(a+' | known TRAIN calibration materials; differences are relative to base',fontsize=14)
    save(fig, '04_material_'+a, '已知 TRAIN 材料与相对基准差值；每行统一色标，不将微小变化放大成新模式。')

# Same exposed public case: truth is an offline reference, never an online choice.
with np.load(RUN/'inputs/physics_R_base.npz',allow_pickle=False) as z:
    phys=z['chi'].copy().reshape(24,24); truth=z['truth'].copy().reshape(24,24)
with np.load(RUN/'inputs/N3.npz',allow_pickle=False) as z:
    neural=z['chi'].copy().reshape(24,24)
fields=[phys,neural,truth]
decomp={d['anchor']:d for d in M['decomposition']}
gt_residual=next(p for p in M['paths'] if p['anchor']=='physics_R_base')['points'][-1]['whitened_residual']
stats=[(decomp[a]['GT_material_MAE'],decomp[a]['anchor_whitened_residual']) for a in ['physics_R_base','N3']]+[(0.,gt_residual)]
fig,axes=plt.subplots(4,3,figsize=(12,13),layout='constrained')
for row in range(4):
    values=([1+f.real for f in fields] if row==0 else [f.imag for f in fields] if row==1
            else [f.real-truth.real for f in fields] if row==2 else [f.imag-truth.imag for f in fields])
    if row<2:lo,hi=min(v.min() for v in values),max(v.max() for v in values);cmap='viridis'
    else:hi=max(np.max(abs(v)) for v in values);lo=-hi;cmap='RdBu_r'
    for col,v in enumerate(values):
        im=axes[row,col].imshow(v,origin='lower',interpolation='nearest',cmap=cmap,vmin=lo,vmax=hi,extent=[-.5,.5,-.5,.5])
        axes[row,col].set_xticks([]);axes[row,col].set_yticks([])
        if row==0:
            mae,res=stats[col]
            axes[row,col].set_title(['physics_R_base','N3','Exposed GT (offline only)'][col]+f'\nMAE {mae:.6f} | residual {res:.3f}',fontsize=11)
    axes[row,0].set_ylabel(['Re(1 + chi)','Im(chi)','Re(chi - GT)','Im(chi - GT)'][row])
    fig.colorbar(im,ax=list(axes[row,:]),shrink=.8,label='Shared color limits within this row')
fig.suptitle('Same public case | OFFLINE material diagnostic',fontsize=16)
fig.supxlabel('Fixed tau = 23.785449. No independently converged physics-only reconstruction is shown.\nExposed GT is an offline comparison reference; no algorithm selects truth.',fontsize=11)
save(fig,'07_publiccase_materials','同一公开案例的离线材料对照；GT 已暴露，仅用于事后比较，不是算法选出的结果；没有独立收敛的纯物理重建。')

post = json.loads((RUN/'POSTHOC_CONSISTENCY.json').read_text())
fig, axes = plt.subplots(1, 3, figsize=(17, 5.5), layout='constrained', gridspec_kw={'width_ratios':[1,1,.85]})
for ax, a in zip(axes[:2], ['250033_weak','250034_ring']):
    rows = sorted([r for r in post['chart_vs_Maxwell'] if r['anchor']==a], key=lambda r:r['t'])
    for key,label,color,style in [
        ('frozen_chart_response',r'Frozen chart: $\|t J Dg\,v\|$','#176b9a','o--'),
        ('actual_material_increment_J_response',r'Material secant: $\|J\,[g(u+tv)-g(u)]\|$','#d95c35','o-'),
        ('actual_full_response',r'Full response: $\|F(g(u+tv))-F(g(u))\|$','#368b70','x:')]:
        ax.semilogy([r['t'] for r in rows],[r[key] for r in rows],style,color=color,label=label,markersize=7)
    ax.set(title=a+' | known TRAIN', xlabel='Prior-coordinate step t', ylabel='Whitened response norm',xticks=[.25,.5,.75,1])
    ax.grid(alpha=.2); ax.legend(fontsize=9,loc='upper left')
    r=rows[-1]
    ax.text(.04,.17,f"t=1: frozen chart {r['frozen_chart_response']:.5g}\nmaterial secant {r['actual_material_increment_J_response']:.5g}; full {r['actual_full_response']:.5g}\nMaxwell-only relative linearization error: {r['Maxwell_only_linearization_error']:.5g}",transform=ax.transAxes,fontsize=9,bbox={'facecolor':'white','edgecolor':'#ccd6df','alpha':.9})
g=M['geometry']; x=np.arange(len(g)); width=.34
for offset,key,label,color in [(-width/2,'Maxwell_curvature_norm','Maxwell curvature','#176b9a'),(width/2,'decoder_curvature_norm','Decoder curvature','#d95c35')]:
    vals=[r[key] for r in g]
    axes[2].bar(x+offset,vals,width,label=label,color=color)
    for xx,yy in zip(x+offset,vals): axes[2].annotate(f'{yy:.3g}',(xx,yy),xytext=(0,5),textcoords='offset points',ha='center',fontsize=9)
axes[2].set(yscale='log',xticks=x,xticklabels=['250033_weak','250034_ring'],ylabel='Saved directional response-curvature norm',title='Local curvature components',ylim=(4e-5,100))
axes[2].legend(fontsize=9); axes[2].grid(axis='y',alpha=.2)
fig.suptitle('Chart curvature versus Maxwell feedback | posthoc decomposition of paid arrays',fontsize=15)
fig.supxlabel('Response norms do not establish vector agreement; the reported Maxwell-only error compares vectors.\nKnown TRAIN calibration only; these curves do not demonstrate Maxwell divergence or inverse discoveries.',fontsize=10)
save(fig,'06_chart_vs_maxwell','已支付数组的事后拆分：比较先验切线、实际材料增量与完整响应；区分 decoder 曲率与 Maxwell 曲率，不声称 Maxwell 发散。')

# Direct attribution by endpoint; shared Jacobian/curvature costs stay explicit.
cost_rows=[]
for entry in json.loads((RUN/'results/ENDPOINTS.json').read_text()):
    name=entry['id']
    events=[e for e in LEDGER if e.get('context',{}).get('endpoint')==name]
    cost_rows.append({'endpoint':name,'states':sum(e['amount'] for e in events if e['event']=='STATE'),
                      'direct_RHS':sum(e['amount'] for e in events if e['event']=='RHS')})
(OUT/'endpoint_costs.json').write_text(json.dumps(cost_rows,indent=2))
cost_table='<section><h2>已支付端点的直接归属成本</h2><p>下面逐端点统计账本中的直接归属；共享初始 J、刷新 J、曲率项与路径导数另计，不是从零生成端点的独立算法成本。完整实验总计仍为 45 states / 1880 RHS。预测与校正按顺序执行，不能作为等成本对照。</p><table><tr><th>Endpoint</th><th>States</th><th>Direct RHS</th></tr>'+''.join(f"<tr><td>{html.escape(r['endpoint'])}</td><td>{r['states']}</td><td>{r['direct_RHS']}</td></tr>" for r in cost_rows)+'</table></section>'
cards='\n'.join(f'<section><h2>{html.escape(n)}</h2><p>{html.escape(s)}</p><a href="figures/{n}.pdf">PDF</a> · <a href="figures/{n}.png">PNG</a><img src="figures/{n}.png" alt="{html.escape(n)}"></section>' for n,s in FIGS)
(RUN/'gallery.html').write_text('''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Finite feedback diagnostic gallery</title><style>body{font:17px/1.65 system-ui;max-width:1250px;margin:40px auto;padding:0 22px;color:#15293b;background:#f6f8fb}h1{line-height:1.2}section{background:white;padding:24px;margin:25px 0;border-radius:12px}img{display:block;width:100%;height:auto;margin-top:14px}a{color:#126292}.scope{padding:20px;background:#e5edf5;border-radius:10px}</style><h1>有限反馈诊断图集</h1><p class="scope">仅重绘已支付的保存结果：45 states / 1880 RHS；本次绘图新增物理成本为零。GT 路径仅作几何诊断；已知 TRAIN 路径仅作校准，不是逆问题发现、解纤维证明、后验样本或多模态发现。原始噪声阈值 tau=23.785449 保持不变。所有图标题为英文，图片和 PDF 可独立下载。</p>'''+cards+cost_table+'</html>',encoding='utf-8')
print(json.dumps({'gallery':str(RUN/'gallery.html'),'figures':[n for n,_ in FIGS], 'new_physics_states':0,'new_RHS':0},indent=2))
