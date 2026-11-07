"""Scientific plots with explicit sample/class/style training exposure labels."""
from pathlib import Path
import html,json,textwrap
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

def load_npz(path):
    with np.load(path,allow_pickle=False)as z:return {k:z[k]for k in z.files}

def exposure(row):
    known=lambda k:'SEEN' if row.get(k) else 'UNSEEN'
    label=f"digit {row['digit']}" if row.get('digit')is not None else row['family']
    font=known('seen_font')if row.get('font')else'N/A'
    return (f"{row['split'].upper()} {row['parent_id']} | image {'TRAINED' if row.get('seen_image') else 'NOT TRAINED'} | "
            f"{label}: {known('seen_class')} | font {font} / rule {known('seen_rule')} | "
            f"A={row['amplitude']:g} | {'clean' if row['noise']==0 else 'noisy'}")

def short(run):
    names={'strict_gn45_anchor':'Strict GN','continue_trf':'Continue TRF',
           'weakreg45_anchor':'Weak-reg TRF','discrepancy45_anchor':'TRF final-stop',
           'discrepancy_allscales45_anchor':'TRF scale-stop','cie_dual_cold55':'CIE / Fourier',
           'recurrent_refreshed_q_seed7_cascade2':'Q-refresh RNN',
           'recurrent_feedback1_seed7_cascade2':'Feedback RNN',
           'recurrent_feedback_hinge_seed7_cascade2':'Hinge RNN',
           'g_weak_seed7_cascade2_const_ft_material':'G-weak + material FT',
           'g_weak_seed7_cascade2_const_ft_fullvie':'G-weak + VIE FT'}
    if run in names:return names[run]
    s=run.replace('_seed7_cascade2','').replace('_seed17_cascade2',' s17').replace('_seed27_cascade2',' s27')
    return s.replace('v_physics','Voxel').replace('j_full','J-full').replace('j_null','J-null').replace('g_weak','G-weak').replace('g_dark','G-dark').replace('soft_dark','Soft-dark').replace('continue_trf','Continue TRF')

def metric_text(row,stage):
    if not row:return 'NOT RUN'
    c=row.get('top1') if stage=='final' else next((c for c in row.get('candidates',[])if c['label']=='raw_1'),None)
    if c is None:return 'NO FEASIBLE TOP1'if stage=='final'else'NO RAW MATERIAL'
    return f"N={c['NRMSE']:.3f}; I={c['IoU_Otsu']:.2f}\nr={c['relative_residual']:.1e}"

def make_plots(scores_path,study,runs=None,tag=None,gallery=True):
    study=Path(study);scores=json.loads(Path(scores_path).read_text());out=study/'figures/final'
    # This renderer selects the predeclared TEST image atlas. VAL selection
    # charts have their own score-only entrypoints and never open this atlas.
    if not (study/'audit/METHOD_FREEZE.json').exists() or not scores['rows'] or any(r['split']!='test'for r in scores['rows']):
        raise ValueError('The image atlas requires frozen TEST scores; use the VAL score-only plotting entrypoint for model selection')
    if tag:out=out/tag
    out.mkdir(parents=True,exist_ok=True)
    records={r['obs_id']:r for r in scores['rows']};by={(r['run'],r['obs_id']):r for r in scores['rows']}
    all_runs=list(dict.fromkeys(r['run']for r in scores['rows']))
    runs=runs or [r for r in all_runs if 'seed17'not in r and 'seed27'not in r and 'posthoc'not in r]
    selection=json.loads((study/'audit/FIGURE_SELECTION.json').read_text())['rows']
    meta=json.loads((study/'data/manifests/evaluation.json').read_text())['records']
    test=[r for r in meta if r['split']=='test'];by_parent_noise={(r['parent_id'],r['noise']):r for r in test}
    mats={};truths={};anchors={}
    def arrays(run,obs):
        key=(run,obs)
        if key not in mats:
            path=study/'evaluation'/run/f'{obs}.npz';mats[key]=load_npz(path)if path.exists()else{}
        return mats[key]
    def truth(r):
        if r['parent_id']not in truths:truths[r['parent_id']]=load_npz(study/'data'/r['truth_path'])['m'].reshape(24,24)
        return truths[r['parent_id']]
    def anchor(r):
        if r['obs_id']not in anchors:anchors[r['obs_id']]=load_npz(study/'anchors'/f"{r['obs_id']}.npz")['m0'].reshape(24,24)
        return anchors[r['obs_id']]
    def material(run,r,stage):
        a=arrays(run,r['obs_id'])
        if stage=='final':return np.asarray(a['top1']).reshape(24,24)if'top1'in a else None
        row=by.get((run,r['obs_id']),{})
        raw=next((c for c in row.get('candidates',[])if c['label']=='raw_1'),None)
        return np.asarray(a['candidate_materials'][raw['index']]).reshape(24,24)if raw is not None else None
    def panel(ax,image,title,caption='',vmax=.55):
        if image is None:ax.text(.5,.5,'NOT RUN /\nNO FEASIBLE TOP1',ha='center',va='center',fontsize=7,transform=ax.transAxes)
        else:
            ax.imshow(image,cmap='viridis',vmin=.02,vmax=vmax,interpolation='nearest')
            saturated=int(np.count_nonzero((image<.02)|(image>vmax)))
            if saturated:caption+=("\n"if caption else"")+f'color-saturated cells: {saturated}'
        ax.set_xticks([]);ax.set_yticks([]);ax.set_title(title,fontsize=8)
        if caption:ax.set_xlabel(caption,fontsize=7)
    produced=[]
    for noise in (0.,.001):
        chosen=[by_parent_noise[(s['parent_id'],noise)]for s in selection]
        for stage in ('raw','final'):
            fig,axs=plt.subplots(len(chosen),2+len(runs),figsize=(2.0*(2+len(runs)),2.2*len(chosen)),squeeze=False,layout='constrained')
            for i,r in enumerate(chosen):
                panel(axs[i,0],truth(r),'GT');panel(axs[i,1],anchor(r),'Strict GN')
                identity_label=r['digit']if r.get('digit')is not None else r['family']
                axs[i,0].set_ylabel(f"{identity_label}\n{r['group']}",fontsize=8)
                axs[i,0].set_title('GT | '+r['parent_id'],fontsize=8)
                c='SEEN'if r.get('seen_class')else'UNSEEN';f=('SEEN'if r.get('seen_font')else'UNSEEN')if r.get('font')else'N/A'
                rule='SEEN'if r.get('seen_rule')else'UNSEEN'
                axs[i,0].set_xlabel(f"image/mother NOT TRAINED\nclass {c}; font {f}\nrule {rule}; A={r['amplitude']:g}",fontsize=7)
                for j,run in enumerate(runs):panel(axs[i,j+2],material(run,r,stage),short(run),metric_text(by.get((run,r['obs_id'])),stage))
            fig.suptitle(f"NEW INDEPENDENT TEST | {'clean'if noise==0 else'noise=1e-3'} | {stage}\nAll displayed images/mothers absent from TRAIN; class/style exposure explicitly labelled. Fixed material scale [.02,.55].",fontsize=13)
            im=plt.cm.ScalarMappable(norm=plt.Normalize(.02,.55),cmap='viridis');fig.colorbar(im,ax=axs.ravel().tolist(),shrink=.3,label='Real material contrast; values beyond range saturate',extend='both')
            file=out/f"test_{'clean'if noise==0 else'noisy'}_{stage}.png";fig.savefig(file,dpi=150);plt.close(fig);produced.append(str(file.relative_to(study)))
    # Full all-case gallery, method chunks kept legible and every failure retained.
    gallery_dir=study/'figures/all_test_cases';gallery_dir.mkdir(parents=True,exist_ok=True);pages=[]
    gallery_runs=all_runs
    for r in test if gallery else []:
        if r['obs_id']not in records:continue
        for start in range(0,len(gallery_runs),5):
            chunk=gallery_runs[start:start+5];fig,axs=plt.subplots(2,2+len(chunk),figsize=(2.4*(2+len(chunk)),5.8),squeeze=False,layout='constrained')
            for i,stage in enumerate(('raw','final')):
                panel(axs[i,0],truth(r),'GT');panel(axs[i,1],anchor(r),'Strict GN')
                axs[i,0].set_ylabel(stage,fontsize=10)
                for j,run in enumerate(chunk):panel(axs[i,j+2],material(run,r,stage),short(run),metric_text(by.get((run,r['obs_id'])),stage))
            fig.suptitle(textwrap.fill(exposure(r),110)+'\nFixed material color scale [.02,.55]; raw and true-physics-qualified online top1.',fontsize=9)
            # Lossless export keeps every rendered scientific pixel while
            # avoiding a very large static-site atlas. Main figures stay PNG.
            file=gallery_dir/f"{r['obs_id']}_p{start//5+1}.webp"
            fig.set_dpi(130);fig.canvas.draw()
            Image.fromarray(np.asarray(fig.canvas.buffer_rgba())).save(file,format='WEBP',lossless=True,method=3)
            plt.close(fig)
            pages.append({'file':file.name,'identity':exposure(r),'parent_id':r['parent_id'],'obs_id':r['obs_id'],'runs':chunk})
    if gallery:(gallery_dir/'INDEX.json').write_text(json.dumps(pages,indent=2)+'\n')
    body=['<!doctype html><meta charset="utf-8"><title>R4 独立测试完整图册</title><style>body{font-family:system-ui;margin:2rem}img{width:100%;max-width:1500px}figure{margin:2rem 0}figcaption{white-space:pre-wrap}</style><h1>R4 独立测试完整图册</h1><p>所有对象均未参与训练。类别、字体及规则是否见过，逐图标明。包含失败与未完成结果；色阶统一为[0.02,0.55]。</p>']
    body += [f'<figure><figcaption>{html.escape(p["identity"])}</figcaption><img loading="lazy" src="{p["file"]}"></figure>'for p in pages]
    if gallery:(gallery_dir/'GALLERY.html').write_text('\n'.join(body))
    # Per-noise quality/data/cost points; seeds are never pooled as independent scenes.
    for noise in (0.,.001):
        fig,axs=plt.subplots(1,3,figsize=(15,4.5),layout='constrained')
        for run in all_runs:
            rs=[r for r in scores['rows']if r['run']==run and r['noise']==noise and r.get('top1')]
            if not rs:continue
            n=np.mean([r['top1']['NRMSE']for r in rs]);res=np.mean([r['top1']['relative_residual']for r in rs]);wall=np.mean([r['wall_seconds']for r in rs]);rhs=np.mean([r['costcounts'].get('actual_complex_RHS',0)for r in rs])
            planned=sum(r['run']==run and r['noise']==noise for r in scores['rows'])
            axs[0].scatter(res,n,label=short(run)+f' ({len(rs)}/{planned} qualified)');axs[1].scatter(wall,n);axs[2].scatter(rhs,n)
        axs[0].set(xlabel='Full-VIE relative data residual (mean)',ylabel='Material NRMSE (mean)',xscale='log');axs[1].set(xlabel='Reported method seconds (scope in ledger)',ylabel='Material NRMSE');axs[2].set(xlabel='Charged method complex RHS (mean)',ylabel='Material NRMSE')
        axs[0].legend(fontsize=6);fig.suptitle(f"NEW TEST | {'clean'if noise==0 else'noisy'} | conditional on qualified top1 (counts shown); training and shared setup separately invoiced")
        file=out/f"quality_data_cost_{'clean'if noise==0 else'noisy'}.png";fig.savefig(file,dpi=160);plt.close(fig);produced.append(str(file.relative_to(study)))
    (out/'FIGURES.json').write_text(json.dumps({'files':produced,'all_case_pages':len(pages),'atlas_format':'lossless WebP; main scientific figures PNG','training_image_claim':'all NEW TEST objects and mothers excluded from TRAIN, identity flags from frozen manifest','selection':'pre-registered first case per family/digit at A=.2; no GT-outcome selection'},indent=2)+'\n')
    return out
