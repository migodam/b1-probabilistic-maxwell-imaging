"""Post-seal display only. Imports no forward model or optimizer.

Run with Inv_SLAM/Gaussian/.venv_nn/bin/python; all outputs beside this file.
Original scored metrics, arrays, images and online selections remain untouched.
"""
import time
CPU_START = time.process_time()
WALL_START = time.monotonic()
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
CONF = ROOT / 'delegated/confirmation'
PUBLIC = CONF / 'selected_public_outputs'
METHODS = ('cie_single','cie_dual','global_dual','a_plus','fullG19','dr','multistart','nn_direct')
LABELS = dict(zip(METHODS, ('CIE','Dual','Global','A+','G19','DR','Multi','NN')))
METRICS_PATH = CONF / 'RESULTS/offline_scoring_v1/per_object.json'
MANIFEST_PATH = PUBLIC / 'COMPACT_PUBLIC_MANIFEST.json'
META_PATH = CONF / 'sealed_truth/confirm24/generator_metadata.json'
metrics = {(r['case_id'],r['method']):r for r in json.loads(METRICS_PATH.read_text())
           if r['stage']=='confirm24' and r['selector']=='primary' and r['view']=='endpoint'}
manifest = json.loads(MANIFEST_PATH.read_text())['rows']
meta = {r['case_id']:r for r in json.loads(META_PATH.read_text())['records']}
paths = {(r['case_id'],r['method']):PUBLIC/str(r['compact_primary']).replace('\\','/')
         for r in manifest if r['case_id'] in meta}
arrays = {}
sources = []
for case in sorted(meta):
    truth_path = CONF/'sealed_truth/confirm24'/(case+'.npz')
    with np.load(truth_path,allow_pickle=False) as z: arrays[(case,'GT')] = z['m'].copy()
    sources.append({'case_id':case,'method':'GT','path':str(truth_path.relative_to(ROOT)),'key':'m'})
    for method in METHODS:
        path = paths[(case,method)]
        with np.load(path,allow_pickle=False) as z: arrays[(case,method)] = z['m'].copy()
        assert arrays[(case,method)].shape == (576,) and np.isfinite(arrays[(case,method)]).all()
        sources.append({'case_id':case,'method':method,'path':str(path.relative_to(ROOT)),'key':'m'})

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none'})
outputs = []
def save(fig,name):
    for ext in ('png','pdf'):
        p = OUT/(name+'.'+ext)
        if p.exists(): raise FileExistsError(p)
        fig.savefig(p,dpi=170,bbox_inches='tight',facecolor='white')
        outputs.append(str(p.relative_to(ROOT)))
    plt.close(fig)

def caption(case,method):
    if method == 'GT':
        m=meta[case]
        return case.rsplit('_',1)[-1]+' '+m['family']+'\n'+m['regime']+' '+str(int(m['high_chi']))
    r=metrics[(case,method)]
    return f"N={r['NRMSE']:.3f}  I={r['IoU']:.2f}"

for family in ('levelsets','branches','mixed_islands','continuous_random'):
    cases = [c for c in sorted(meta) if meta[c]['family']==family]
    fig,axes=plt.subplots(len(cases),9,figsize=(14.5,1.52*len(cases)+.7),squeeze=False)
    fig.subplots_adjust(left=.06,right=.965,bottom=.025,top=.95,wspace=.055,hspace=.35)
    for i,case in enumerate(cases):
        for j,method in enumerate(('GT',)+METHODS):
            ax=axes[i,j];ax.set_xticks([]);ax.set_yticks([])
            ax.imshow(arrays[(case,method)].reshape(24,24),vmin=0,vmax=20,cmap='viridis',origin='lower',interpolation='nearest')
            ax.set_xlabel(caption(case,method),fontsize=6.3 if method=='GT' else 7,color='#202020')
            if i==0:ax.set_title('GT' if method=='GT' else LABELS[method],fontsize=10,color='#202020')
    fig.suptitle('confirm24 / primary / structural review; N=amplitude NRMSE, I=chi>1 IoU; neutral labels',fontsize=11,y=.993)
    bar=fig.add_axes([.976,.2,.008,.58])
    fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,20),cmap='viridis'),cax=bar,label='contrast chi')
    save(fig,'confirm24_primary_'+family+'_neutral')

case = 'r12_confirm24_0009'
for windows,name in [((20,),case+'_primary_detail_full0_20'),((20,3),case+'_primary_detail_full_and_weak')]:
    fig=plt.figure(figsize=(10.6*len(windows),10.4))
    subs=fig.subfigures(1,len(windows),squeeze=False).ravel()
    for sub,vmax in zip(subs,windows):
        axes=sub.subplots(3,3)
        sub.subplots_adjust(left=.045,right=.89,bottom=.06,top=.91,wspace=.15,hspace=.29)
        for ax,method in zip(axes.flat,('GT',)+METHODS):
            ax.imshow(arrays[(case,method)].reshape(24,24),vmin=0,vmax=vmax,cmap='viridis',origin='lower',interpolation='nearest')
            ax.set_xticks([]);ax.set_yticks([])
            ax.set_title('GT' if method=='GT' else LABELS[method],fontsize=13,color='#202020')
            ax.set_xlabel(caption(case,method) if method!='GT' else '0009 branches / mixed chi 18',fontsize=10,color='#202020')
        title='FULL contrast scale [0,20]' if vmax==20 else 'WEAK-STRUCTURE display [0,3]\nvalues above 3 saturated for display only'
        sub.suptitle(title,fontsize=16,y=.985)
        bar=sub.add_axes([.92,.14,.016,.69])
        sub.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,vmax),cmap='viridis'),cax=bar,label='contrast chi')
    fig.supxlabel('Same saved primary arrays and original metrics; no optimization, no new success classification.',fontsize=12,y=.007)
    save(fig,name)

receipt = {'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
    'purpose':'Post-hoc structural display requested by user; old threshold not met does not mean no structural recovery.',
    'new_FW':0,'new_optimizer_runs':0,'online_GT_selection_changes':0,'new_success_classifications':0,
    'original_scores_and_images_modified':False,'primary_success_statistics_modified':False,
    'metrics_source':str(METRICS_PATH.relative_to(ROOT)),
    'compact_manifest':str(MANIFEST_PATH.relative_to(ROOT)),
    'metadata_source':str(META_PATH.relative_to(ROOT)),
    'array_sources':sources,'outputs':outputs,
    'display':{'primary_family_scale':[0,20],'colormap':'viridis','origin':'lower','interpolation':'nearest',
        'labels':'All captions neutral #202020; same NRMSE/IoU values from original scored rows.',
        'detail_case':case,'weak_window':[0,3],
        'weak_window_warning':'Post-hoc fixed display saturation only, not a new metric or algorithm. Full [0,20] alongside. New window does not change arrays.'},
    'plotting_CPU_seconds':time.process_time()-CPU_START,'wall_seconds':time.monotonic()-WALL_START}
(OUT/'RECEIPT.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ['new_FW','new_optimizer_runs','plotting_CPU_seconds','wall_seconds']},indent=2),flush=True)
