"""Explicit posthoc gallery after all32 seal and completed approved scoring.
This separate helper may read only the two declared private truth materials.
It never imports a solver, changes predictions, or selects online endpoints.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parents[1]

def read(p):return json.loads(Path(p).read_text())
def material(path):
    with np.load(path,allow_pickle=False) as z:return np.asarray(z['material'],np.complex128).copy()
def rms(a,b):return float(np.linalg.norm(a-b)/24)
def run(posthoc_approved):
    if not posthoc_approved:raise RuntimeError('Explicit root posthoc-gallery authorization required')
    proto=read(P/'config/possibility_volume_protocol.json')
    stats=read(P/'config/possibility_statistics_frozen.json')
    root=Path(proto['output']);score=root/'offline_statistics';seal=read(root/'PREDICTIONS_SEALED.json')
    receipt=read(score/'SCORING_COMPLETE.json');approval=read(P/'config/possibility_volume_offline_score_approval.json')
    if receipt.get('status')!='COMPLETE' or seal.get('status')!='ALL32_ONLINE_PREDICTIONS_SEALED':raise RuntimeError('All32 seal and scoring completion required')
    if approval.get('approved_by_root') is not True or approval.get('action')!='OFFLINE_TRUTH_SCORING_AFTER_ALL32_SEALED':raise RuntimeError('Root offline truth approval required')
    if approval['protocol_id']!=proto['id'] or approval['statistics_id']!=stats['id'] or float(approval['utc'])<float(seal['utc']):raise RuntimeError('Approval/seal identity or order')
    summary=read(score/'SUMMARY.json');rows=read(score/'PARENT_PREFIX_SCORES.json')
    fixed='FIXED_Q8_RAY';eps=str(stats['epsilon_primary'])
    chosen=[r for r in rows if r['method']==fixed and r['K_direction_budget']==16]
    failures=[r for r in chosen if not r['hits'][eps]]
    empty=[r for r in failures if r['best_material_MSE'] is None]
    if empty:worst=sorted(empty,key=lambda r:r['parent'])[0];worst_rule='POSTHOC_FIRST_EMPTY_FIXED_CLOUD_FAILURE (right-censored nearest error)'
    elif failures:worst=max(failures,key=lambda r:(r['best_material_MSE'],r['parent']));worst_rule='POSTHOC_LARGEST_FIXED_CLOUD_NEAREST_MSE_AMONG_COVERAGE_FAILURES'
    else:worst=None;worst_rule='NO_FIXED_COVERAGE_FAILURE'
    cases=[(proto['parents'][0]['parent'],'PREDECLARED_FIRST_PARENT')]
    if worst is not None and worst['parent']!=cases[0][0]:cases.append((worst['parent'],worst_rule))
    elif worst is not None:cases[0]=(cases[0][0],'PREDECLARED_FIRST_PARENT_AND_'+worst_rule)
    fixed_coverage=next(r for r in summary['coverage'] if r['method']==fixed and r['K']==16 and r['epsilon']==stats['epsilon_primary'])
    anchor_coverage=next(r for r in summary['coverage'] if r['method']==fixed and r['K']==0 and r['epsilon']==stats['epsilon_primary'])
    sealed={f['path']:f['bytes'] for p in seal['parents'] for f in p['files']}
    def public_path(p):
        p=Path(p)
        if str(p) not in sealed or p.stat().st_size!=sealed[str(p)]:raise RuntimeError('Unsealed prediction array '+str(p))
        return p
    out=P/'figures/e5';out.mkdir(parents=True,exist_ok=True);provenance=[]
    for parent,rule in cases:
        done=read(public_path(root/'parents'/parent/'COMPLETE.json'))
        anchor=material(public_path(done['anchor_path']))
        ast=read(public_path(root/'parents'/parent/'anchor/ANCHOR_STATUS.json'))
        records=sorted([r for r in done['rows'] if r['arm']==fixed and r['status']=='EXECUTED' and r['primary_qualified']],key=lambda r:r['slot'])
        cloud=[];nonanchor=[]
        if ast['legal'] and ast['noise_ball']:cloud.append(anchor)
        threshold=stats['uniqueness']['material_RMS_threshold']
        for r in records:
            m=material(public_path(r['selected_path']))
            if not cloud or min(rms(m,k) for k in cloud)>threshold:
                cloud.append(m);nonanchor.append((r['slot'],m))
        if len(nonanchor)>=3:reps=[nonanchor[0],nonanchor[len(nonanchor)//2],nonanchor[-1]]
        else:reps=nonanchor+[None]*(3-len(nonanchor))
        # Only now read this declared case's truth; posthoc rule already fixed.
        truth_path=Path(proto['private_output'])/'truth'/f'{parent}.npz'
        with np.load(truth_path,allow_pickle=False) as z:
            if str(z['parent'])!=parent:raise RuntimeError('Truth identity mismatch')
            truth=np.asarray(z['material'],np.complex128).copy()
        if truth.shape!=(576,) or not np.isfinite(truth).all():raise RuntimeError('Truth shape/finite')
        mean=np.mean(np.stack(cloud),axis=0) if cloud else anchor.copy()
        panels=[('GT (offline only)',truth),('Online anchor',anchor)]
        for name,rep in zip(['First qualified fixed','Middle-slot fixed','Last qualified fixed'],reps):
            panels.append(('No additional qualified member',None) if rep is None else (name+f' / slot{rep[0]}',rep[1]))
        panels.append(('Uniform cloud mean\nheuristic, not validated member' if cloud else 'Anchor fallback\nempty qualified cloud',mean))
        row=next(r for r in chosen if r['parent']==parent)
        fig,axes=plt.subplots(2,6,figsize=(16,6.3),layout='constrained')
        for col,(title,m) in enumerate(panels):
            for rr,(part,low,high,label) in enumerate([(np.real,.1,1.9,'Re susceptibility'),(np.imag,.04,.22,'Im susceptibility')]):
                ax=axes[rr,col];ax.set_xticks([]);ax.set_yticks([])
                if m is None:ax.set_facecolor('#eeeeee');ax.text(.5,.5,'No member',ha='center',va='center',transform=ax.transAxes)
                else:
                    im=ax.imshow(part(m).reshape(24,24),origin='lower',vmin=low,vmax=high,cmap='viridis')
                    if rr==0:ax.set_title(title+f'\noffline RMS={rms(m,truth):.4f}',fontsize=8.5)
                if col==0:ax.set_ylabel(label,fontsize=10)
        for rr,(low,high) in enumerate([(.1,1.9),(.04,.22)]):
            sm=plt.cm.ScalarMappable(cmap='viridis',norm=matplotlib.colors.Normalize(low,high))
            fig.colorbar(sm,ax=axes[rr,:].tolist(),shrink=.78,pad=.01)
        fig.suptitle(f"E5 posthoc material gallery — fixed tube coverage {fixed_coverage['successes']}/32; anchor {anchor_coverage['successes']}/32 (epsilon=.075)\n"
                     f"{parent}: {rule}; this parent's fixed-cloud hit={row['hits'][eps]}, qualified unique={len(cloud)}\n"
                     'Samples use first/middle/last online slot order, never nearest-GT selection; mean is not a posterior mean.',
                     fontsize=10)
        fname=f'gallery_{parent}'
        fig.savefig(out/(fname+'.png'),dpi=170);fig.savefig(out/(fname+'.pdf'));plt.close(fig)
        provenance.append(dict(parent=parent,case_selection=rule,representative_selection='first/middle/last distinct qualified nonanchor in online slot order',
            representative_slots=[None if r is None else r[0] for r in reps],qualified_unique=len(cloud),
            fixed_epsilon075_hit=row['hits'][eps],mean_semantics='uniform finite qualified cloud; may leave P1 decoder support; no physics requalification',
            truth_read_after_case_rule_fixed=True,online_predictions_changed=False,new_Maxwell_RHS=0))
    (out/'GALLERY_POSTHOC_PROVENANCE.json').write_text(json.dumps(dict(utc=time.time(),root_posthoc_authorized=True,
        scoring_complete=str(score/'SCORING_COMPLETE.json'),cases=provenance,plot_scope='Explanatory posthoc cases, not population evidence or online selection'),indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--posthoc-approved',action='store_true');a=p.parse_args();run(a.posthoc_approved)
