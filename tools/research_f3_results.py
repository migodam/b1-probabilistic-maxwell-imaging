"""Saved-array accounting, offline material scores and scientific plots. No solver."""
from pathlib import Path
import json,sys,time,collections,os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from finite_feedback.common import dump,rms
from prior_constrained import prior
P=ROOT/'research/research_f3_quotient_exploration';OLD=Path(os.environ.get('B1_F3_LEGACY_ROOT',str(ROOT.parent/'b1-research-f/research/research_f_intrinsic_transport')))

def arr(row):
 with np.load(ROOT/row['array']) as z:return z['chi'],z['prediction']
def valid(r):return r['data_feasible'] and r['likelihood_competitive'] and r['task_displacement'] and r['structure_check'] and r['posterior_supported_positive_density']
def save(fig,name):
 (P/'figures').mkdir(exist_ok=True);fig.savefig(P/f'figures/{name}.png',dpi=175,bbox_inches='tight');fig.savefig(P/f'figures/{name}.pdf',bbox_inches='tight');plt.close(fig)
def main():
 if not (P/'STAGE_CLOSED.json').exists():raise RuntimeError('no offlineGT until physical execution closes')
 rows=json.loads((P/'results/RECORDS.json').read_text());paths=json.loads((P/'results/PATH_RECORDS.json').read_text());ledger=[json.loads(x) for x in (P/'RESOURCE_LEDGER.jsonl').read_text().splitlines()]
 costs={k:sum(r.get('amount',0) for r in ledger if r['event']==event) for k,event in [('states','STATE'),('complex_RHS','RHS')]}
 costs.update(controlled_wall_seconds=sum(r.get('controlled_wall_seconds',0) for r in ledger if r['event']=='JOB_END'),full_state_forward_RHS=sum(r.get('amount',0) for r in ledger if r['event']=='RHS' and r.get('role','').startswith('full')),
  adjoint_RHS=sum(r.get('amount',0) for r in ledger if r['event']=='RHS' and r.get('role','').endswith('_adjoint')),Hessian_RHS=sum(r.get('amount',0) for r in ledger if r['event']=='RHS' and 'F3_Hessian_current' in r.get('role','')),
  LU_factorizations=sum(r.get('counts',{}).get('LU_factorizations',0) for r in ledger if r['event']=='ACTION'),caps=dict(states=64,complex_RHS=2048,wall_seconds=1200),GPU_training=0,Mie=0,new3D=0,endpoint_files=len(list((P/'results/endpoints').glob('*.npz'))))
 bymethod=[]
 for case in sorted(set(r['case'] for r in rows)):
  for method in ['A_frozen','B_updated','C_native','D_intrinsic_feedback','native_random','MH']:
   rr=[r for r in rows if r['case']==case and r['method']==method]
   bymethod.append(dict(case=case,method=method,endpoint_records=len(rr),data_feasible=sum(r['data_feasible'] for r in rr),competitive=sum(r['likelihood_competitive'] for r in rr),task_positive_prior_competitive=sum(valid(r) for r in rr),task_target_matched=sum(r.get('target_match',False) for r in rr),evaluated_MH_proposals_not_samples=method=='MH'))
 eligible=[r for r in rows if valid(r)];unique=[]
 for r in eligible:
  x,_=arr(r)
  if all(rms(x-y['chi'])>1e-9 for y in unique):unique.append(dict(row=r,chi=x))
 pairwise=np.array([[rms(a['chi']-b['chi']) for b in unique] for a in unique])
 dump(P/'results/ELIGIBLE_ENDPOINTS.json',dict(unique_material_endpoints=[dict(case=z['row']['case'],method=z['row']['method'],q=z['row']['q'],array=z['row']['array'],delta_Phi=z['row']['delta_Phi'],log_prior=z['row']['log_prior'],RMS=z['row']['material_RMS']) for z in unique],dedup_RMS_tolerance=1e-9,pairwise_RMS=pairwise,posterior_modes=None,posterior_mass=None))
 dump(P/'RESOURCE_SUMMARY.json',costs)
 # Two contexts of the same exposed parent, shared color scales, GT only here.
 evaluations=[];panels=[]
 for case in sorted(set(r['case'] for r in rows)):
  f=json.loads((P/f'results/F15/{case}_FINAL.json').read_text())
  for path in sorted((P/'results/endpoints').glob(case+'_F15_*.npz')):
   with np.load(path) as z:
    if np.array_equal(z['q'],f['q']):base=z['chi'];basepred=z['prediction'];break
  oldrow=next(r for r in json.loads((OLD/'results/F1/SUMMARY.json').read_text()) if r['case']==case)
  with np.load(OLD/f'results/F1/{case}/state_11.npz') as z:oldbase=z['chi']
  with np.load(OLD/f'inputs/offline_only/{case}.npz') as z:gt=z['truth']
  candidates=[r for r in rows if r['case']==case and r['method']=='native_random' and r['posterior_supported_positive_density']]
  best=min(candidates,key=lambda r:r['Phi']);image,_=arr(best)
  dfirst=next(r for r in rows if r['case']==case and r['method']=='D_intrinsic_feedback');dimage,_=arr(dfirst)
  panels.append((case,[oldbase,base,image,dimage,gt],['Old F1','F3 likelihood anchor','Native candidate','Intrinsic/feedback','Offline truth'],best,dfirst))
  scores={name:dict(complex_MAE=float(np.mean(abs(x-gt))),complex_RMS=rms(x-gt),real_MAE=float(np.mean(abs((x-gt).real))),imag_MAE=float(np.mean(abs((x-gt).imag)))) for name,x in zip(['old_F1','F3_anchor','native_candidate','intrinsic_candidate'],[oldbase,base,image,dimage])}
  structures=[]
  for z in unique:
   if z['row']['case']!=case:continue
   structures.append(dict(method=z['row']['method'],**prior.structural_distance(prior.rho_field(base,prior.background('weak' if 'weak' in case else 'ring')),prior.rho_field(z['chi'],prior.background('weak' if 'weak' in case else 'ring')))))
  evaluations.append(dict(case=case,known_family_global_geometry=True,same_parent_two_backgrounds=True,scores=scores,native_selection='minimumPhi among positiveprior-density native candidates; noGT',native_candidate_deltaPhi=best['delta_Phi'],intrinsic_candidate_deltaPhi=dfirst['delta_Phi'],structures=structures))
 real_limits=(min(x.real.min() for _,ims,*_ in panels for x in ims),max(x.real.max() for _,ims,*_ in panels for x in ims));imag_limits=(min(x.imag.min() for _,ims,*_ in panels for x in ims),max(x.imag.max() for _,ims,*_ in panels for x in ims))
 for case,images,titles,best,dfirst in panels:
  fig,axes=plt.subplots(2,5,figsize=(16,6.5),layout='constrained')
  for j,(x,title) in enumerate(zip(images,titles)):
   for i,(part,lim) in enumerate([(x.real,real_limits),(x.imag,imag_limits)]):im=axes[i,j].imshow(part.reshape(24,24),origin='lower',vmin=lim[0],vmax=lim[1],cmap='viridis');axes[i,j].set_xticks([]);axes[i,j].set_yticks([])
   axes[0,j].set_title(title)
  axes[0,0].set_ylabel('Re chi');axes[1,0].set_ylabel('Im chi')
  for i,lim in enumerate([real_limits,imag_limits]):fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(*lim),cmap='viridis'),ax=axes[i,:],fraction=.025)
  fig.suptitle(case+' | candidate-only, no posterior masses\nNative deltaPhi=%.3f ; intrinsic deltaPhi=%.3f ; same color scale, GT offline only'%(best['delta_Phi'],dfirst['delta_Phi']))
  save(fig,'01_material_'+('weak' if 'weak' in case else 'ring'))
 # DeltaPhi scatter and matched-radius bars.
 methods=['A_frozen','B_updated','C_native','D_intrinsic_feedback','native_random','MH'];markers=['o','s','^','D','x','+'];fig,axes=plt.subplots(1,2,figsize=(14,5),layout='constrained')
 for ax,case in zip(axes,sorted(set(r['case'] for r in rows))):
  for method,marker in zip(methods,markers):
   rr=[r for r in rows if r['case']==case and r['method']==method];ax.scatter([r['material_RMS'] for r in rr],[r['delta_Phi'] for r in rr],label=method,marker=marker,s=55)
  ax.axhline(2,color='black',ls='--',label='competitive deltaPhi<=2');ax.axvline(.00527555,color='grey',ls=':');ax.set(xlabel='actual material RMS from F3 anchor',ylabel='deltaPhi vs fixed old F1',title=case.split('_')[1]);ax.set_ylim(-.5,85);ax.legend(fontsize=8)
 fig.suptitle('Noise-ball feasibility does not imply likelihood competitiveness. MH points are proposals, not all retained samples.');save(fig,'02_RMS_likelihood')
 fig,axes=plt.subplots(1,2,figsize=(13,4.5),layout='constrained')
 for ax,case in zip(axes,sorted(set(r['case'] for r in rows))):
  chain=json.loads((P/f'results/{case}_MH_CHAIN.json').read_text());ch=chain['steps'];ax.step(range(1,9),[r['material_RMS_from_initial'] for r in ch],where='post',label='retained chain RMS')
  ax.axhline(.00527555,color='black',ls='--',label='task scale');ax.set(xlabel='MH step (8 only)',ylabel='retained material RMS',title=case.split('_')[1]+': 3/8 accepted');ax.legend()
 fig.suptitle('Correct numerical-target MH kernel; eight steps do not establish mixing or calibrated probabilities.');save(fig,'03_MH_chain')
 fig,ax=plt.subplots(figsize=(9,5),layout='constrained')
 for pid in sorted(set(r['path_id'] for r in paths)):
  rr=[r for r in paths if r['path_id']==pid];plan=next(r for r in json.loads((P/'configs/PATH_VALIDATION.json').read_text())['paths'] if r['id']==pid)
  ax.plot([0]+[r['t'] for r in rr]+[1],[0]+[r['delta_Phi'] for r in rr]+[plan['endpoint_deltaPhi']],marker='o',label=f"path{pid}: {plan['method']}")
 ax.axhline(2,color='black',ls='--');ax.set(xlabel='linear quotient-coordinate interpolation t',ylabel='deltaPhi vs fixed reference',title='All sampled midpoints noise-feasible; none likelihood-competitive');ax.legend();save(fig,'04_path_obstruction')
 metadata=dict(scope='OFFLINE_ONLY_CLOSED_STAGE',new_physics=0,new_RHS=0,material_evaluation=evaluations,unique_eligible_endpoint_count=len(unique),pairwise_eligible_RMS=pairwise,prior_audit='PRIOR_AUDIT.json',method_summary=bymethod,path_midpoints=dict(count=len(paths),data_feasible=sum(r['data_feasible'] for r in paths),competitive=sum(r['likelihood_competitive'] for r in paths)),resources=costs,
  posterior_mass=None,calibration='NOT_ESTIMATED',outside_model_mass=None,support='conditional original smoothP0, not general 2D/3D')
 dump(P/'METRICS.json',metadata);dump(P/'OFFLINE_EVALUATION.json',dict(evaluations=evaluations,new_physics=0,new_RHS=0))
 print(json.dumps(dict(costs=costs,unique_eligible=len(unique),pairwise_RMS=pairwise.tolist(),methods=bymethod)))
if __name__=='__main__':main()

def extra_plots():
    eligible=json.loads((P/'results/ELIGIBLE_ENDPOINTS.json').read_text())['unique_material_endpoints'];case='p250073_weak_a1_snr20'
    f=json.loads((P/f'results/F15/{case}_FINAL.json').read_text())
    for path in sorted((P/'results/endpoints').glob(case+'_F15_*.npz')):
        with np.load(path) as z:
            if np.array_equal(z['q'],f['q']):base=z['chi'];break
    differences=[]
    for r in eligible:
        x,_=arr(r);differences.append(x-base)
    fig,axes=plt.subplots(2,len(eligible),figsize=(12,6),layout='constrained')
    for i,part in enumerate([np.real,np.imag]):
        cap=max(np.max(abs(part(d))) for d in differences)
        for j,(d,r) in enumerate(zip(differences,eligible)):
            im=axes[i,j].imshow(part(d).reshape(24,24),origin='lower',vmin=-cap,vmax=cap,cmap='RdBu_r');axes[i,j].set_xticks([]);axes[i,j].set_yticks([])
            if i==0:axes[i,j].set_title(r['method']+'\ndeltaPhi=%.3f'%r['delta_Phi'])
        axes[i,0].set_ylabel('delta Re chi' if i==0 else 'delta Im chi');fig.colorbar(im,ax=axes[i,:],fraction=.025)
    fig.suptitle('Three nearby endpoints, not three modes: pairwise RMS .000639-.001100; no topology-count changes')
    save(fig,'05_candidate_differences')
    rr=json.loads((P/'results/RECORDS.json').read_text());fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    for ax,bg in zip(axes,['ring','weak']):
        rs=[next(r for r in rr if bg in r['case'] and r['method']==method and r['phase']=='F2_C1') for method in ['A_frozen','B_updated','C_native','D_intrinsic_feedback']]
        y=[r['delta_Phi'] for r in rs];bars=ax.bar(['Frozen','Updated','Native','Intrinsic'],y,color=['#3e7e9a','#e5ac52','#419678','#c76a64']);ax.bar_label(bars,fmt='%.3f');ax.axhline(2,ls='--',color='black');ax.set(title=bg+' at matched s',ylabel='deltaPhi (lower is better)');ax.set_ylim(0,max(y)*1.2)
    fig.suptitle('Same material RMS .00527556: coordinate/retraction gain; no stable high-order increment')
    save(fig,'06_matched_distance')
if __name__=='__main__':extra_plots()
