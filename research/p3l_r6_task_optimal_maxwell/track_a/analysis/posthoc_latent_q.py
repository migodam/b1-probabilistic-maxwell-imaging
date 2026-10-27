"""POSTHOC_FUNCTIONAL_ANALYSIS: existing GH nodes + exact prior fibre, zeroRHS."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
from pathlib import Path
import json,csv,time,resource,functools
import numpy as np
from scipy.special import ndtr,betaincinv
from numpy.polynomial.legendre import leggauss
HERE=Path(__file__).resolve().parent;AOUT=HERE.parent
LO,HI=-4.,4.;Z0=ndtr(HI)-ndtr(LO);A=1/1.3;PH=(ndtr(-A)-ndtr(LO))/Z0;ULO,UHI=ndtr(LO),ndtr(HI);RLO,RHI=.8+ULO,.8+UHI;ILO,IHI=.03+.15*ULO,.03+.15*UHI;AREA=(RHI-RLO)*(IHI-ILO)
def H(s):
 t=np.clip((np.asarray(s)+1)/2,0,1);return t*t*t*(10+t*(-15+6*t))
@functools.lru_cache(None)
def gl(n):return leggauss(n)
def fibre(R,I,n):
 lo=max(R/RHI,I/IHI);hi=min(R/RLO,I/ILO,1.)
 if not 0<lo<hi:return dict(density=0.,transition=0.,pure=0.,q_mean=0.,q_second=0.,q_variance=0.,deltaR_mean=0.,q_deltaR_cov=0.,phase_gt_half=0.)
 zl=(2*betaincinv(3,3,lo)-1)/1.3;zh=(2*betaincinv(3,3,hi)-1)/1.3;x,w=gl(n);z=zl+(zh-zl)*(x+1)/2;q=H(1.3*z);weight=w*(zh-zl)/2*np.exp(-z*z/2)/np.sqrt(2*np.pi)/Z0/q**2/AREA;transition=weight.sum();pure=PH/AREA if RLO<R<RHI and ILO<I<IHI else 0.;d=transition+pure;m1=(weight@q+pure)/d;m2=(weight@(q*q)+pure)/d;inv=(weight@(1/q)+pure)/d;dr=R*inv
 if zh<=0:phase=0.
 elif zl>=0:phase=1.
 else:
  z=zh*(x+1)/2;q2=H(1.3*z);wp=w*zh/2*np.exp(-z*z/2)/np.sqrt(2*np.pi)/Z0/q2**2/AREA;phase=(wp.sum()+pure)/d
 return dict(density=float(d),transition=float(transition),pure=float(pure),q_mean=float(m1),q_second=float(m2),q_variance=float(m2-m1*m1),deltaR_mean=float(dr),q_deltaR_cov=float(R-m1*dr),phase_gt_half=float(phase))
def write(name,v):(HERE/name).write_text(json.dumps(v,indent=2,default=lambda v:v.item() if isinstance(v,np.generic) else v.tolist())+'\n')
def main():
 tick=time.process_time();obs=np.load(AOUT/'OBSERVATIONS.npz');raw=np.load(AOUT/'raw/S3_source_16_9.npz');host_response=raw['prediction'][0];assert raw['V'][0]==0 and raw['layer'][0]==0
 local=json.loads((AOUT/'LOCAL_REFINEMENT_RESULTS.json').read_text())['rows'];records=[];nodes={}
 for case in ['REGULAR','PURE_OVERLAP']:
  ix=obs['names'].tolist().index(case);y=obs['observation'][ix];sigma=obs['sigma'][ix];hostlogZ=np.log(PH)-np.sum(abs(host_response-y)**2)/sigma**2
  for order in [15,21]:
   saved=np.load(AOUT/f'raw/{case}_material_GH{order}.npz');V=saved['V'];w96=saved['weights'];ref=next(v for v in local if v['observation']==case and v['order']==order);byorder={n:[fibre(*v,n) for v in V] for n in [96,128]}
   for nf in [96,128]:
    f=byorder[nf];d=np.array([v['density'] for v in f]);d96=np.array([v['density'] for v in byorder[96]]);weights=w96*d/d96;p=weights/weights.sum();logZ=ref['log_evidence']+np.log(weights.sum()/w96.sum());logtotal=np.logaddexp(logZ,hostlogZ);hp=float(np.exp(hostlogZ-logtotal));nonhost=float(np.exp(logZ-logtotal));mq=np.array([v['q_mean'] for v in f]);vqin=np.array([v['q_variance'] for v in f]);m1=nonhost*(p@mq);m2=nonhost*(p@np.array([v['q_second'] for v in f]));variance=m2-m1*m1;within=nonhost*(p@vqin);between=nonhost*(p@(mq-m1)**2)+hp*m1*m1;pure=nonhost*(p@np.array([v['pure']/v['density'] for v in f]));phase=nonhost*(p@np.array([v['phase_gt_half'] for v in f]));dr=nonhost*(p@np.array([v['deltaR_mean'] for v in f]))+hp*(RLO+RHI)/2;matR=nonhost*(p@V[:,0]);covQD=matR-m1*dr;mode=fibre(*saved['mode'],nf)
    rec=dict(case=case,GH_order=order,fibre_GL=nf,E_Q=float(m1),E_Q2=float(m2),Var_Q=float(variance),within_fibre_Q_variance=float(within),between_material_Qmean_variance=float(between),total_variance_identity_gap=float(variance-within-between),P_Q_equals1=float(pure),P_Q_gt_half=float(phase),E_native_deltaR=float(dr),Cov_Q_deltaR=float(covQD),material_mean_R=float(matR),mode_R=float(saved['mode'][0]),mode_I=float(saved['mode'][1]),mode_fibre_E_Q=mode['q_mean'],mode_fibre_Var_Q=mode['q_variance'],mode_fibre_P_Q_equals1=mode['pure']/mode['density'],mode_fibre_Cov_Q_deltaR=mode['q_deltaR_cov'],log_host_evidence=float(hostlogZ),log_local_nonhost_evidence=float(logZ),log_host_to_local_posterior_mass=float(hostlogZ-logtotal),numeric_host_mass=hp,host_contribution_rule='Q=0 atom retained via logevidence; underflow numeric0 is not proofzero; deltaR host mean1.3',minimum_possible_Q_on_GH_nodes=float(min(max(R/RHI,I/IHI) for R,I in V)),maximum_possible_Q_on_GH_nodes=float(max(min(R/RLO,I/ILO,1.) for R,I in V)),max_inner_density_relative_change_96to128=float(max(abs(v['density']/u['density']-1) for v,u in zip(byorder[128],byorder[96]))),scope='POSTHOC_FUNCTIONAL_ANALYSIS; fixedS3 conditionalsubmodel; localGH Q1/Q2, no globaltailcertificate')
    records.append(rec);nodes[f'{case}_GH{order}_fibre{nf}_Qmean']=mq;nodes[f'{case}_GH{order}_fibre{nf}_Qvar']=vqin;nodes[f'{case}_GH{order}_fibre{nf}_weights']=nonhost*p
 comparisons=[]
 for case in ['REGULAR','PURE_OVERLAP']:
  lo=next(r for r in records if r['case']==case and r['GH_order']==15 and r['fibre_GL']==128);hi=next(r for r in records if r['case']==case and r['GH_order']==21 and r['fibre_GL']==128);inner=next(r for r in records if r['case']==case and r['GH_order']==21 and r['fibre_GL']==96)
  comparisons.append(dict(case=case,GH15to21_abs_E_Q=abs(lo['E_Q']-hi['E_Q']),GH15to21_abs_Var_Q=abs(lo['Var_Q']-hi['Var_Q']),GH15to21_abs_P_Q_equals1=abs(lo['P_Q_equals1']-hi['P_Q_equals1']),GL96to128_abs_E_Q=abs(inner['E_Q']-hi['E_Q']),GL96to128_abs_Var_Q=abs(inner['Var_Q']-hi['Var_Q']),qualification='observed refinement errors only; no inherited primarygate or globaltail certificate'))
 np.savez_compressed(HERE/'POSTHOC_LATENT_Q_NODE_MOMENTS.npz',**nodes)
 with (HERE/'POSTHOC_LATENT_Q.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
 write('POSTHOC_LATENT_Q.json',dict(status='POSTHOC_FUNCTIONAL_ANALYSIS',records=records,comparisons=comparisons,new_physics_RHS=0,new_fits=0,A_total_RHS_unchanged=157696,global_tail_certificate=False,material_reference_qualification_not_automatically_inherited=True,CPU_seconds=time.process_time()-tick,CPU_threads=1,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
 print(json.dumps(dict(records=[r for r in records if r['GH_order']==21 and r['fibre_GL']==128],comparisons=comparisons),indent=2))
if __name__=='__main__':main()
