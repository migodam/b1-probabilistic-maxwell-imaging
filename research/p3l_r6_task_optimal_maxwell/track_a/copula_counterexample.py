"""ROOT-EXECUTED bounded Full-physics copula control; do not import for a run.

At fixed saved material mode V, sample the prior fibre jointly, then apply one
fixed permutation of Q relative to D. This is a quadrature-conditioned coupling
control, never a posterior sample/calibration/coverage experiment. Max528RHS;
no adjoint/J calls. The implementing worker only syntax-checks this script.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[_key]='1'
import sys,json,time,resource
from pathlib import Path
import numpy as np
from scipy.special import ndtr,ndtri,betaincinv
from numpy.polynomial.legendre import leggauss
OUT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(OUT))
from shared.physics import Book,Target,SIGMA,write,save
TRACK=OUT/'track_a';RAW=TRACK/'raw/copula_counterexample'
MAX_CALLS=66;MAX_RHS=528;N=32;SEED=202610270612
LO,HI=-4.,4.;Z0=ndtr(HI)-ndtr(LO);ZEDGE=1/1.3;PH=(ndtr(-ZEDGE)-ndtr(LO))/Z0
ULO,UHI=ndtr(LO),ndtr(HI);RLO,RHI=.8+ULO,.8+UHI;ILO,IHI=.03+.15*ULO,.03+.15*UHI;AREA=(RHI-RLO)*(IHI-ILO)
def H(s):
 t=np.clip((np.asarray(s)+1)/2,0,1);return t*t*t*(10+t*(-15+6*t))
def conditional_fibre(v,n=128):
 R,I=v.real,v.imag;ql=max(R/RHI,I/IHI);qh=min(R/RLO,I/ILO,1.)
 if not 0<ql<qh:raise ValueError('saved nonzero mode outside conditional fibre image')
 zl=(2*betaincinv(3,3,ql)-1)/1.3;zh=(2*betaincinv(3,3,qh)-1)/1.3;x,w=leggauss(n);z=zl+(zh-zl)*(x+1)/2;q=H(1.3*z);wt=w*(zh-zl)/2*np.exp(-z*z/2)/np.sqrt(2*np.pi)/Z0/q**2/AREA
 pure=PH/AREA if RLO<R<RHI and ILO<I<IHI else 0.;density=wt.sum()+pure
 Q=np.r_[q,1.];prob=np.r_[wt,pure]/density;mom={str(k):float(prob@(Q**k)) for k in [-2,-1,1,2]};m1,m2,im1,im2=[mom[str(k)] for k in [1,2,-1,-2]];factor=m1*im1;shift=v*(factor-1)
 return dict(z=z,Q=Q,prob=prob,density=float(density),pure_probability=float(prob[-1]),moments=mom,exact_marginal_product_factor=float(factor),CS_excess=float(factor-1),conditional_Cov_Q_deltaR=float(v.real*(1-factor)),exact_product_mean_increment=[float((v*factor).real),float((v*factor).imag)],exact_product_mean_bias=[float(shift.real),float(shift.imag)],exact_product_complex_material_variance=float(abs(v)**2*(m2*im2-factor**2)),scope='GL128 prior conditional fibre at fixed material V; not complete posterior')
def draw_pairs(v,fibre,seed):
 rng=np.random.default_rng(seed);index=rng.choice(len(fibre['Q']),size=N,replace=True,p=fibre['prob']);Q=fibre['Q'][index].copy();z=np.empty(N);transition=index<len(fibre['z']);z[transition]=fibre['z'][index[transition]];n_pure=int((~transition).sum())
 if n_pure:z[~transition]=ndtri(ndtr(ZEDGE)+rng.uniform(size=n_pure)*(ndtr(HI)-ndtr(ZEDGE)))
 delta=v/Q;native=np.zeros((N,66));native[:,0]=z;native[:,64]=delta.real-.8;native[:,65]=(delta.imag-.03)/.15
 # Exactly one fixed random permutation, with fixed points allowed; no resampling.
 permutation=rng.permutation(N);factorized=native.copy();factorized[:,0]=z[permutation]
 return native,factorized,Q,delta,permutation,index

def main():
 protocol=TRACK/'COPULA_COUNTEREXAMPLE_PROTOCOL.json';result_path=TRACK/'COPULA_COUNTEREXAMPLE_RESULTS.json'
 if result_path.exists():raise FileExistsError('existing result; no overwrite or implicit physics replay')
 modes={}
 for case in ['REGULAR','PURE_OVERLAP']:
  with np.load(TRACK/f'raw/{case}_material_GH21.npz') as z:modes[case]=complex(*z['mode'])
 before=Book.totals();A0=int(before.get('stages',{}).get('A',{}).get('RHS',0));caps=json.loads((OUT/'configs/AUTHORIZATION.json').read_text())
 if before['active_jobs'] or before['active_actions']:raise RuntimeError('root must wait until existing physics job closes; no concurrent worker')
 if A0+MAX_RHS>caps['stage_caps']['A']:raise RuntimeError('insufficient remaining A reservation; no budget expansion')
 write(protocol,dict(status='FROZEN_BEFORE_PHYSICS',role='CONTROLLED_FIXED_MATERIAL_COPULA_COUNTEREXAMPLE',saved_mode_sources={c:str(TRACK/f'raw/{c}_material_GH21.npz') for c in modes},modes={c:[v.real,v.imag] for c,v in modes.items()},GT_read=False,Dtruth_read=False,seed=SEED,case_seed_rule='202610270612+caseindex, REGULAR0/PURE_OVERLAP1',pairs_per_case=32,fibre='transitionGL128 phi(z)/Z0/Q²/rectangle_area + purePH/rectangle_area; host absent sincefixedVnonzero',pure_z='conditional standardNormal truncated[zedge,4]; allotherfieldcoordinates0',factorization='one sameRNG permutation of sampled Q/z relative to sampled D; empirical-marginal coupling control, not iid product posterior draws',max_reference_calls=2,max_factorized_calls=64,max_forward_calls=66,max_new_RHS=528,adjoint_calls=0,exact_material_cache='byte-identical decoded complex128 material only; no rounding or phantomdiscount',A_RHS_before=A0,stage_A_cap=caps['stage_caps']['A'],response_metric='complexL2(Ffactor-Fjoint)/sigma; RMS and sqrt2-whiteL2 also saved',scope='neither Γcalibration nor fullposterior sampling; no Full66claim; no primaryAthreshold changes'))
 b=Book('A','fixedmaterial_fibre_copula_counterexample',role='CONTROLLED_COPULA_PHYSICS',max_wall=120);status='FAILED';error=None;calls=0;tick=time.perf_counter();rows=[]
 try:
  target=Target(b);roi=target.prior.mask
  if not np.allclose(target.prior.A[roi,0]/target.prior.temperature,1.3) or target.cfg['field_mean']!=0:raise ValueError('fixed uniform-field P1decoder contract changed')
  if not np.allclose(target.prior.real_interval,[.8,1.8]) or not np.allclose(target.prior.imag_interval,[.03,.18]):raise ValueError('frozen property intervals changed')
  def paid(q,purpose):
   nonlocal calls
   used=int(Book.totals()['stages']['A']['RHS'])-A0
   if calls>=MAX_CALLS or used+8>MAX_RHS:raise RuntimeError('66call/528RHS local hard stop before nextquery')
   calls+=1;return target.evaluate(q,purpose)
  for case_index,(case,v) in enumerate(modes.items()):
   fibre=conditional_fibre(v);joint,factorized,Q,delta,perm,index=draw_pairs(v,fibre,SEED+case_index);jm=[];fm=[]
   for q in np.r_[joint,factorized]:
    if not target.prior.legal(q) or abs(q[0])>4+1e-12 or np.max(abs(ndtri(q[64:])))>4+1e-9:raise ValueError('control left frozen truncated-prior support')
   for q in joint:jm.append(target.prior.decode(q)[0])
   for q in factorized:fm.append(target.prior.decode(q)[0])
   jm=np.asarray(jm);fm=np.asarray(fm);expected=target.prior.host+v*roi;gauge_error=float(np.max(abs(jm-expected)))
   if gauge_error>2e-13:raise RuntimeError('joint fibre does notdecode requestedfixedmaterial')
   # One jointreference; remaining jointstates are decoder-only same-material controls.
   ref=paid(joint[0],f'{case}_joint_material_reference');refpred=ref.state.data.copy();refchi=ref.state.chi.copy();target.persist_cache(ref,RAW/f'{case}_REFERENCE_CACHE.npz');cache={refchi.tobytes():(refpred,ref.key)};del ref
   predictions=[];endpoints=[];did_query=[]
   for k,q in enumerate(factorized):
    key=fm[k].tobytes()
    if key in cache:prediction,endpoint=cache[key];did_query.append(False)
    else:
     point=paid(q,f'{case}_factorized_QD_{k:02d}');prediction=point.state.data.copy();endpoint=point.key
     if not np.array_equal(point.state.chi,fm[k]):raise RuntimeError('predecoded/query material mismatch')
     cache[key]=(prediction,endpoint);did_query.append(True);del point
    predictions.append(prediction);endpoints.append(endpoint)
   predictions=np.asarray(predictions);errors=predictions-refpred;norm=np.sqrt(np.sum(abs(errors)**2,axis=(1,2)))/SIGMA;rms=np.sqrt(np.mean(abs(errors)**2,axis=(1,2)))/SIGMA;jv=np.mean(jm[:,roi]-target.prior.host,axis=1);fv=np.mean(fm[:,roi]-target.prior.host,axis=1);jmean=jv.mean();fmean=fv.mean();jvar=float(np.mean(abs(jv-jmean)**2));fvar=float(np.mean(abs(fv-fmean)**2));sample_product=Q.mean()*np.mean(1/Q);sample_target=v*sample_product
   save(RAW/f'{case}.npz',q_joint_native=joint,q_factorized_native=factorized,Q_joint=Q,Q_factorized=Q[perm],D_joint_complex=delta,permutation=perm,conditional_GL_indices=index,joint_material=jm,factorized_material=fm,joint_reference_q=joint[0],joint_reference_material=refchi,joint_reference_prediction=refpred,factorized_prediction=predictions,prediction_difference=errors,prediction_difference_over_sigma=errors/SIGMA,complex_L2_over_sigma=norm,complex_RMS_over_sigma=rms,real_whitened_L2=norm*np.sqrt(2),factorized_query_performed=np.array(did_query),factorized_endpoint=np.array(endpoints),conditional_Q_nodes=fibre['Q'],conditional_Q_prob=fibre['prob'],saved_material_mode=np.array([v.real,v.imag]))
   row=dict(case=case,seed=SEED+case_index,pairs=32,reference_calls=1,factorized_fresh_calls=sum(did_query),factorized_exact_material_cache_hits=32-sum(did_query),joint_gauge_max_abs_decoder_error=gauge_error,joint_actual_mean_increment=[float(jmean.real),float(jmean.imag)],factorized_empirical_mean_increment=[float(fmean.real),float(fmean.imag)],empirical_mean_material_shift=[float((fmean-jmean).real),float((fmean-jmean).imag)],empirical_mean_shift_ROI_RMS=float(abs(fmean-jmean)),joint_empirical_complex_material_variance=jvar,factorized_empirical_complex_material_variance=fvar,empirical_material_variance_defect=fvar-jvar,empirical_marginal_product_factor=float(sample_product),empirical_product_of_marginals_mean_increment=[float(sample_target.real),float(sample_target.imag)],single_permutation_mean_vs_empirical_product=[float((fmean-sample_target).real),float((fmean-sample_target).imag)],Full_response_complex_L2_over_sigma_mean=float(norm.mean()),Full_response_complex_L2_over_sigma_p90=float(np.quantile(norm,.9)),Full_response_complex_L2_over_sigma_max=float(norm.max()),Full_response_complex_RMS_over_sigma_mean=float(rms.mean()),conditional_prior_quadrature={k:w for k,w in fibre.items() if k not in ['z','Q','prob']},qualification='fixed savedmaterial mode, quadrature-conditioned priorfibre/marginal-copula control; not Γcalibration, posterior samples, independentworlds, or Full66claim',source_data=str(RAW/f'{case}.npz'),reference_cache=str(RAW/f'{case}_REFERENCE_CACHE.npz'))
   rows.append(row);write(TRACK/'COPULA_COUNTEREXAMPLE_PARTIAL.json',dict(status='PARTIAL',rows=rows,actual_new_RHS=int(Book.totals()['stages']['A']['RHS'])-A0,attempted_forward_calls=calls));print(json.dumps(dict(case=case,calls=calls,charged_RHS=int(Book.totals()['stages']['A']['RHS'])-A0)),flush=True)
   if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>2_000_000_000:raise MemoryError('2GB local resource stop')
  new=int(Book.totals()['stages']['A']['RHS'])-A0
  if new>MAX_RHS or b.counts.get('actual_derivative_RHS',0)!=0 or new!=8*calls:raise RuntimeError('unexpected derivativeorRHS accounting')
  status='COMPLETE';write(result_path,dict(status='COMPLETE_CONTROLLED_COPULA_COUNTEREXAMPLE',rows=rows,A_RHS_before=A0,A_RHS_after=int(Book.totals()['stages']['A']['RHS']),actual_new_RHS=new,attempted_forward_calls=calls,max_forward_calls=66,max_RHS=528,adjoint_calls=0,wall_seconds=time.perf_counter()-tick,GT_read=False,Dtruth_read=False,fullposterior_sampling_or_Gamma_calibration=False,source_model='fixed activeS3 truncatedprior only; V fromsavedGH21 likelihoodmode',primaryAthresholds_modified=False,final_scientific_review='parent-owned'))
 except BaseException as exc:error=repr(exc);write(TRACK/'COPULA_COUNTEREXAMPLE_FAILURE.json',dict(status='FAILED',error=error,attempted_calls=calls,actual_new_RHS=int(Book.totals()['stages']['A']['RHS'])-A0,rows_completed=rows));raise
 finally:b.close(status,error)
if __name__=='__main__':main()
