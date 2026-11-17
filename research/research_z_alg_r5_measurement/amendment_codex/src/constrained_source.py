"""Codex-only targeted development. Same sources/power cap, solved constraints.
No GT in design. Old data and discovered fold/root candidates only.
"""
from complete_quotient import *
from scipy.optimize import brentq
from experiment import solve_trf
import subprocess

def constrained_match(d,Q,cap=4.):
 def source(lam):
  a=np.linalg.solve(Q+lam*np.eye(len(Q)),d.conj());power=np.real(np.vdot(a,Q@a));return a/np.sqrt(power)
 a=source(0.);initial=np.linalg.norm(a);lam=0.
 if initial>cap:
  hi=1e-10
  while np.linalg.norm(source(hi))>cap:hi*=10
  lam=brentq(lambda z:np.linalg.norm(source(z))-cap,0.,hi,xtol=1e-16,rtol=1e-14);a=source(lam)
 return dict(a=a,lagrange_ratio=lam,power=float(np.real(np.vdot(a,Q@a))),feed_norm=float(np.linalg.norm(a)),unconstrained_feed=float(initial),objective=float(abs(d@a)))

def left_inverse_certificate(J):
 jm=np.array([[float(z.mid) for z in row] for row in J]);C=np.linalg.pinv(jm);M=[]
 for i in range(C.shape[0]):
  M.append([mp.iv.mpf(i==j)-sum((mp.iv.mpf(float(C[i,k]))*J[k][j] for k in range(len(J))),mp.iv.mpf(0)) for j in range(C.shape[0])])
 q=max(sum((upper(z) for z in row),mp.iv.mpf(0)) for row in M)
 return dict(certified_full_material_rank=bool(q<1),left_inverse_contraction_upper=float(q.b),midpoint_sigma_min=float(np.linalg.svd(jm,compute_uv=False)[-1]),left_inverse_infinity_norm=float(np.linalg.norm(C,np.inf)),scope='All materials in reported small convex box; this is not whole noisy feasible-set coverage')

def run():
 begin=time.process_time();book=Ledger('constrained_source_development');mp.iv.prec=128;g=geometry('T0');pos,Q=dictionary('T0',g);p=Physics(g,pos,book);y=np.load(BASE/'inputs/T0.npz')['old_noisy'];old=np.eye(p.L)[:,:p.s];roots=[np.array(m) for m in json.loads((BASE/'results/T0/noisy/Q.json').read_text())['accepted_roots']];fold=np.array(json.loads((OUT/'results/minimax_fold.json').read_text())['m']);_,Jold=p.response(fold,old,True);u=np.linalg.svd(Jold,full_matrices=False)[2][-1];_,Jdict=p.response(fold,np.eye(p.L),True);Jc=Jdict.reshape(p.L,2,p.n);d_fold=(Jc[:,0]+1j*Jc[:,1])@u;d_pair=(p.dual(roots[0])-p.dual(roots[1]))@p.B
 designs={key:constrained_match(d,Q) for key,d in [('pair_constrained',d_pair),('fold_constrained',d_fold)]}
 # Strictly fixed before obtaining new synthetic measurement. Both are separate policies.
 save(OUT/'results/constrained_source_design_freeze.json',dict(case='T0',exposure='POST_RESULT_DEVELOPMENT_NOT_NEW_BLIND',deployment_policy='fold_constrained',designs=designs,fold=fold,weak_material_direction=u,root_materials=roots,GT_used=False))
 traces={};intervals=[mp.iv.mpf([float(v-1e-6),float(v+1e-6)]) for v in fold]
 F,JI,_,_=oldcert._interval_forward_jacobian(intervals,p.H.tolist(),p.B.tolist(),p.R.tolist(),[0j]*p.L)
 for key,design in designs.items():
  a=design['a'];req=OUT/'results'/(key+'_request.npz');ob=OUT/'results'/(key+'_observed.npz');np.savez(req,a=a)
  subprocess.run([sys.executable,str(BASE/'src/broker.py'),'T0','noisy',key,str(req),str(ob)],check=True)
  yp=np.load(ob)['y'].reshape(-1);A=np.column_stack([old,a]);yy=np.r_[y,yp]
  qfit=solve_trf(p,A,yy,method=key+'_candidate_joint',budget=10.,initials=[np.full(p.n,15.25)]+roots)
  independent=solve_trf(p,A,yy,method=key+'_TRF64',budget=30.)
  Jnew=[]
  for i in range(p.n):Jnew.append(sum((ivc(al)*mp.iv.mpc(JI[2*l][i],JI[2*l+1][i]) for l,al in enumerate(a)),mp.iv.mpc(0)))
  Jaug=JI[:2*p.s]+[[v.real for v in Jnew],[v.imag for v in Jnew]];proof=left_inverse_certificate(Jaug)
  traces[key]=dict(observation=yp,source=design,candidate_joint=qfit,same_source_TRF64=independent,fold_rank_proof=proof,new_measurements=1,prediction_gap=abs(d_pair@a),source_precision='exact binary64 weights over exact rounded source dictionary in interval rank proof')
  save(OUT/'results'/(key+'_trace.json'),traces[key])
 save(OUT/'results/constrained_source_predictions_frozen.json',dict(files={str(f.relative_to(OUT)):digest(f) for f in (OUT/'results').glob('*constrained*') if f.is_file()},GT_not_read=True))
 # Separate post-prediction score, explicitly exposed development.
 truth=np.load(BASE/'private/truth.npz')['T0'];bank=np.load(BASE/'private/response_bank.npz');score={}
 for key,tr in traces.items():
  score[key]=dict(candidate_joint_RMS=float(np.linalg.norm(tr['candidate_joint']['best']-truth)/np.sqrt(p.n)),same_source_TRF64_RMS=float(np.linalg.norm(tr['same_source_TRF64']['best']-truth)/np.sqrt(p.n)),model_discrepancy=float(abs((bank['T0_shifted']-bank['T0_clean'])@tr['source']['a'])),gap=tr['prediction_gap'],fold_rank_proof=tr['fold_rank_proof'])
 save(OUT/'results/constrained_source_scores.json',dict(score=score,exposure='EXPOSED_DEVELOPMENT_ONLY',CPU_seconds=time.process_time()-begin,receipt=book.receipt()))
 print(serial(score),flush=True)
if __name__=='__main__':run()
