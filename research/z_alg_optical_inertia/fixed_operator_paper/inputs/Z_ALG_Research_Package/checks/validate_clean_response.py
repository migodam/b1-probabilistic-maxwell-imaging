from pathlib import Path
import numpy as np,json
from scipy.linalg import eigvalsh
z=np.load(Path(__file__).resolve().parent / 'active_tm_example.npz');Q=z['Q'];F=z['F'];olds=z['olds'];K=z['K'];Bd=z['Bd'];new=int(z['new']);p=F.shape[0];d=Bd.shape[1];s=len(olds);nu=int(np.sum(eigvalsh(K)<0));V=np.eye(d)[:,olds];Ytrue=F[:,olds]
rng=np.random.default_rng(4419);N=rng.normal(size=Ytrue.shape)+1j*rng.normal(size=Ytrue.shape);r=.02*np.linalg.norm(Ytrue);N=N*r/np.linalg.norm(N);Y=Ytrue+N
Z=np.vstack((V,Y));H=Z.conj().T@Q@Z;H=(H+H.conj().T)/2;t=r*r;Pt=H-t*np.eye(s);kold=int(np.sum(eigvalsh(Pt)>0));Qc=Q-Q@Z@np.linalg.solve(Pt,Z.conj().T@Q);Qc=(Qc+Qc.conj().T)/2;P=Qc[d:,d:];ee=eigvalsh(P);a=np.eye(d)[:,new];f=Qc[d:,:d]@a;mu=-np.linalg.solve(P,f);beta=float(np.real(np.vdot(f,np.linalg.solve(P,f))-np.vdot(a,Qc[:d,:d]@a)));y=F[:,new];q=float(np.real(np.vdot(np.r_[a,y],Qc@np.r_[a,y])));truequad=float(np.real(np.vdot(y-mu,P@(y-mu))));scales=[]
for lam in [1,10,100]:
 Za=np.column_stack((Z,lam*np.r_[a,y]));Ha=Za.conj().T@Q@Za;Ha=(Ha+Ha.conj().T)/2;ev=eigvalsh(Ha)[::-1];tail=float(np.maximum(ev[nu:],0).sum());scales.append({'lambda':lam,'tail':tail,'old_budget_squared':float(t),'conditional_schur':lam*lam*q-t})
rep={'scope':'Floating-point physical TM check of clean-response ellipse with nonzero old noise; no interval certification','nu':nu,'kold':kold,'relative_old_noise':.02,'r_old':float(r),'t_old':float(t),'lambda_min_receiver':float(ee[0]),'beta_clean_no_added_t':beta,'true_response_quad':truequad,'clean_conditional_q':q,'clean_ellipse_max_radius':float(np.sqrt(beta/ee[0])),'future_true_norm':float(np.linalg.norm(y)),'hypothetical_scaling_checks':scales}
open(Path(__file__).resolve().parent / 'active_clean_results.json','w').write(json.dumps(rep,indent=2));print(json.dumps(rep,indent=2))
np.savez(Path(__file__).resolve().parent / 'active_clean_example.npz',Q=Q,Y=Y,V=V,Qc=Qc,mu=mu,beta=beta,r=r,y=y,a=a,K=K,F=F)
