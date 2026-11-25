from pathlib import Path
import json
import numpy as np
from scipy.linalg import eigvalsh
z=np.load(Path(__file__).resolve().parent / 'active_tm_example.npz')
G,S,Bd,K,Q,x,F,olds,new,t,C,center,beta=[z[k] for k in ['G','S','Bd','K','Q','x','F','olds','new','t','C','center','beta']]
n=len(x);p=S.shape[0];d=Bd.shape[1];s=len(olds);nu=int(np.sum(eigvalsh(K)<0));V=np.eye(d)[:,olds];Y=F[:,olds];Z=np.vstack((V,Y));aa=np.eye(d)[:,int(new)]
H=(Z.conj().T@Q@Z);H=(H+H.conj().T)/2;P=H-t*np.eye(s);Cyy=C[d:,d:]
ev,U=np.linalg.eigh(Cyy);inverse_sqrt=(U/np.sqrt(ev))@U.conj().T
rng=np.random.default_rng(9182); mismatches=0;checkmax=0
for i in range(500):
 v=rng.normal(size=p)+1j*rng.normal(size=p);v=v/np.linalg.norm(v)
 scale=.5 if i%2 else 1.5; yy=center+scale*np.sqrt(beta)*inverse_sqrt@v
 zz=np.r_[aa,yy]; aug=np.column_stack((Z,zz));Ha=aug.conj().T@Q@aug-t*np.eye(s+1);Ha=(Ha+Ha.conj().T)/2
 count=int(np.sum(eigvalsh(Ha)>1e-9)); schur=float(np.real(np.vdot(zz,C@zz)-t)); predicted=nu+int(schur>0)
 mismatches+=int(count!=predicted)
 ell=float(np.real(np.vdot(yy-center,Cyy@(yy-center))-beta));checkmax=max(checkmax,abs(ell-schur))
# Exact no-information count test for fewer initial sources and batches within rank deficit.
oldsmall=olds[:nu-1];Vs=np.eye(d)[:,oldsmall];Zs=np.vstack((Vs,F[:,oldsmall]));Hs=Zs.conj().T@Q@Zs;kold=int(np.sum(eigvalsh(Hs)>1e-9));slack=nu-kold;vacuous_errors=0
for i in range(100):
 yy=100*(rng.normal(size=p)+1j*rng.normal(size=p));zz=np.r_[aa,yy];Za=np.column_stack((Zs,zz));Ha=Za.conj().T@Q@Za
 if int(np.sum(eigvalsh(Ha)>1e-7))>nu:vacuous_errors+=1
report={'model':'16-cell 2D TM pulse-cell Green; circular self term; 0.7 wavelengths square; 12 receivers half aperture at 2 wavelengths; 20 plane-wave dictionary directions', 'claim_scope':'Floating-point algebra/mechanism check only; no interval-rounding certificate, continuum-error claim, or competing-class separation result', 'box':[0,3], 'sector_weight':1.0,'optical_weight':0.,'n':n,'p':p,'dictionary_size':d,'old_source_indices':olds.tolist(),'new_source_index':int(new),'s':s,'nu':nu,'lambda_min_K':float(eigvalsh(K)[0]),'min_abs_eigenvalue_K':float(np.min(abs(eigvalsh(K)))),'threshold_t':float(t),'joint_noise_radius':float(np.sqrt(t)),'nplus_old_threshold':int(np.sum(eigvalsh(P)>0)),'old_threshold_min_abs_gap':float(np.min(abs(eigvalsh(P)))),'lambda_min_cond_receiver':float(ev[0]),'beta':float(beta),'max_ellipse_radius':float(np.sqrt(beta/ev[0])),'true_new_response_norm':float(np.linalg.norm(F[:,int(new)])),'true_quad':float(np.real(np.vdot(F[:,int(new)]-center,Cyy@(F[:,int(new)]-center)))),'future_ellipse_inertia_checks':500,'future_ellipse_inertia_mismatches':mismatches,'max_square_completion_error':checkmax,'old_small_sources':len(oldsmall),'old_small_positive_index':kold,'old_small_slack':slack,'arbitrary_huge_reading_checks':100,'no_information_mismatches':vacuous_errors}
open(Path(__file__).resolve().parent / 'active_tm_results.json','w').write(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
