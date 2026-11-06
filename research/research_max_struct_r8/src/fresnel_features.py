"""Frozen ten-feature NN transform for source-specific measured receiver maps.
No labels, nominal geometry, calibration gain fitting or heldout scattering.
"""
import time
import numpy as np
from scipy.linalg import svd
from research.research_max_struct_r8.src.warm_features import FEATURE_NAMES

def features(model,Y,centers,halfwidth=.08):
 wall,cpu=time.perf_counter(),time.process_time();E,S=model.E,model.Ss;n,ns=E.shape;J=np.zeros_like(E);ranks=[]
 num=np.sum(E.conj()*model.receive_adjoint(Y),axis=1)
 den=np.sum(abs(E)**2*np.sum(abs(S)**2,axis=0),axis=1);born=num/np.maximum(den,1e-300)
 for s in range(ns):
  u,q,vh=svd(np.ascontiguousarray(S[:,:,s]),full_matrices=False,check_finite=False);rank=int(np.count_nonzero(q>1e-3*q[0]));ranks.append(rank)
  J[:,s]=vh[:rank].conj().T@((u[:,:rank].conj().T@Y[:,s])/q[:rank])
 Q=E+model.G@J;read=np.sum(Q.conj()*J,axis=1)/np.maximum(np.sum(abs(Q)**2,axis=1),1e-20)
 energy=np.mean(abs(E)**2,axis=1);je=np.mean(abs(J)**2,axis=1)/np.maximum(energy,1e-20);qe=np.mean(abs(Q)**2,axis=1)/np.maximum(energy,1e-20)
 x=[np.arcsinh(born.real),np.arcsinh(born.imag),np.arcsinh(read.real),np.arcsinh(read.imag),np.log1p(je),np.log(np.maximum(qe,1e-12)),np.ones(n),np.full(n,ns/12),centers[:,0]/halfwidth,centers[:,1]/halfwidth]
 return np.stack(x).reshape(10,*model.shape).astype(np.float32),dict(cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,receiver_SVD=ns,dpcs_ranks=ranks,Green_actions_columns=ns,receiver_adjoint_columns=ns,feature_names=FEATURE_NAMES,coordinates_halfwidth_m=halfwidth,source_aperture_fraction=1.,truth_access='NONE',heldout_scattering_access='NONE',qualification='OUT_OF_DISTRIBUTION measured scalar-TM transfer; frozen synthetic weights')
