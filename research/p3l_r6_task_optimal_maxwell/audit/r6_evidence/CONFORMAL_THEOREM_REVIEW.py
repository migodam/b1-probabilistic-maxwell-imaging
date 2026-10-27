"""Saved-vector-only review of THEORY §2. Imports no Target/solver/trainer."""
from pathlib import Path
import json,time
import numpy as np
ROOT=Path.cwd();OUT=ROOT/'audit/r6_evidence';start=time.perf_counter()
with np.load(ROOT/'track_a/raw/COMPENSATION_VECTORS.npz',allow_pickle=False) as z:J=z['J_material_white'].copy();c=z['c_white'].copy();q=z['q'].copy()
with np.load(ROOT/'track_a/raw/COMPENSATION_BASE_CACHE.npz',allow_pickle=False) as z:baseJ=z['J'].copy();baseq=z['q'].copy()
assert np.array_equal(q,baseq);assert J.shape==(512,2) and c.shape==(512,)
assert np.array_equal(J[:,0],baseJ[:,64]) and np.array_equal(J[:,1],baseJ[:,65]/.15)
def rotate(v):
 x=v.reshape(8,64);return np.concatenate((-x[:,32:],x[:,:32]),axis=1).ravel()
cr=float(np.linalg.norm(J[:,1]-rotate(J[:,0]))/np.linalg.norm(J[:,0]));assert cr<1e-12
D=complex(1.2,.0975);kappa=2.74625;rebuilt=-kappa*(D.real*J[:,0]+D.imag*J[:,1]);cerr=float(np.linalg.norm(c-rebuilt)/np.linalg.norm(c));assert cerr<1e-12

def measure(A,coef):
 gram=A.T@A;scale=np.trace(gram)/2;u,_=np.linalg.qr(A);one=A[:,0]/np.linalg.norm(A[:,0]);return dict(gram=gram.tolist(),gram_relative_deviation=float(np.linalg.norm(gram-scale*np.eye(2))/scale),column_norm_squared=float(scale),rank=int(np.linalg.matrix_rank(A)),two_property_normal_ratio=float(np.linalg.norm(coef-u@(u.T@coef))/np.linalg.norm(coef)),one_real_property_normal_ratio=float(np.linalg.norm(coef-one*(one@coef))/np.linalg.norm(coef)))
def select(sources,receivers):return np.array([s*64+offset+r for s in sources for offset in (0,32) for r in receivers])
rows=[]
for sources,receivers,label in [(range(n),range(32),f'first{n}sources_allreceivers') for n in (1,2,4,8)]+[(range(8),range(0,32,2),'all8_even_receivers'),([0,2],[1,7,19],'sources02_receivers1719')]:
 ids=select(sources,receivers);m=measure(J[ids],c[ids]);m.update(label=label,real_components=len(ids),complex_channels=len(ids)//2);assert m['gram_relative_deviation']<1e-12 and m['two_property_normal_ratio']<1e-12 and abs(m['one_real_property_normal_ratio']-abs(D.imag)/abs(D))<1e-12;rows.append(m)
# Deliberately not complex-compatible: double every real component only.
weights=np.tile(np.r_[np.full(32,2.),np.ones(32)],8);improper=measure(weights[:,None]*J,weights*c)
# Complex-compatible unequal channel precision preserves conformality.
channel=np.linspace(.3,1.7,256).reshape(8,32);sameweights=np.c_[channel,channel].ravel();proper=measure(sameweights[:,None]*J,sameweights*c);assert proper['gram_relative_deviation']<1e-12
result=dict(status='SAVED_VECTOR_ALGEBRA_PASS_CONDITIONAL_THEOREM_SCOPE',CR_rotation_relative_error=cr,coefficient_reconstruction_relative_error=cerr,subset_checks=rows,theoretical_one_property_ratio=abs(D.imag)/abs(D),improper_real_imag_scaling_diagnostic=improper,proper_channel_weighting_diagnostic=proper,
 qualifications=['Read only already-paid COMPENSATION_VECTORS and BASE_CACHE; original vector acquisition40RHS remains historical paid evidence.','Conformal Gram and compensation residual follow mathematically from the same holomorphic derivative; their repetition is not independent physical evidence or branch coverage.','Vector-valued holomorphic F is locally conformal as a real2 immersion only when a!=0; not a global injectivity/distance-preserving map.','Complex-compatible proper whitening can be correlated/nonisotropic; isotropic noise is sufficient, not necessary.','Interior propertybox required for exact D/Q compensation; boundaries can block the compensating path.'],new_Maxwell_RHS=0,new_training_updates=0,wall_seconds=time.perf_counter()-start)
(OUT/'CONFORMAL_THEOREM_REVIEW.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
