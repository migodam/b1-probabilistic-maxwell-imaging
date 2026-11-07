"""Exposed-state current curvature and reset-step diagnostic; no physical forward/GT."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import json,sys,time
from collections import Counter
import numpy as np
from scipy.linalg import svd,eigh
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
from research.research_max_struct_r6.src.fft_tsom_ops import NativeFourierAPCS
from research.research_max_struct_r10.src.lowrank_denominator import LowRankDenominator

def load(p):
    with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def scalar(x):return float(np.real(x))
def norm2(a,w):return scalar(np.sum(abs(a)**2*w[None,:]))
def main():
    cpu,wall=time.process_time(),time.perf_counter();cost=Counter();rows=[]
    inventory=json.loads((ROOT/'delegated/diagnostics/STATE_INVENTORY.json').read_text())
    out=ROOT/'RESULTS/current_curvature_exposed';out.mkdir(parents=True,exist_ok=True)
    for obj in inventory['objects'][:3]:
        pub={}
        for p in obj['portable_public_files']:pub.update(load(p))
        G,E,S,Y=(pub[k] for k in ('G','E','S','Y'))
        source=load(Path(obj['folder'])/'SOURCE_METRICS.npz');wd,ws=source['data_weights'],source['state_weights']
        basis=load(Path(obj['folder'])/'DENOMINATOR_BASIS.npz')
        _,_,Vh=svd(S,full_matrices=False,check_finite=False);cost['receiver_SVD']+=1
        for st in obj['stages']:
            data=load(st['path']);m,J,Q=(data[k] for k in ('m','J','Q'));meta=st['metadata']
            beta=float(meta['beta_current']);mf=int(meta['basis_info']['mf']);r=meta['receiver_rank']
            U=NativeFourierAPCS((24,24),mf,Vh[:r].conj().T,cost)
            raw=U.apply(np.eye(U.k,dtype=complex));left,sing,_=svd(raw,full_matrices=False,check_finite=False)
            rank=int(np.sum(sing>sing[0]*max(raw.shape)*np.finfo(float).eps));basisU=left[:,:rank]
            cost['actual_APCS_state_basis_SVD']+=1;cost['diagnostic_dense_G_columns']+=rank
            base=basisU-m[:,None]*(G@basisU);SU=S@basisU
            op=LowRankDenominator(G,basis['U'],'galerkin',cost);op.B=basis['B'];factor=op.factor(m,beta)
            e=J-m[:,None]*Q;D=S@J-Y
            H0=beta*base/(1+beta*m[:,None]);HK=beta*op.solve(factor,base)
            R0=beta*e/(1+beta*m[:,None]);RK=beta*op.solve(factor,e)
            arrays=dict(m=m,APCS_singular_values=sing);summaries={};directions={}
            for name,H,R in [('scalar',H0,R0),('green',HK,RK)]:
                dataGram=SU.conj().T@SU;stateGram=H.conj().T@H
                spectra=[]
                for j in range(Y.shape[1]):
                    gram=wd[j]*dataGram+ws[j]*stateGram
                    spectra.append(eigh(gram,eigvals_only=True,check_finite=False))
                    cost['source_current_Hessian_eigh']+=1
                spectrum=np.asarray(spectra);small=scalar(spectrum.min());large=scalar(spectrum.max())
                grad=SU.conj().T@(D*wd)+H.conj().T@(R*ws)
                direction=-grad;den=norm2(SU@direction,wd)+norm2(H@direction,ws)
                alpha=scalar(np.vdot(grad,grad))/den;step=direction*alpha
                directions[name]=step
                summaries[name]=dict(min_eigenvalue=small,max_eigenvalue=large,
                    positive_definite=bool(small>0),condition=large/small if small>0 else None,
                    gradient_norm=float(np.linalg.norm(grad)),reset_line_alpha=alpha,
                    current_step_Frobenius=float(np.linalg.norm(step)),
                    old_objective=.5*(norm2(D,wd)+norm2(R,ws)),
                    new_objective=.5*(norm2(D+SU@step,wd)+norm2(R+H@step,ws)))
                arrays[name+'_spectrum']=spectrum;arrays[name+'_current_step']=basisU@step
            x,y=directions['scalar'],directions['green']
            cosine=scalar(np.vdot(x,y))/(np.linalg.norm(x)*np.linalg.norm(y))
            sid=f"{obj['case_id']}_s{st['stage']}"
            rows.append(dict(case_id=obj['case_id'],stage=st['stage'],beta_current=beta,mf=mf,
                independent_complex_current_DOF_per_source=rank,nominal_Fourier_columns=U.k,
                comparison='Same m/J/Q, same actual APCS state space; orthonormal coordinates for curvature, reset PR diagnostic only',
                scalar=summaries['scalar'],green=summaries['green'],reset_step_cosine=cosine,
                max_curvature_ratio=summaries['green']['max_eigenvalue']/summaries['scalar']['max_eigenvalue'],
                interpretation='No GT or next material readout; orthonormal reset step is not the saved native PR next step'))
            np.savez_compressed(out/(sid+'.npz'),**arrays)
    (out/'METRICS.json').write_text(json.dumps(rows,indent=2)+'\n')
    (out/'RECEIPT.json').write_text(json.dumps(dict(status='COMPLETE',truth_access='NONE',physical_forward=0,physical_adjoint=0,
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,auxiliary_counts=dict(cost),
        scope='Six exposed states; current-Hessian curvature in independent state coordinates, not a proof of recovery'),indent=2)+'\n')
    print(json.dumps([dict(case_id=r['case_id'],stage=r['stage'],ratio=r['max_curvature_ratio'],cos=r['reset_step_cosine'],
        condition_scalar=r['scalar']['condition'],condition_green=r['green']['condition']) for r in rows],indent=2))
if __name__=='__main__':main()
