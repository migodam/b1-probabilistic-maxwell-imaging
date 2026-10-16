"""Reproduce saved-array KKT/box diagnostics without physics imports or LU."""
from pathlib import Path
import sys,json,time
import numpy as np
from scipy.optimize import lsq_linear
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from research_f.common import OUT,Prior,whiten,response_jac,dump

def main():
    start=time.perf_counter();records=[]
    for row in json.loads((OUT/'results/F1/SUMMARY.json').read_text()):
        case=row['case'];p=np.array(row['final_p']);meta=json.loads((OUT/f'inputs/online/{case}.json').read_text());P=Prior(meta)
        with np.load(ROOT/row['final_array']) as z:field=z['field'];prediction=z['prediction']
        with np.load((ROOT/row['final_array']).with_name('state_11_JAC.npz')) as z:A=z['A'];D=z['D'];T=z['transfer'];C=z['C']
        with np.load(OUT/f'inputs/online/{case}.npz') as z:y=z['observation'];sigma=float(z['sigma'])
        r=whiten(prediction-y,sigma);g=A.T@r;step=-np.linalg.pinv(A,rcond=1e-7)@r
        DM=np.r_[D.real,D.imag]/24;U,s,V=np.linalg.svd(DM,full_matrices=False);E=V[:7].T/s[:7];B=A@E
        fit=lsq_linear(A,-r,bounds=(1e-10-p,1-1e-10-p),tol=1e-13,lsq_solver='exact',max_iter=1000)
        v=np.random.default_rng(261017).normal(size=8);direct=whiten(T@(field*(D@v)[:,None]),sigma)
        Ccheck=np.einsum('qp,ps->qsp',T,field).reshape(256,576)
        phys=P.chart.lo+(P.chart.hi-P.chart.lo)*p
        rec=dict(case=case,new_physics=0,new_RHS=0,physical_weights=phys[6:],weight_ratio=phys[6]/phys[7],
            gradient=g,TRF_scaled_gradient_inf=float(np.max(abs(g*np.where(g<0,1-1e-10-p,p-1e-10)))),
            unconstrained_parameter_step_norm=float(np.linalg.norm(step)),unconstrained_material_step_RMS=float(np.linalg.norm(DM@step)),
            material_frame_GN_norm=float(np.linalg.norm(np.linalg.pinv(B,rcond=1e-7)@r)),
            bounded_linear_material_step_RMS=float(np.linalg.norm(DM@fit.x)),bounded_linear_residual=float(np.linalg.norm(r+A@fit.x)),
            bounded_linear_reduction=float(np.linalg.norm(r)-np.linalg.norm(r+A@fit.x)),bounded_linear_optimality=float(fit.optimality),
            C_reconstruction_relative=float(np.linalg.norm(C-Ccheck)/np.linalg.norm(C)),Jv_contraction_relative=float(np.linalg.norm(direct-A@v)/np.linalg.norm(direct)),
            actual_nonlinear_endpoint_evaluated=False,changes_frozen_gate=False)
        records.append(rec)
    dump(OUT/'results/STATIONARITY_REPLAY.json',dict(role='saved linear algebra diagnostic, never qualification revision',records=records,new_states=0,new_RHS=0,wall_seconds=time.perf_counter()-start))
    print(json.dumps([dict(case=r['case'],material_step=r['bounded_linear_material_step_RMS'],linear_gain=r['bounded_linear_reduction']) for r in records]))
if __name__=='__main__':main()
