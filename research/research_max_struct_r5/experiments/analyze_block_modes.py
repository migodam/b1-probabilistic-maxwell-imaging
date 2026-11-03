"""Geometry-fixed modal audit of saved residuals; no forward, optimizer or truth."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[name]='1'
from pathlib import Path
import json, csv, time
from datetime import datetime, timezone
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]

def main():
    wall,cpu=time.perf_counter(),time.process_time()
    config=json.loads((ROOT/'configs/trajectory_replay_frozen.json').read_text())
    # Rows in the fixed block ordering [top-left, top-right, bottom-left, bottom-right].
    H=np.array([[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]],float)/2
    rows=[]; operator_rows=[]
    for run in config['runs']:
        public={}
        for filename in run['public_files']:
            with np.load(REPO/filename,allow_pickle=False) as z:
                public.update({k:z[k] for k in z.files})
        G=public['G']; shape=tuple(public['shape']); grid=np.arange(len(G)).reshape(shape)
        ids=np.array([grid[i:i+2,j:j+2].ravel() for i in range(0,shape[0],2) for j in range(0,shape[1],2)])
        B=G[ids[:,:,None],ids[:,None,:]].copy(); B[:,np.arange(4),np.arange(4)]=0
        # Uniform-background equal-cell scalar Green depends only on separation.
        ge,gd=G[ids[0,0],ids[0,1]],G[ids[0,0],ids[0,3]]
        lam=np.array([2*ge+gd,-gd,-gd,-2*ge+gd])
        if run['method']=='global':
            for beta in (1.,.2):
                for m in (float(public['host']),float(public['material_bounds'][1])):
                    a=m/(1+(beta-G[0,0])*m)
                    operator_rows.append(dict(case=run['case_id'],beta=beta,material=m,
                        eigenvalues=[[x.real,x.imag] for x in lam],
                        plus_weights=np.abs(1+a*lam)**2,minus_weights=np.abs(1-a*lam)**2))
        for path in sorted((ROOT/'RESULTS/trajectory_diagnostics'/run['run_id']).glob('checkpoint*.npz')):
            with np.load(path,allow_pickle=False) as z:
                m,u,v,ws=z['m_before'],z['u'],z['Lu'],z['ws']
            meta=json.loads(path.with_suffix('.json').read_text());beta=meta['beta_material']
            a=m/(1+(beta-np.diag(G))*m)
            um=np.einsum('ab,nbs->nas',H,u[ids])
            vm=np.einsum('ab,nbs->nas',H,v[ids])
            energy=np.sum(np.abs(um)**2*ws[None,None,:],axis=(0,2))
            cross=2*np.real(np.sum(um.conj()*vm*ws[None,None,:],axis=(0,2)))
            # Actual heterogeneous-material block in modal coordinates, not a uniform surrogate.
            Lmodal=np.einsum('ab,nbc,cd->nad',H,a[ids][:,:,None]*B,H.T)
            diagonal=np.diagonal(Lmodal,axis1=1,axis2=2)
            diag_v=diagonal[:,:,None]*um
            diag_cross=2*np.real(np.sum(um.conj()*diag_v*ws[None,None,:],axis=(0,2)))
            off=Lmodal.copy();off[:,np.arange(4),np.arange(4)]=0
            rows.append(dict(run_id=run['run_id'],case=run['case_id'],method=run['method'],checkpoint=path.stem,
                stage=meta['stage'],iteration=meta['iteration'],beta=beta,
                coherent_energy_fraction=float(energy[0]/energy.sum()),
                mode_energy_fractions=(energy/energy.sum()).tolist(),mode_cross=cross.tolist(),
                coherent_cross=float(cross[0]),other_modes_cross=float(cross[1:].sum()),
                diagonal_modal_cross=float(diag_cross.sum()),offdiagonal_modal_cross=float((cross-diag_cross).sum()),
                modal_mixing_frobenius_fraction=float(np.linalg.norm(off)/np.linalg.norm(Lmodal)),
                total_cross=float(cross.sum()),truth_access='NONE'))
    out=ROOT/'RESULTS/block_modes';out.mkdir(exist_ok=True)
    convert=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item() if isinstance(x,np.generic) else str(x)
    payload=dict(classification='POST_SCORE_DEVELOPMENT_NO_TRUTH_INPUT',
        hypothesis='Fixed square-block Walsh modes separate coherent from alternating residuals; homogeneous material is diagonal, heterogeneous material mixes them.',
        mode_order=['constant','x alternating','y alternating','checkerboard'],
        scientific_scope='Residual metric only, not physical Jacobian spectrum, adaptive source choice, or material recovery predictor.',
        rows=rows,uniform_material_reference=operator_rows)
    (out/'MODAL_ANALYSIS.json').write_text(json.dumps(payload,indent=2,default=convert))
    with (out/'MODAL_ANALYSIS.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    receipt=dict(status='COMPLETED',timestamp=datetime.now(timezone.utc).isoformat(),truth_access='NONE',
        new_physical_LU=0,new_complex_RHS=0,optimizer_calls=0,identity_error_floor_tests=0,
        wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,
        diagnostic_checkpoints=len(rows),block_size=4,blocks_per_checkpoint=144,
        cost_scope='read saved arrays and fixed 4x4 matrix products; no forward; import startup outside timer',
        money=None)
    (ROOT/'audit/BLOCK_MODES_RECEIPT.json').write_text(json.dumps(receipt,indent=2))
    print(json.dumps(receipt))

if __name__=='__main__':main()
