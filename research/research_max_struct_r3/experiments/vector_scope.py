"""New scope question: rank-one tensor dyad versus rank-three isotropic cell.
Uses the qualified R1 finite DDA; not scalar TM relabelled as vector Maxwell.
"""
from pathlib import Path
import sys,json,time
from datetime import datetime,timezone
import numpy as np
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r1.src.vector_backend import TensorDDA

def main():
    tick,cpu=time.perf_counter(),time.process_time();out=ROOT/'results/vector_scope';out.mkdir(parents=True,exist_ok=True)
    if(out/'RECEIPT.json').exists():raise FileExistsError('preserve scope run')
    points=np.stack(np.meshgrid([-.03,.03],[-.03,.03],[-.03,.03],indexing='ij'),axis=-1).reshape(-1,3)
    backend=TensorDDA(points,.06**3);eps=np.broadcast_to(np.array([[4.,.3,0],[.3,2.5,.2],[0,.2,3.2]]),(8,3,3)).copy()
    dirs=np.array([[0,0,1],[1,0,0],[0,1,0]]);pol=np.array([[1,0,0],[0,1,0],[0,0,1]])
    recv=.6*np.vstack((np.eye(3),-np.eye(3)));S=backend.external_electric_operator(recv)
    v=np.array([1.,.4,.2]);v/=np.linalg.norm(v);rows=[];preds={};fields={};spectra={}
    base=backend.solve(eps,dirs,pol);Y0=S@base.j
    fields['base']=dict(E=base.E,H=base.H,j=base.j,epsilon=eps);rows.append(dict(family='base',t=0,solver=base.receipt))
    for family,pert in [('dyad',np.outer(v,v)),('isotropic',np.eye(3))]:
        changes=[]
        for t in [1.,2.,4.]:
            material=eps.copy();material[0]+=t*pert;state=backend.solve(material,dirs,pol);Y=S@state.j
            changes.append((Y-Y0).ravel());preds[f'{family}_{int(t)}']=Y
            fields[f'{family}_{int(t)}']=dict(E=state.E,H=state.H,j=state.j,epsilon=material)
            rows.append(dict(family=family,t=t,min_epsilon_eigenvalue=float(np.linalg.eigvalsh(material).min()),max_epsilon_eigenvalue=float(np.linalg.eigvalsh(material).max()),
                branch_material_frobenius=float(np.linalg.norm(material-eps)),port_gap=float(np.linalg.norm(Y-Y0)),solver=state.receipt))
        mat=np.stack(changes,axis=1);s=np.linalg.svd(mat,compute_uv=False);spectra[family]=dict(singular_values=s.tolist(),relative_to_first=(s/s[0]).tolist())
    payload=dict(points=points,receivers=recv,directions=dirs,polarizations=pol,v=v,Y0=Y0)
    payload.update({name+'_'+k:value for name,d in fields.items()for k,value in d.items()});payload.update(preds)
    np.savez_compressed(out/'arrays.npz',**payload)
    rec=dict(status='COMPLETED',classification='FINITE_3D_VECTOR_SCOPE_CONTROL_NOT_IMAGING',physical_convention='exp(-i omega t),lambda0=mu=backgroundeps=1',
        geometry='8cubiccells centers±.03, side.06, totalbodywidth.12lambda',mesh=[2,2,2],sources=3,receivers=6,complex_observations=54,
        backend='unchanged R1TensorDDA CM+RR dyadicGreen',truth_access='DISCLOSED_ANCHOR',prior=None,noise_model='none',
        scope='rank-one material tensor perturbation can have scalar Mobius response; isotropic3D cell generally rank3, scalarTM theorem not automatic',
        physical_cost=dict(LU=7,actual_complex_RHS=21,source_RHS=21),auxiliary_cost=dict(vector_Goff_assembly=1,external_electric_Green_build=1,near_electric_actions=7,
            internal_magnetic_Green_build=7,internal_magnetic_actions=7,response_SVD=2,response_SVD_shape=[54,3]),
        vector_setup_receipt=backend.setup_receipt,states=rows,response_change_spectra=spectra,
        continuum_accuracy='NOT_CLAIMED; no new imaging or grid-convergence advantage; inherited R1Mie qualification distinct from this small finiteDDA scope',
        wall_seconds=time.perf_counter()-tick,cpu_seconds=time.process_time()-cpu,finished_utc=datetime.now(timezone.utc).isoformat(),
        expected_rank_gate=None,gate_reason='normalized singulars reported, no arbitrary rank acceptance threshold')
    (out/'RECEIPT.json').write_text(json.dumps(rec,indent=2));print(json.dumps(spectra))
if __name__=='__main__':
    with threadpool_limits(limits=1):main()
