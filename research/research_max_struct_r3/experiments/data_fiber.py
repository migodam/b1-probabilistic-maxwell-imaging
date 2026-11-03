"""One inherited nonpermutation fiber pair, one prescribed physical loss.
No re-solving or re-certifying the old atlas; new lossy response and tangent only.
"""
from pathlib import Path
import sys,json,time
from datetime import datetime,timezone
import numpy as np
from scipy.linalg import lu_solve,null_space
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_z_alg_r4_diversity.src.physics import build_operators
from research.research_max_struct_r3.src.passive_continuation import PassiveModel
from research.research_max_struct_r2.experiments.run_blind import clean

def main():
    t,cpu=time.perf_counter(),time.process_time();old=REPO/'research/research_z_alg_asymmetric_validation'
    cfg=json.loads((old/'EXPERIMENT_CONTRACT.json').read_text());geo=cfg['E1']
    node='results/E1/P1/B/RADIAL_2_08/ITERATIONS.json'
    m1=np.array(json.loads((old/node).read_text())[-1]['m']);m0=np.array(geo['public_reference_materials'][1])
    ops=build_operators(geo['centers'],geo['sources'],geo['receivers'],cfg['disk_radius'],cfg['k']);ops['grid_shape']=[2,3]
    states={};models=[];jac={};out=ROOT/'results/fiber';out.mkdir(parents=True,exist_ok=True)
    if(out/'RECEIPT.json').exists():raise FileExistsError('existing evidence')
    for gamma in [0.,1.]:
        model=PassiveModel(ops,gamma);models.append(model)
        states[(gamma,1)]=model.state(m1)
        states[(gamma,0)]=st=model.state(m0)
        with model._physics.operation('RECEIVER_TRANSFER',rhs=model.nr,receiver_transpose_RHS=model.nr):
            O=lu_solve(st['lu'],model.S.T,trans=1,check_finite=False).T
        J=np.einsum('ri,is->sri',O,st['Q']).reshape(-1,model.n)
        jac[gamma]=np.stack((J.real,J.imag),axis=1).reshape(-1,model.n)
    N=null_space(jac[0.]);u,s,vh=np.linalg.svd(jac[1.]@N,full_matrices=False);h=N@vh[0]
    gaps={str(g):float(np.linalg.norm(states[(g,1)]['Y']-states[(g,0)]['Y'])) for g in [0.,1.]}
    keys=set().union(*(m.cost for m in models));counts={k:sum(m.cost.get(k,0)for m in models)for k in keys}
    rec=dict(status='COMPLETED',classification='INHERITED_PUBLIC_ZALG_PAIR_NEW_LOSS_DIAGNOSTIC',started_source=node,
        time_convention='exp(-i omega t)',material_class='real[.5,4] + common i gamma',gamma=[0,1],mesh='six nonoverlapping radius.04 disk pulses',sources=2,receivers=1,
        old_pair_gap=gaps['0.0'],new_pair_gap=gaps['1.0'],
        pair_lipschitz_ratio=gaps['1.0']/gaps['0.0'],material_rms=float(np.linalg.norm(m1-m0)/np.sqrt(6)),
        old_real_jacobian_singular_values=np.linalg.svd(jac[0.],compute_uv=False).tolist(),
        old_tangent_norm=float(np.linalg.norm(jac[0.]@h)),new_tangent_norm=float(np.linalg.norm(jac[1.]@h)),
        max_state_defect=max(st['state_relative_residual']for st in states.values()),
        interpretation='EMPIRICAL finite model stable-transform obstruction; near equality is not proof against arbitrary discontinuous transform; theorem has separate exact-fiber assumptions',
        truth_access='DISCLOSED_EXISTENCE_ANCHOR_NOT_IMAGING',physical_cost=counts,Green_setups=1,Green_kernel_evaluations=54,
        wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-cpu,finished_utc=datetime.now(timezone.utc).isoformat(),
        old_atlas_revalidation=False,code_hash=None,config_hash=None)
    np.savez_compressed(out/'arrays.npz',G=ops['G'],E0=ops['E0'],S=ops['S'],m0=m0,m1=m1,J0=jac[0.],J1=jac[1.],null_tangent=h,
        Y00=states[(0.,0)]['Y'],Y01=states[(0.,1)]['Y'],Y10=states[(1.,0)]['Y'],Y11=states[(1.,1)]['Y'])
    (out/'RECEIPT.json').write_text(json.dumps(clean(rec),indent=2));print(json.dumps(rec))
if __name__=='__main__':
    with threadpool_limits(limits=1):main()
