from pathlib import Path
import json,sys,shutil,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from material_fiber.core import OUT,OLD,CAPS,write,save
from p1_material.prior import MaterialPrior

def main():
    if (OUT/'configs/AUTHORIZATION.json').exists():raise RuntimeError('campaign already initialized; no replacement')
    cfg=json.loads((ROOT/'research/p1_material_target_20261018/configs/ONLINE_TARGET.json').read_text())
    write(OUT/'configs/TARGET.json',cfg)
    write(OUT/'configs/AUTHORIZATION.json',dict(campaign_id='P2-material-fiber-20261019',explicit_user_authorization=True,user_instruction='PLEASE IMPLEMENT THIS PLAN; independent 4000 states /132000 complex RHS /8h CPU physics; GPU training <=8h',target_id=cfg['target_id'],caps={k:dict(states=v[0],complex_RHS=v[1],controlled_wall_seconds=v[2]) for k,v in CAPS.items()},total=dict(states=4000,complex_RHS=132000,controlled_wall_seconds=28800,GPU_training_seconds=28800),old_balances_used=False,new_Mie=0,new_3D=0))
    points=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points'];P=MaterialPrior(cfg,points)
    (OUT/'inputs/background').mkdir(parents=True,exist_ok=True);shutil.copy2(OLD/'inputs/background/weak.npz',OUT/'inputs/background/weak.npz')
    rows=[]
    for split,n,code in [('train',12,0),('validation',4,1),('heldout',8,2)]:
        for i in range(n):
            parent=f'P2_{split}_{i:03d}';rng=np.random.default_rng(np.random.SeedSequence([2026101901,code,i,0]));q=P.draw(rng)[0]
            nr=np.random.default_rng(np.random.SeedSequence([2026101901,code,i,1]));noise=cfg['noise']['sigma']/np.sqrt(2)*(nr.normal(size=(32,8))+1j*nr.normal(size=(32,8)))
            save(OUT/f'inputs/offline_only/{parent}.npz',q=q,noise=noise)
            rows.append(dict(parent=parent,split=split,index=i,seed_derivation=[2026101901,code,i],GT_online=False))
    write(OUT/'configs/PARENTS.json',rows)
    write(OUT/'configs/PROTOCOL.json',dict(campaign='P2-material-fiber-20261019',target_id=cfg['target_id'],master_seed=2026101901,fresh_splits={'train':12,'validation':4,'heldout':8},phase_A_methods=['fixed_LM','adaptive_TR','material_TR','coarse_TR'],phase_B_methods=['random_matched','native','hard_prior','soft_Maxwell','soft_frozen','soft_dynamic','soft_adaptive','full_J'],radii=[.02,.04],stress_radius=.08,material_svd_threshold=1e-8,soft_beta=1,soft_p=1,lambda_cutoff='1/s^2',radius_relative_tolerance=.02,signed_progress_fraction=.5,likelihood_competitive_DeltaPhi=2,posterior_mass=None,methods_frozen_before_heldout_scoring=True,controlled_states_include_rejections=True,source_commits={'B1':'8e4519e687acfc5de7782c7b7ca0d055b013929f','A30':'717d3d26a290248bd2b9aeb55f3d03694b407b84','C1':'83056881e89329212b16f5be47a4d0dbc4465943'},historical_assets_read_only=True))
    # Reuse paid, non-GT B endpoint and its paid receiver transfer.
    with np.load(OLD/'results/physical/endpoints/state_0048.npz') as z:d={k:z[k].copy() for k in z.files}
    found=None
    for p in (OLD/'results/physical/jacobians').glob('*.npz'):
        with np.load(p) as z:
            if np.array_equal(z['q'],d['q']):found=z['T'].copy();break
    if found is None:raise RuntimeError('paid B Jacobian missing')
    save(OUT/'cache/OLD_B.npz',**d,receiver_T=found)
    (OUT/'inputs/online').mkdir(parents=True,exist_ok=True);shutil.copy2(OLD/'inputs/online/P1_obs_B.npz',OUT/'inputs/online/P1_obs_B.npz')
    write(OUT/'SOURCE_INVENTORY.json',dict(baseline='8e4519e',old_paid_states=26,old_paid_RHS=976,old_ledgers_modified=False,background=str(OLD/'inputs/background/weak.npz'),old_B_endpoint=str(OLD/'results/physical/endpoints/state_0048.npz'),new_physics_in_preparation=0,Semantic_Scholar_status='2026-10-19 focused collect returned HTTP429; no citations or novelty conclusions inferred',agent_status='native gpt-6.1-sol medium spawn attempted; existing thread slots full; parent continues'))
    print(json.dumps(dict(parents=len(rows),new_physics=0,output=str(OUT))))
if __name__=='__main__':main()
