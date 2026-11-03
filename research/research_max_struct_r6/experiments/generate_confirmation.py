"""R6 data custodian only. Never import from inversion. Write once, no GT exposure."""
from pathlib import Path
from datetime import datetime, timezone
import json, sys, time, uuid
import numpy as np
from scipy.stats import chi2
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from research.research_z_alg_r4_diversity.src.physics import Model

def dump(p,obj):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f: json.dump(obj,f,indent=2,allow_nan=False)

def rel(p): return str(p.relative_to(REPO))

PRIOR=dict(host=.1,amplitude=12.,material_bounds=[.02,20.],smooth_weight=.0001,tv_weight=.0001,tv_epsilon=.01,prior_width=19.5)

def geometry(family,rng):
    # Coordinates in cell widths, integrated 16x16 per cell. Rigid perturbations
    # preserve analytic minimum widths/gaps; no font/template or reduced inversion.
    q=(np.arange(24*16)+.5)/16-12
    y,x=np.meshgrid(q,q,indexing='ij')
    dx,dy=rng.uniform(-.22,.22,2); theta=rng.uniform(-.06,.06)
    X=(x-dx)*np.cos(theta)+(y-dy)*np.sin(theta)
    Y=-(x-dx)*np.sin(theta)+(y-dy)*np.cos(theta)
    s=float(rng.uniform(.985,1.015)); X/=s; Y/=s
    r=np.hypot(X,Y)
    if family=='eccentric_thin_ring':
        inside=(r<8.8)&(np.hypot(X-.55,Y+.25)>5.65)
        width,gap=2.54*s,None
    elif family=='slender_curved_bar':
        a=np.arctan2(Y,X); inside=(r>5.3)&(r<8.1)&(a>-.45)&(a<2.8)
        width,gap=2.8*s,None
    elif family=='near_double_island':
        inside=(np.hypot(X+4.55,Y+.3)<3.25)|(np.hypot(X-4.55,Y-.3)<3.35)
        width,gap=6.5*s,(np.hypot(9.1,.6)-6.6)*s
    elif family=='slit_shell':
        inside=(r>5.35)&(r<8.45)&~((X>0)&(np.abs(Y)<1.05))
        width,gap=3.1*s,2.1*s
    elif family=='two_layer_cavity':
        inside=((r<9.2)&(r>6.6))|((r<4.6)&(r>2.15))
        width,gap=2.45*s,2*s
    elif family=='three_unequal_islands':
        discs=[(-5.,-4.,2.2),(5.,-3.,2.7),(0.,5.,3.)]
        inside=np.zeros_like(X,dtype=bool)
        for cx,cy,rr in discs: inside|=np.hypot(X-cx,Y-cy)<rr
        width=4.4*s
        gap=min(np.hypot(a[0]-b[0],a[1]-b[1])-a[2]-b[2] for i,a in enumerate(discs) for b in discs[i+1:])*s
    else: raise ValueError(family)
    alpha=inside.reshape(24,16,24,16).mean((1,3))
    return alpha,dict(dx_cells=float(dx),dy_cells=float(dy),rotation_radians=float(theta),scale=s,minimum_thickness_cells=float(width),minimum_gap_cells=None if gap is None else float(gap))

def main():
    t,c=time.perf_counter(),time.process_time()
    pub=ROOT/'data_public'; private=ROOT/'sealed_truth'
    if (pub/'cases.json').exists() or (ROOT/'audit/GENERATION_RECEIPT.json').exists(): raise FileExistsError('Preserve write-once attempt')
    records=[]; truth=[]
    metadata=json.loads((REPO/'research/research_fbg_r3_experiments/data/e3_eval_metadata_portable.json').read_text())
    historic=[r for r in metadata if r['amplitude']==12 and r['noise']=='clean']
    assert len(historic)==12 and len({r['object_id'] for r in historic})==12
    for r in historic:
        base=REPO/'research/research_fbg_r2_adaptive_feedback'
        obs=base/'data_contract/stage4_public'/f"{r['case_id']}.npz"
        op=base/'data_contract/stage4_public/operators'/f"{r['layout']}.npz"
        assert obs.exists() and op.exists()
        records.append(dict(case_id=r['case_id'],cohort='DEV_ALL_DISCLOSED',observation_file=rel(obs),operators_file=rel(op),grid_shape=[24,24],**PRIOR))
        truth.append(dict(source='FBG-R3 historical clean amplitude12',**r,truth_path=rel(base/'truth_after_freeze'/r['truth_file'])))
    extras=[('obs_b20d5a55677d4eca','research_max_struct_r2','data_public','sealed_truth'),('obs_42cb233e4f224bb3','research_max_struct_r4','data_public/confirmation','sealed_truth/confirmation'),('obs_63abafa750fc483c','research_max_struct_r4','data_public/confirmation','sealed_truth/confirmation')]
    for ident,study,dp,tp in extras:
        obs=REPO/'research'/study/dp/(ident+'.npz'); tr=REPO/'research'/study/tp/(ident+'.npz')
        assert obs.exists() and tr.exists()
        op=obs if study.endswith('r2') else obs.parent/'operators.npz'
        records.append(dict(case_id=ident,cohort='DEV_ALL_DISCLOSED',observation_file=rel(obs),operators_file=rel(op),grid_shape=[24,24],**PRIOR))
        truth.append(dict(case_id=ident,truth_path=rel(tr),source=study))
    assert len({r['case_id'] for r in records})==15
    dump(private/'scoring/historical_mapping_PRIVATE.json',truth)
    operator=REPO/'research/research_z_alg_r3_digits/public_data/operators.npz'
    with np.load(operator,allow_pickle=False) as z: ops={k:z[k] for k in z.files}
    assert (ops['G'].shape,ops['E0'].shape,ops['S'].shape)==((576,576),(576,12),(48,576))
    model=Model(ops); rng=np.random.default_rng(20261103)
    families=['eccentric_thin_ring','slender_curved_bar','near_double_island','slit_shell','two_layer_cavity','three_unequal_islands']; rng.shuffle(families)
    costs=[]; mapping=[]
    fresh=pub/'confirmation'; fresh.mkdir(parents=True,exist_ok=True)
    private_new=private/'confirmation'; private_new.mkdir(parents=True,exist_ok=True)
    with threadpool_limits(limits=1):
        for family in families:
            seed=int(rng.integers(0,2**63-1)); local=np.random.default_rng(seed)
            ident='obs_'+uuid.uuid4().hex[:16]
            alpha,geom=geometry(family,local)
            assert geom['minimum_thickness_cells']>=2 and (geom['minimum_gap_cells'] is None or geom['minimum_gap_cells']>=1.5)
            m=.1+12*alpha.ravel(); a,b=time.perf_counter(),time.process_time()
            state=model.forward(m,label='R6_NEW_GENERATION'); clean=state['y']
            sigma=.001*np.linalg.norm(clean)/np.sqrt(clean.size)
            noise=sigma/np.sqrt(2)*(local.normal(size=clean.shape)+1j*local.normal(size=clean.shape))
            radius=float(np.sqrt(chi2.ppf(.99,2*clean.size)/2))
            np.savez_compressed(private_new/(ident+'.npz'),m=m,alpha=alpha,mask=alpha>=.5,Y_clean=clean,noise=noise)
            np.savez_compressed(fresh/(ident+'.npz'),Y=clean+noise,sigma=sigma,noise_radius=radius,grid_shape=[24,24],**PRIOR)
            records.append(dict(case_id=ident,cohort='NEW_CONFIRMATION_UNSCORED',observation_file=rel(fresh/(ident+'.npz')),operators_file=rel(operator),grid_shape=[24,24],noise_relative_complex_RMS=.001,**PRIOR))
            mapping.append(dict(case_id=ident,family=family,seed_PRIVATE=seed,geometry_PRIVATE=geom,truth_path=rel(private_new/(ident+'.npz'))))
            costs.append(dict(case_id=ident,LU=1,source_complex_RHS=12,state_relative_residual=state['state_relative_residual'],wall_seconds=time.perf_counter()-a,cpu_seconds=time.process_time()-b))
    dump(private_new/'mapping_PRIVATE.json',dict(master_seed_PRIVATE=20261103,records=mapping,release='Only independent scorer after all predictions sealed; never solver'))
    dump(pub/'cases.json',dict(status='PUBLIC_NO_GT',path_base='repository_root',historical_count=15,new_count=6,records=records,noise_new='CN: E|noise|^2=sigma^2; sigma=.001*norm(Yclean)/sqrt(576)',prior_scope='identical prior for all methods; historical physical observations kept unchanged',inverse_crime='NEW same-grid forward/inverse discretization; no cross-grid claim'))
    dump(ROOT/'audit/GENERATION_RECEIPT.json',dict(status='GENERATED_PRIVATE_TRUTH_HELD',created_utc=datetime.now(timezone.utc).isoformat(),backend='unchanged research.research_z_alg_r4_diversity.src.physics.Model',operator_source=rel(operator),operator_assembly='NONE; existing disk-pulse arrays reused',new_cases=6,historical_cases=15,new_physical_counts=dict(model.counts),records=costs,failed_attempts=[dict(stage='historical_inventory_before_forward',error='TypeError duplicate case_id in dictionary constructor; corrected before any data generated',LU=0,source_complex_RHS=0,process_wall_seconds=.692821625,charged=True)],wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-c,scope='including historical inventory/operator load/six complete forwards/noise/serialization; imports and final receipt write excluded',truth_exposed_to_solver=False,inverse_crime=True))
    print(json.dumps(dict(status='COMPLETE',historical=15,new=6,cost=dict(model.counts),manifest=str(pub/'cases.json'))))

if __name__=='__main__': main()
