"""Write-once private truth + allowlisted public observations. No inversion.

Run only after configs/blind_cases_frozen.json is written and parent authorizes.
No material/shape/seed values are printed. Do not import into an inverse runner.
"""
from pathlib import Path
import sys,json,time,secrets,uuid,traceback
from datetime import datetime,timezone
import numpy as np
from scipy.stats import chi2
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from research.research_z_alg_r4_diversity.src.physics import Model

def utc():return datetime.now(timezone.utc).isoformat()

def dump_new(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x') as f:json.dump(value,f,indent=2,allow_nan=False)

def shape_alpha(family,rng):
    n=192;y,x=np.mgrid[:n,:n];x=(x+.5)/n*2-1;y=(y+.5)/n*2-1
    angle=rng.uniform(-.25,.25);c,s=np.cos(angle),np.sin(angle)
    u,v=c*x-s*y-rng.uniform(-.07,.07),s*x+c*y-rng.uniform(-.07,.07)
    if family=='ring':
        a,b=rng.uniform(.55,.71),rng.uniform(.45,.62)
        outer=(u/a)**2+(v/b)**2<=1
        inner=((u-rng.uniform(-.12,.12))/(a*rng.uniform(.42,.58)))**2+((v-rng.uniform(-.12,.12))/(b*rng.uniform(.42,.58)))**2<1
        support=outer&~inner
    elif family=='curved_bar':
        radius=np.sqrt((u-.03)**2+(v+.01)**2);phi=np.mod(np.arctan2(v+.01,u-.03),2*np.pi)
        support=(abs(radius-rng.uniform(.49,.58))<rng.uniform(.10,.15))&(phi>.25)&(phi<5.25)
    elif family=='two_islands':
        left=((u+.38)/rng.uniform(.20,.29))**2+((v+.12)/rng.uniform(.27,.38))**2<=1
        right=((u-.34)/rng.uniform(.23,.34))**2+((v-.21)/rng.uniform(.17,.29))**2<=1
        support=left|right
    else:raise ValueError(family)
    return support.reshape(24,8,24,8).mean(axis=(1,3))

def generate():
    started,cpu=time.perf_counter(),time.process_time()
    config_path=ROOT/'configs/blind_cases_frozen.json';cfg=json.loads(config_path.read_text())
    if cfg['status']!='FROZEN_BEFORE_GENERATION':raise RuntimeError('Frozen config required')
    public,private=ROOT/'data_public',ROOT/'sealed_truth'
    if any(public.glob('obs_*.npz')) or (public/'manifest.json').exists() or any(private.glob('*.npz')) or (private/'mapping_PRIVATE.json').exists():
        raise FileExistsError('Existing realization must be preserved; generator is write-once')
    public.mkdir(parents=True,exist_ok=True);private.mkdir(parents=True,exist_ok=True)
    receipt=dict(run_id='BLIND_GENERATION',status='RUNNING',started_utc=utc(),
      classification='FRESH_BLIND_GENERATOR_ONLY',code_hash=None,config_hash=None,
      hash_reason='No hash checks per user; source and frozen config paths retained',
      source='experiments/generate_blind.py',frozen_config='configs/blind_cases_frozen.json',
      operator_source=cfg['operator_source'],truth_access='GENERATOR_ONLY_NOT_MAIN_OR_INVERSION',
      failed_attempts=[],operators_reassembled=False,historical_assembly_seconds=None,
      historical_assembly_reason='Frozen public operators loaded; historical assembly not rerun or priced',
      records=[])
    model=None
    try:
        with np.load(REPO/cfg['operator_source'],allow_pickle=False) as data:
            ops={k:data[k] for k in ('G','E0','S','centers','sources','receivers','k','radius')}
        assert (ops['G'].shape,ops['E0'].shape,ops['S'].shape)==((576,576),(576,12),(48,576))
        seed=secrets.randbits(128);rng=np.random.default_rng(seed)
        families=list(cfg['private_shape_families']);rng.shuffle(families)
        model=Model(ops);pending=[];mapping=[]
        with threadpool_limits(limits=1):
            for family in families:
                ident='obs_'+uuid.uuid4().hex[:16]
                alpha=shape_alpha(family,rng);material=.1+12*alpha.ravel()
                tick,ct=time.perf_counter(),time.process_time()
                state=model.forward(material,label='fresh_blind_generation')
                clean=state['y'];M=clean.size;sigma=1e-3*np.linalg.norm(clean)/np.sqrt(M)
                noise=sigma/np.sqrt(2)*(rng.normal(size=clean.shape)+1j*rng.normal(size=clean.shape))
                Y=clean+noise;radius=float(np.sqrt(chi2.ppf(.99,2*M)/2))
                np.savez_compressed(private/f'{ident}.npz',m=material,alpha=alpha,mask=alpha>=.5,Y_clean=clean,noise=noise)
                payload=dict(G=ops['G'],E0=ops['E0'],S=ops['S'],shape=np.array([24,24]),grid_shape=np.array([24,24]),
                  coords=ops['centers'],sources=ops['sources'],receivers=ops['receivers'],Y=Y,
                  k=ops['k'],radius=ops['radius'],sigma=np.array(sigma),noise_radius=np.array(radius),
                  material_bounds=np.array([.02,20.]),amplitude=np.array(12.),host=np.array(.1),
                  smooth_weight=np.array(1e-4),tv_weight=np.array(1e-4),tv_epsilon=np.array(.01),prior_width=np.array(19.5))
                pending.append((ident,payload));mapping.append(dict(observation_id=ident,family=family,truth_file=ident+'.npz',amplitude=12.,host=.1))
                receipt['records'].append(dict(observation_id=ident,physical_LU=1,forward_complex_RHS=12,
                   linear_relative_residual=state['state_relative_residual'],wall_seconds=time.perf_counter()-tick,cpu_seconds=time.process_time()-ct,
                   timer_scope='forward,noise,private npz output; public output counted only in process timer'))
        for ident,payload in sorted(pending):np.savez_compressed(public/f'{ident}.npz',**payload)
        dump_new(private/'mapping_PRIVATE.json',dict(seed_PRIVATE=str(seed),records=mapping,
          policy='Open only after all frozen prediction dispositions sealed and parent scoring authorization; do not upload to inverse process'))
        dump_new(public/'manifest.json',dict(status='FRESH_UNSCORED',created_utc=utc(),
          records=[dict(observation_id=ident,path=ident+'.npz')for ident,_ in sorted(pending)],count=3,
          shape=[24,24],material_unknowns=576,source_count=12,receiver_count=48,
          public_payload_keys=sorted(pending[0][1]),noise='Circular complex Gaussian; real/imag sigma^2/2; sigma and99%whitened radius in each file',
          objective='Normalized physical data + unchanged historical prior',truth_access='NONE_TO_INVERSION',
          family_mapping='Private until all predictions sealed; all three high amplitude12',frozen_config='configs/blind_cases_frozen.json'))
        receipt.update(status='GENERATED_PRIVATE_TRUTH_HELD',finished_utc=utc(),observations=3,
          physical_counts=dict(model.counts),wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu,
          timer_scope='Inside generate after imports, includes operator load,3 forward states,noise and public/private serialization; receipt own write excluded')
    except Exception as exc:
        receipt.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc(),
          physical_counts=dict(model.counts)if model is not None else {},wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu)
        dump_new(ROOT/'audit/GENERATION_RECEIPT.json',receipt);raise
    dump_new(ROOT/'audit/GENERATION_RECEIPT.json',receipt)
    print(json.dumps(dict(status=receipt['status'],public_manifest=str(public/'manifest.json'),observations=3,
                         physical_LU=receipt['physical_counts']['LU'],forward_complex_RHS=receipt['physical_counts']['source_RHS'],
                         wall_seconds=receipt['wall_seconds'],truth_exposed_to_inversion=False)))

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--execute-after-freeze',action='store_true')
    args=parser.parse_args()
    if not args.execute_after_freeze:parser.error('Fresh generation requires explicit parent freeze authorization')
    generate()
