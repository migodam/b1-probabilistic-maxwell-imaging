"""Generator-only disclosed R2 geometry replay and optional32-grid generation.

--verify-geometry-only replays existing RNG and compares existing24 alpha only.
--generate writes one new disclosed-development observation; never a blind case.
No inverse module imports this file. No installation, backend or remote control.
"""
from pathlib import Path
from copy import deepcopy
from datetime import datetime,timezone
import argparse,json,sys,time,traceback
import numpy as np
from scipy.stats import chi2
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parent
STUDY=ROOT.parent
REPO=STUDY.parents[1]
sys.path.insert(0,str(REPO))
from research.research_max_struct_r2.experiments.generate_blind import shape_alpha as original_shape_alpha
from research.research_z_alg_r4_diversity.src.physics import Model,build_operators


def shape_support192(family,rng):
    # Exact original192 geometry logic. Only return changes: binary support,
    # allowing both8x8(old24) and6x6(new32) cell area averages of the SAME raster.
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
    return support


def recover_geometry(cfg):
    """No Maxwell call/new data: validate family shuffle and old24 geometry exactly."""
    r2=REPO/'research/research_max_struct_r2'
    mapping=json.loads((r2/'sealed_truth/mapping_PRIVATE.json').read_text())
    if int(mapping['seed_PRIVATE'])!=cfg['original_generator_seed']:
        raise RuntimeError('Original disclosed seed differs; refuse guessed replay')
    rng=np.random.default_rng(cfg['original_generator_seed'])
    families=['ring','curved_bar','two_islands'];rng.shuffle(families)
    if families!=[row['family'] for row in mapping['records']]:
        raise RuntimeError('Original family shuffle/order not reproduced')
    target=None;checks=[]
    for family,row in zip(families,mapping['records']):
        state=deepcopy(rng.bit_generator.state)
        support=shape_support192(family,rng)
        alpha24=support.reshape(24,8,24,8).mean(axis=(1,3))
        ref_rng=np.random.default_rng();ref_rng.bit_generator.state=state
        original_alpha=original_shape_alpha(family,ref_rng)
        if not np.array_equal(alpha24,original_alpha) or rng.bit_generator.state!=ref_rng.bit_generator.state:
            raise RuntimeError('Copied geometry or RNG consumption differs from original generator')
        if row['observation_id']==cfg['original_observation_id']:
            if family!='curved_bar':raise RuntimeError('Original ID/family mismatch')
            with np.load(r2/'sealed_truth'/row['truth_file'],allow_pickle=False) as truth:
                old_alpha=truth['alpha']
                if not np.array_equal(alpha24,old_alpha):
                    raise RuntimeError('Same-object24 alpha identity FAIL; stop without new32 data')
            target=support.copy()
            checks.append(dict(observation_id=row['observation_id'],family=family,
                original24_alpha_equal=True,max_abs_alpha_error=float(np.max(abs(alpha24-old_alpha))),
                copied_generator_and_rng_state_equal=True))
        # Original forward is deterministic and consumes no RNG. Original noise
        # consumes two complete normal arrays, even on objects before the target.
        rng.normal(size=(48,12));rng.normal(size=(48,12))
    if target is None:raise RuntimeError('Original target absent; stop')
    return target,dict(status='EXACT_EXISTING24_GEOMETRY_REPLAY',classification=cfg['classification'],
        family_order=families,normal_arrays_consumed_per_object=2,normal_array_shape=[48,12],checks=checks,
        new_data_generated=False,physical_forward_calls=0)


def dump(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x',encoding='utf-8') as stream:json.dump(obj,stream,indent=2,allow_nan=False)


def generate(cfg):
    wall,cpu=time.perf_counter(),time.process_time()
    support,replay=recover_geometry(cfg)
    public,truth,audit=ROOT/'data_public',ROOT/'development_truth',ROOT/'audit/generation'
    if any(path.exists() for path in (public,truth,audit)):
        raise FileExistsError('Preserve existing diagnostic/attempt; no regeneration')
    for path in (public,truth,audit):path.mkdir(parents=True)
    receipt=dict(status='RUNNING',classification=cfg['classification'],truth_access='GENERATOR_ONLY_DISCLOSED_R2_DEVELOPMENT',
        geometry_replay=replay,failed_attempts=[],operator_setup={},physical_counts={})
    model=None
    try:
        dump(audit/'FROZEN_CONFIG.json',cfg)
        # Read old PUBLIC geometry/acquisition only; no old Y used by the solve.
        with np.load(REPO/cfg['original_public_path'],allow_pickle=False) as data:
            sources,receivers,k=data['sources'],data['receivers'],float(data['k'])
            old_centers=data['coords'];old_radius=float(data['radius'])
        if (not np.allclose([old_centers.min(),old_centers.max(),old_radius],[-.46,.46,.018],rtol=0,atol=1e-14)
            or sources.shape!=(12,2) or receivers.shape!=(48,2)):
            raise ValueError('Original geometry/acquisition not the declared disk grid')
        side=32;h=.96/side;axis=(np.arange(side)+.5)*h-.48
        yy,xx=np.meshgrid(axis,axis,indexing='ij');centers=np.c_[xx.ravel(),yy.ravel()]
        tick,ct=time.perf_counter(),time.process_time()
        with threadpool_limits(limits=1):
            ops=build_operators(centers,sources,receivers,.45*h,k)
            receipt['operator_setup']=dict(wall_seconds=time.perf_counter()-tick,cpu_seconds=time.process_time()-ct,
                builds=1,n=1024,spacing=h,radius=.45*h,box=[-.48,.48])
            alpha=support.reshape(32,6,32,6).mean(axis=(1,3));m=.1+12*alpha.ravel()
            model=Model(ops);state=model.forward(m,label='disclosed_development_grid32_generation')
        clean=state['y'];sigma=.001*np.linalg.norm(clean)/np.sqrt(clean.size)
        rng=np.random.default_rng(cfg['new_noise_seed'])
        noise=sigma/np.sqrt(2)*(rng.normal(size=clean.shape)+1j*rng.normal(size=clean.shape))
        noise_radius=float(np.sqrt(chi2.ppf(.99,2*clean.size)/2))
        np.savez_compressed(truth/'grid32_curved_bar_truth.npz',alpha=alpha,m=m,mask=alpha>=.5,
                            support192=support,Y_clean=clean,noise=noise)
        payload={key:ops[key] for key in ('G','E0','S')}
        payload.update(Y=clean+noise,grid_shape=np.array([32,32]),shape=np.array([32,32]),coords=centers,
            sources=sources,receivers=receivers,k=np.array(k),radius=np.array(.45*h),host=np.array(.1),
            amplitude=np.array(12.),sigma=np.array(sigma),noise_radius=np.array(noise_radius),
            material_bounds=np.array([.02,20.]),smooth_weight=np.array(1e-4),tv_weight=np.array(1e-4),
            tv_epsilon=np.array(.01),prior_width=np.array(19.5))
        np.savez_compressed(public/'grid32_curved_bar.npz',**payload)
        dump(public/'manifest.json',dict(status='DISCLOSED_DEVELOPMENT_OBSERVATION',classification=cfg['classification'],
            path='grid32_curved_bar.npz',shape=[32,32],case_count=1,source_count=12,receiver_count=48,
            payload_keys=sorted(payload),truth_access='NONE_TO_INVERSION'))
        receipt.update(status='GENERATED_DISCLOSED_DEVELOPMENT_GRID32',physical_counts=dict(model.counts),
            linear_relative_residual=float(state['state_relative_residual']))
    except Exception as exc:
        receipt.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc())
        receipt['failed_attempts'].append(dict(error=repr(exc),charged=True));raise
    finally:
        receipt.update(wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,
            physical_counts=dict(model.counts) if model else {},timer_scope='replay/existing geometry read/operator assembly/one full forward/new noise/all writes; final receipt write excluded')
        dump(audit/'GENERATION_RECEIPT.json',receipt)
    print(json.dumps(dict(status=receipt['status'],classification=cfg['classification'],physical_counts=receipt['physical_counts'])))


if __name__=='__main__':
    parser=argparse.ArgumentParser();mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--verify-geometry-only',action='store_true');mode.add_argument('--generate',action='store_true')
    args=parser.parse_args();cfg=json.loads((ROOT/'config.json').read_text())
    if args.verify_geometry_only:
        _,record=recover_geometry(cfg)
        path=ROOT/'audit/GEOMETRY_REPLAY_PREPARATION.json'
        dump(path,record);print(json.dumps(record))
    else:generate(cfg)
