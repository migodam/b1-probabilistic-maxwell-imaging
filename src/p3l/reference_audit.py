"""Read-only diagnostics for four explicit saved P3-L chains.

No target/backend/book or truth is imported. Samples retain rejection repeats.
Passing finite-chain gates never qualifies a reference without a separate
cross-mechanism comparison. Outputs describe only the supplied chain files.
"""
from pathlib import Path
import json

import numpy as np

from p1_material.prior import saturation
from .diagnostics import material_diagnostics


PROJECTION_SEED = 2026102250
QUANTILE_LEVELS = np.array([.025, .5, .975])
RHAT_LIMIT = 1.01
MATERIAL_ESS_MIN = 1000.
MEAN_RMS_MCSE_MAX = .001
N_PIXELS = 576
DECODE_BATCH = 1024
PIXEL_DIAGNOSTIC_BATCH = 32


def _jsonable(value):
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, np.generic):
        return _jsonable(value.item())
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def decode_material_batches(q, prior):
    """Vectorized original MaterialPrior law, without calls per draw.

    Uses the supplied original prior's field map/config, original saturation,
    bounded property intervals, host and ROI. Accepts (...,66) and returns
    (...,576); no Gaussian approximation or posterior independence is assumed.
    """
    q = np.asarray(q, float)
    if q.ndim < 2 or q.shape[-1] != 66 or not np.all(np.isfinite(q)):
        raise ValueError('finite q with shape (...,66) required')
    if np.any((q[...,64:] < 0) | (q[...,64:] > 1)):
        raise ValueError('chain contains an illegal uniform property')
    field_map = np.asarray(prior.A, float)
    mask = np.asarray(prior.mask, bool)
    if (field_map.shape != (576,64) or mask.shape != (576,)
            or np.count_nonzero(mask) != 256):
        raise ValueError('original P1 576-pixel/256-ROI field contract required')
    if not np.all(np.isfinite(field_map)) or not np.isfinite(prior.temperature) or prior.temperature <= 0:
        raise ValueError('invalid supplied original prior field map/halfwidth')
    real, imag = np.asarray(prior.real_interval,float), np.asarray(prior.imag_interval,float)
    if real.shape != (2,) or imag.shape != (2,) or not np.all(np.isfinite(np.r_[real,imag])):
        raise ValueError('invalid supplied property intervals')
    flat = q.reshape(-1,66)
    material = np.empty((len(flat),576), complex)
    for start in range(0,len(flat),DECODE_BATCH):
        values = flat[start:start+DECODE_BATCH]
        field = float(prior.config['field_mean'])+values[:,:64]@field_map.T
        fraction = saturation(field/float(prior.temperature))*mask
        increment = real[0]+(real[1]-real[0])*values[:,64]+1j*(imag[0]+(imag[1]-imag[0])*values[:,65])
        material[start:start+len(values)] = complex(prior.host)+increment[:,None]*fraction
    if not np.all(np.isfinite(material)):
        raise ValueError('decoder produced nonfinite material')
    return material.reshape(q.shape[:-1]+(576,))


def _finite_max(values):
    a = np.asarray(values,float)
    return float(np.max(a)) if a.size and np.all(np.isfinite(a)) else None


def _finite_min(values):
    a = np.asarray(values,float)
    return float(np.min(a)) if a.size and np.all(np.isfinite(a)) else None


def _summary(diag):
    status = np.asarray(diag['status'])
    names, counts = np.unique(status,return_counts=True)
    return dict(status_counts=dict(zip(names.tolist(),counts.tolist())),
                maximum_rhat=_finite_max(diag['rhat']),
                minimum_bulk_ess=_finite_min(diag['bulk_ess']),
                minimum_tail_ess=_finite_min(diag['tail_ess']),
                minimum_mean_ess=_finite_min(diag['mean_ess']),
                all_diagnostics_assessed=bool(np.all(status=='DIAGNOSTICS_ONLY')))


def _material_functions(material, roi):
    indices = np.flatnonzero(roi)
    # Evenly spaced positions in sorted fixed ROI indices, independent of data.
    selected = indices[np.rint(np.linspace(0,len(indices)-1,16)).astype(int)]
    rng = np.random.default_rng(PROJECTION_SEED)
    weights = rng.normal(size=(8,len(indices)))
    weights /= np.linalg.norm(weights,axis=1,keepdims=True)
    values = material[...,indices]
    blocks = [values.real.mean(axis=-1),values.imag.mean(axis=-1),
              (values.real**2).mean(axis=-1),(values.imag**2).mean(axis=-1)]
    names = ['ROI_REAL_MEAN','ROI_IMAG_MEAN','ROI_REAL_MEAN_SQUARE','ROI_IMAG_MEAN_SQUARE']
    for pixel in selected:
        blocks.extend([material[...,pixel].real,material[...,pixel].imag])
        names.extend([f'PIXEL_{pixel}_REAL',f'PIXEL_{pixel}_IMAG'])
    projected = values@weights.T
    for k in range(8):
        blocks.extend([projected[...,k].real,projected[...,k].imag])
        names.extend([f'ROI_PROJECTION_{k}_REAL',f'ROI_PROJECTION_{k}_IMAG'])
    return np.stack(blocks,axis=-1),names,selected,weights


def audit_reference(paths, prior, warmup, output_dir):
    """Audit exactly four saved CHAIN.npz files, then write JSON/NPZ.

    Equal chain lengths are required; no thinning, truncating to shortest,
    deduplication, truth access or adaptive function choice occurs. The files'
    saved warmup, when present, must agree with the explicit argument. Metadata
    target-family absence prevents gate qualification; mismatched families
    raise. File identity and independence of seeds remain caller obligations.
    """
    paths = [Path(p).expanduser().resolve() for p in paths]
    if len(paths)!=4 or len(set(paths))!=4:
        raise ValueError('exactly four distinct saved chain paths required')
    if not isinstance(warmup,(int,np.integer)) or warmup<0:
        raise ValueError('warmup must be a nonnegative integer')
    qs, phis, metadata = [],[],[]
    for path in paths:
        if path.name!='CHAIN.npz' or any(p.lower() in ('offline_only','locked_truth','truth') for p in path.parts):
            raise ValueError('only explicit saved CHAIN.npz inference outputs may be read')
        with np.load(path,allow_pickle=False) as data:
            if not {'q','Phi'} <= set(data.files):
                raise ValueError('chain requires q and Phi')
            q = np.asarray(data['q'],float)
            phi = np.asarray(data['Phi'],float)
            saved_warmup = int(data['warmup']) if 'warmup' in data else None
            family = str(data['target_family'].item()) if 'target_family' in data else None
        if q.ndim!=2 or q.shape[1]!=66 or phi.shape!=(len(q),):
            raise ValueError('saved q/Phi chain layout mismatch')
        if len(q)-warmup<4 or not np.all(np.isfinite(q)) or not np.all(np.isfinite(phi)):
            raise ValueError('need finite chains and at least four retained production draws')
        if saved_warmup is not None and saved_warmup!=warmup:
            raise ValueError('explicit warmup disagrees with saved metadata')
        qs.append(q[warmup:].copy());phis.append(phi[warmup:].copy())
        metadata.append(dict(path=str(path),saved_total_draws=len(q),saved_warmup=saved_warmup,target_family=family))
    if len({len(q) for q in qs})!=1:
        raise ValueError('chain lengths differ; no implicit truncation permitted')
    families = {r['target_family'] for r in metadata if r['target_family'] is not None}
    if len(families)>1:
        raise ValueError('cannot pool different target families')
    q,phi = np.stack(qs),np.stack(phis)
    material = decode_material_batches(q,prior)
    roi = np.asarray(prior.mask,bool)
    functions,names,selected,weights = _material_functions(material,roi)
    material_diag = material_diagnostics(functions)
    latent_diag = material_diagnostics(q[:,:,:64])
    property_diag = material_diagnostics(q[:,:,64:])
    phi_diag = material_diagnostics(phi[:,:,None])
    # Full pixel MCSE is computed in modest blocks, preserving constants.
    pixel_diag = {}
    for start in range(0,576,PIXEL_DIAGNOSTIC_BATCH):
        block = material_diagnostics(material[:,:,start:start+PIXEL_DIAGNOSTIC_BATCH])
        for key,value in block.items():
            if isinstance(value,np.ndarray):
                pixel_diag.setdefault(key,[]).append(value)
    pixel_diag = {key:np.concatenate(value,axis=0) for key,value in pixel_diag.items()}
    mean = material[0,0]+np.mean(material-material[0,0],axis=(0,1))
    variance = np.mean(np.abs(material-mean)**2,axis=(0,1))
    quantile_real = np.quantile(material.real,QUANTILE_LEVELS,axis=(0,1))
    quantile_imag = np.quantile(material.imag,QUANTILE_LEVELS,axis=(0,1))
    pixel_mcse = pixel_diag['mean_mcse']
    # Full-grid RMS: sum both components, divide by all576; outside host is
    # deterministic, but only ROI assessment can establish finite-chain gates.
    rms_mcse = float(np.sqrt(np.sum(pixel_mcse[roi]**2)/576)) if np.all(np.isfinite(pixel_mcse[roi])) else None
    roi_status_ok = bool(np.all(pixel_diag['status'][roi]=='DIAGNOSTICS_ONLY'))
    diagnostics = dict(latent=latent_diag,properties=property_diag,Phi=phi_diag,material_functions=material_diag)
    rhat_ok = all(np.all(np.isfinite(v['rhat'])) and np.all(v['rhat']<RHAT_LIMIT) and np.all(v['status']=='DIAGNOSTICS_ONLY') for v in diagnostics.values())
    roi_rhat_ok = bool(np.all(np.isfinite(pixel_diag['rhat'][roi])) and np.all(pixel_diag['rhat'][roi]<RHAT_LIMIT))
    ess_ok = bool(np.all(np.isfinite(material_diag['bulk_ess'])) and np.all(material_diag['bulk_ess']>=MATERIAL_ESS_MIN) and np.all(np.isfinite(material_diag['tail_ess'])) and np.all(material_diag['tail_ess']>=MATERIAL_ESS_MIN))
    mcse_ok = rms_mcse is not None and rms_mcse<=MEAN_RMS_MCSE_MAX
    metadata_ok = all(m['target_family']=='FULL' for m in metadata)
    gates = dict(target_metadata_full=metadata_ok,all_named_rhat_below_limit=bool(rhat_ok),
                 roi_pixel_rhat_below_limit=roi_rhat_ok,all_roi_components_assessed=roi_status_ok,
                 all_material_function_bulk_and_tail_ess_at_least_minimum=ess_ok,
                 full_grid_mean_rms_mcse_at_most_limit=bool(mcse_ok))
    within_gate = all(gates.values())
    outside = pixel_diag['status'][~roi]
    stats = dict(status='WITHIN_MECHANISM_GATES_PASS_CROSSCHECK_REQUIRED' if within_gate else 'UNQUALIFIED_WITHIN_MECHANISM_DIAGNOSTICS',
                 reference_qualified=False,cross_mechanism_status='NOT_EVALUATED_SEPARATE_INDEPENDENT_MECHANISM_REQUIRED',
                 warmup_removed=int(warmup),chain_count=4,production_draws_per_chain=q.shape[1],
                 chain_sources=metadata,target_family=next(iter(families)) if families else 'UNKNOWN',
                 thresholds=dict(rhat_max_exclusive=RHAT_LIMIT,minimum_material_function_bulk_and_tail_ess=MATERIAL_ESS_MIN,full_grid_mean_rms_mcse_max_inclusive=MEAN_RMS_MCSE_MAX),
                 full_grid_mean_rms_mcse=rms_mcse,gates=gates,
                 summaries={key:_summary(value) for key,value in diagnostics.items()},
                 named_diagnostics=diagnostics,material_function_names=names,
                 pixel_diagnostics=pixel_diag,
                 fixed_host_pixel_count=int(np.sum(~roi)),
                 fixed_host_component_status='KNOWN_DETERMINISTIC_HOST_EXCLUDED_FROM_MIXING_GATES_ONLY',
                 outside_raw_constant_component_count=int(np.sum(outside=='CONSTANT_UNASSESSED')),
                 roi_constant_component_count=int(np.sum(pixel_diag['status'][roi]=='CONSTANT_UNASSESSED')),
                 caution='Finite-chain diagnostics, not global mixing/coverage or GT-risk evidence; all constants explicitly retained; quantiles are marginal components.')
    output = Path(output_dir)
    output.mkdir(parents=True,exist_ok=True)
    arrays = dict(material_mean=mean,material_variance=variance,quantile_levels=QUANTILE_LEVELS,
                  quantile_real=quantile_real,quantile_imag=quantile_imag,roi_mask=roi,
                  selected_roi_indices=selected,projection_weights=weights,projection_seed=np.array(PROJECTION_SEED),
                  material_functions=functions,material_function_names=np.array(names),
                  production_q=q,production_Phi=phi)
    for group,diag in {**diagnostics,'pixel':pixel_diag}.items():
        for key,value in diag.items():
            if isinstance(value,np.ndarray):arrays[f'{group}_{key}']=value
    np.savez_compressed(output/'REFERENCE_ARRAYS.npz',**arrays)
    clean = _jsonable(stats)
    (output/'REFERENCE_STATS.json').write_text(json.dumps(clean,indent=2,allow_nan=False)+'\n')
    return clean
