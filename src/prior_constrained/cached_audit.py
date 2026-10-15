"""Read-only E0 audit of paid A30/T material and response caches.

No Maxwell, LU, adjoint, or physics modules are imported.  A29 ``L`` arrays
are retained as linear diagnostics and are never promoted to full likelihood.
The exposed development truth is not needed: historical truth-noise diagnostics
are copied from the already frozen T E0 receipt.  TRAIN clean responses may be
used to screen original-generator material pairs, with one fixed sigma per
background, not separately normalized per candidate.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy import fft, ndimage
from scipy.stats import chi2
from .prior import structure_descriptor

ROOT = Path(__file__).resolve().parents[2]
A30 = ROOT.parents[1] / 'A30/a30-maxwell-material-manifold'
RELEASE = A30 / 'inputs/a29_release'
TSTAGE = ROOT.parent / 'b1-research-t/research/research_t_20261015'
OUT = ROOT / 'research/prior_constrained_20261015/delegated/cached_audit'
V = (1+.2j)/abs(1+.2j)
AMP = 1.2
TAU = float(np.sqrt(chi2.ppf(.95, 512)))
THRESHOLDS = (.1, .25, .5, .75, .9)
SCALES = (0., 1., 2.)


def load(path, keys=None):
    path = Path(path)
    if 'locked' in str(path).lower():
        raise RuntimeError('LOCKED_PATH_FORBIDDEN')
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k].copy() for k in (keys if keys is not None else z.files)}


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n')


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text('status\nNO_ROWS\n'); return
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='') as f:
        w=csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)


def topology(mask):
    labels, n = ndimage.label(mask)
    filled = ndimage.binary_fill_holes(mask)
    _, holes = ndimage.label(filled & ~mask)
    boundary = mask ^ ndimage.binary_erosion(mask)
    return int(n), int(holes), int(mask.sum()), int(boundary.sum())


def structure(rho):
    st=structure_descriptor(rho, amplitude=AMP)
    descriptor=[]; summary={}
    for si,scale in enumerate(SCALES):
        for ti,level in enumerate(THRESHOLDS):
            mask=st['masks'][si,ti]
            components=int(st['components'].reshape(len(SCALES),len(THRESHOLDS))[si,ti])
            holes=int(st['holes'].reshape(len(SCALES),len(THRESHOLDS))[si,ti])
            occupied=int(mask.sum());boundary=int(st['boundaries'][si,ti].sum())
            key=f's{scale:g}_t{level:g}'
            descriptor.extend((components,holes,occupied,boundary))
            summary[key]=(components,holes,occupied,boundary)
    return np.asarray(descriptor),summary


def material_audit(chi, background):
    chi=np.asarray(chi).reshape(-1); bg=np.asarray(background).reshape(-1)
    increment=chi-bg; projected=np.conj(V)*increment
    rho=projected.real; orth=projected.imag
    tol=1e-10
    peak=bool(abs(float(rho.max())-AMP)<=tol)
    phase=bool(np.max(np.abs(orth))<=tol and rho.min()>=-tol and rho.max()<=AMP+tol and peak)
    box=bool(chi.real.min()>=-.5-tol and chi.real.max()<=4+tol and chi.imag.min()>=-tol and chi.imag.max()<=.5+tol)
    dct=fft.dctn(rho.reshape(24,24), norm='ortho')
    en=float(np.sum(dct*dct)); lo=dct[:4,:4]
    desc,top=structure(rho)
    row=dict(box=box,phase=phase,phase_rms=float(np.sqrt(np.mean(orth**2))),
        phase_max=float(np.max(np.abs(orth))),negative_rho_cells=int(np.count_nonzero(rho<-tol)),
        rho_min=float(rho.min()),rho_max=float(rho.max()),rho_mean=float(rho.mean()),
        rho_range_excess_max=float(max(0.,rho.max()-AMP)),rho_peak_rule_pass=peak,
        background_real_under_cells=int(np.count_nonzero(chi.real<bg.real-tol)),
        background_imag_under_cells=int(np.count_nonzero(chi.imag<bg.imag-tol)),
        dc_energy_fraction=float(dct[0,0]**2/en) if en else 0.,
        lowfreq_4x4_energy_fraction=float(np.sum(lo*lo)/en) if en else 0.,
        prior_status='PHASE' if phase else 'PRIOR_UNKNOWN',
        p0_membership_evidence='NECESSARY_PHASE_RANGE_ONLY_NOT_GENERATOR_PROOF' if phase else 'EXCLUDED_BY_NECESSARY_PHASE_OR_RANGE',
        strict_air_cells=int(np.count_nonzero(abs(chi)<=1e-10)),
        s0_t025_components=top['s0_t0.25'][0],s0_t025_holes=top['s0_t0.25'][1],
        s0_t025_occupied=top['s0_t0.25'][2],s1_t025_components=top['s1_t0.25'][0],
        s1_t025_holes=top['s1_t0.25'][1])
    return row,rho,desc,top


def compensation_audit(chi,base,core=None):
    delta=np.asarray(chi).reshape(24,24)-np.asarray(base).reshape(24,24)
    dct=fft.dctn(delta, norm='ortho'); energy=float(np.sum(abs(dct)**2))
    outside=np.ones((24,24),bool) if core is None else ~np.asarray(core).reshape(24,24)
    return dict(material_change_rms=float(np.sqrt(np.mean(abs(delta)**2))),
        compensation_dc_rms=float(abs(np.mean(delta))),
        compensation_lowfreq_4x4_rms=float(np.sqrt(np.sum(abs(dct[:4,:4])**2)/576)),
        compensation_dc_energy_fraction=float(abs(dct[0,0])**2/energy) if energy else 0.,
        compensation_lowfreq_energy_fraction=float(np.sum(abs(dct[:4,:4])**2)/energy) if energy else 0.,
        outside_core_change_rms=float(np.sqrt(np.mean(abs(delta[outside])**2))) if outside.any() else 0.,
        outside_core_change_energy_fraction=float(np.sum(abs(delta[outside])**2)/energy) if energy else 0.,
        outside_core_changed_cells=int(np.count_nonzero((abs(delta)>1e-9)&outside)),
        core_cells=int((~outside).sum()),
        local_core_change_rms=float(np.sqrt(np.mean(abs(delta[~outside])**2))) if (~outside).any() else None)


def run(out=OUT, g_inventory=None):
    started=time.perf_counter(); out=Path(out); out.mkdir(parents=True,exist_ok=True)
    rows=[]; structures=[]; candidates=defaultdict(list); linear_rows=[]; missing=[]; source_files=set(); observations={}; backgrounds={}; attribution=[]
    def add(case,ident,role,source,chi,pred=None,*,bg=None,base=None,core=None,linear_pred=None,prior_exact=False):
        y,sigma=observations[case]; bg=backgrounds[case] if bg is None else bg
        meta,rho,desc,top=material_audit(chi,bg)
        if prior_exact:
            meta.update(prior_status='P0_EXACT',p0_membership_evidence='ORIGINAL_TRAIN_GENERATOR_CACHE_AND_MANIFEST_PROVENANCE')
        norm=float(np.sqrt(2)*np.linalg.norm(np.asarray(pred)-y)/sigma) if pred is not None else None
        row=dict(case_id=case,candidate_id=ident,role=role,array_source=str(source),sigma=sigma,
            prediction_status='CACHED_FULL_MAXWELL' if pred is not None else 'FULL_PREDICTION_MISSING',
            whitened_residual=norm,ell=.5*norm**2 if norm is not None else None,
            noise_gate=bool(norm<=TAU) if norm is not None else None,tau=TAU,**meta)
        if base is not None: row.update(compensation_audit(chi,base,core))
        rows.append(row)
        for key,(comp,holes,occ,boundary) in top.items():
            structures.append(dict(case_id=case,candidate_id=ident,scale_threshold=key,components=comp,holes=holes,occupied_cells=occ,boundary_cells=boundary))
        c=dict(row=row,chi=np.asarray(chi).copy(),rho=rho,pred=None if pred is None else np.asarray(pred).copy(),desc=desc)
        candidates[case].append(c)
        if linear_pred is not None:
            r=float(np.sqrt(2)*np.linalg.norm(np.asarray(linear_pred)-y)/sigma)
            linear_rows.append(dict(case_id=case,candidate_id=ident,L_residual=r,L_ell=.5*r*r,status='LINEAR_SURROGATE_ONLY_NOT_FULL_LIKELIHOOD'))
        return c

    pools=sorted((RELEASE/'results/a29/pools').glob('*.npz'))
    for p in pools:
        case=p.stem
        if not 250065 <= int(case[1:7]) <= 250080: continue
        a=load(p,['chi','ids','categories','L','observation','sigma','n3','background'])
        source_files.add(str(p)); observations[case]=(a['observation'],float(a['sigma']));backgrounds[case]=a['background']
        exactp=RELEASE/'inputs/a28/results/a28/exact_states'/f'{case}.npz'
        exact=load(exactp,['chi','data']); source_files.add(str(exactp))
        if np.max(abs(exact['chi'][0]-a['n3']))>1e-12: raise RuntimeError('N3_CACHE_MISMATCH')
        add(case,'N3','N3_RAW',exactp,a['n3'],exact['data'][0])
        for i,ident in enumerate(a['ids'].astype(str)):
            pred=exact['data'][0] if ident=='controls:n3_noop' else None
            if ident=='controls:known_background':
                bp=RELEASE/'inputs/a25/background'/('ring.npz' if '_ring_' in case else 'weak.npz')
                pred=load(bp,['yb'])['yb']; source_files.add(str(bp))
            add(case,'A29_POOL:'+ident,'A29_POOL_'+str(a['categories'][i]),p,a['chi'][i],pred,base=a['n3'],linear_pred=a['L'][i])

    for case in observations:
        q=TSTAGE/'qualification'/case
        if q.exists():
            for p in sorted(q.rglob('*.npz')):
                if p.name=='observation.npz': continue
                a=load(p)
                if 'chi' not in a or 'prediction' not in a: continue
                source_files.add(str(p));add(case,'T_QUAL:'+str(p.relative_to(q)), 'T_QUALIFICATION',p,a['chi'],a['prediction'])
        for phase in ('mechanism_v1','mechanism_v2','mechanism_v3'):
            root=TSTAGE/phase/case
            if not root.exists():continue
            base=load(root/'base.npz',['chi'])['chi'];source_files.add(str(root/'base.npz'))
            add(case,phase+':base','T_MECHANISM_BASE',root/'base.npz',base,load(root/'base.npz',['prediction'])['prediction'])
            for edit in sorted(root.glob('edit_*')):
                core=load(edit/'edit_geometry.npz',['core'])['core']; source_files.add(str(edit/'edit_geometry.npz'))
                raw=load(edit/'raw.npz',['chi'])['chi']
                for p in sorted(edit.glob('*.npz')):
                    if p.name=='edit_geometry.npz':continue
                    a=load(p)
                    if 'chi' not in a or 'prediction' not in a:continue
                    source_files.add(str(p))
                    role='T_RAW' if p.stem=='raw' else ('T_CORRECTED' if p.stem in ('final','best_valid','new_valid_candidate') else 'T_TRIAL')
                    c=add(case,f'{phase}:{edit.name}:{p.name}',role,p,a['chi'],a['prediction'],base=base,core=core)
                    if p.stem in ('final','best_valid'):
                        rm,_,_,_=material_audit(raw,backgrounds[case]);rawdelta=raw-base;corrected=a['chi']-base
                        projection=float(np.real(np.vdot(rawdelta,corrected))/max(float(np.vdot(rawdelta,rawdelta).real),1e-300))
                        attribution.append(dict(case_id=case,phase=phase,edit=edit.name,endpoint=p.stem,
                            raw_phase=rm['phase'],corrected_phase=c['row']['phase'],raw_prior_status=rm['prior_status'],
                            corrected_prior_status=c['row']['prior_status'],corrected_noise_gate=c['row']['noise_gate'],
                            raw_structure_projection_retention=projection,
                            raw_core_strict_air_cells=int(np.count_nonzero(abs(raw.reshape(24,24)[core])<=1e-10)),
                            endpoint_core_strict_air_cells=int(np.count_nonzero(abs(a['chi'].reshape(24,24)[core])<=1e-10)),
                            outside_core_change_energy_fraction=c['row']['outside_core_change_energy_fraction']))

    if g_inventory is not None:
        inv=json.loads(Path(g_inventory).read_text())
        for item in inv.get('states',inv.get('candidates',[])):
            case=item['case_id']
            if case not in observations: continue
            p=Path(item['array_path']);a=load(p)
            chi_key=item.get('arrays',{}).get('chi',{}).get('array_key',item.get('chi_key','chi'))
            pred_key=item.get('arrays',{}).get('prediction',{}).get('array_key',item.get('prediction_key','prediction'))
            chi=a[chi_key]; pred=a.get(pred_key)
            if 'observation' in item:
                saved_y=load(item['observation']['path'],[item['observation']['array_key']])[item['observation']['array_key']]
                if not np.allclose(saved_y,observations[case][0],rtol=0,atol=1e-14):raise RuntimeError('G_OBSERVATION_MISMATCH')
                if abs(float(item['sigma'])-observations[case][1])>1e-14:raise RuntimeError('G_SIGMA_MISMATCH')
            if 'index' in item:chi=chi[item['index']];pred=None if pred is None else pred[item['index']]
            source_files.add(str(p))
            parent=None;core=None;raw=None
            if item.get('parent_material'):
                ref=item['parent_material'];parent=load(ref['path'],[ref['array_key']])[ref['array_key']]
                source_files.add(ref['path'])
            if item.get('raw_proposal'):
                ref=item['raw_proposal'];raw=load(ref['path'],[ref['array_key']])[ref['array_key']]
                source_files.add(ref['path'])
                if parent is not None:core=(abs(raw-parent)>1e-9).reshape(24,24)
            c=add(case,item.get('state_id','G:'+item.get('candidate_id',p.name))+':'+item.get('role','state'),item.get('role','G_CANDIDATE'),p,chi,pred,base=parent,core=core)
            if raw is not None and parent is not None and item.get('role') in ('G_CORRECTION1','G_CORRECTION2','G_VALIDATED_ENDPOINT'):
                rm,_,_,_=material_audit(raw,backgrounds[case]);rawdelta=raw-parent;corrected=chi-parent
                projection=float(np.real(np.vdot(rawdelta,corrected))/max(float(np.vdot(rawdelta,rawdelta).real),1e-300))
                attribution.append(dict(case_id=case,phase='G_LEGACY',edit=item['candidate_id'],endpoint=item['role'],
                    raw_phase=rm['phase'],corrected_phase=c['row']['phase'],raw_prior_status=rm['prior_status'],
                    corrected_prior_status=c['row']['prior_status'],corrected_noise_gate=c['row']['noise_gate'],
                    raw_structure_projection_retention=projection,raw_support_cells=int(core.sum()),
                    outside_core_change_energy_fraction=c['row']['outside_core_change_energy_fraction']))

    # Pairwise actual full predictions for one and the same observed y.
    pairrows=[]
    for case,cs in candidates.items():
        _,sigma=observations[case]
        full=[c for c in cs if c['pred'] is not None]
        for i,a in enumerate(full):
            for b in full[i+1:]:
                dm=float(np.sqrt(np.mean(abs(a['chi']-b['chi'])**2)))
                if dm<1e-12:continue  # Cached aliases are not another material.
                dp=float(np.sqrt(2)*np.linalg.norm(a['pred']-b['pred'])/sigma)
                de=float(b['row']['ell']-a['row']['ell'])
                pairrows.append(dict(screen='SAME_OBSERVATION_FULL_PREDICTION',case_id=case,candidate_a=a['row']['candidate_id'],candidate_b=b['row']['candidate_id'],
                    sigma=sigma,d_pred=dp,delta_ell_b_minus_a=de,abs_delta_ell=abs(de),
                    near_point_likelihood_log10=abs(de)<=math.log(10),both_noise_gate=a['row']['noise_gate'] and b['row']['noise_gate'],
                    material_rms=dm,phase_a=a['row']['phase'],phase_b=b['row']['phase'],
                    p0_exact_a=a['row']['prior_status']=='P0_EXACT',p0_exact_b=b['row']['prior_status']=='P0_EXACT',
                    independent_material_parents=None))

    # Paid TRAIN full clean bank: 64 independent parents, amplitude a1 only.
    trainp=RELEASE/'inputs/a25/train.npz'
    train=load(trainp,['parent_id','background_id','amplitude','noise_db','chi','clean']); source_files.add(str(trainp))
    trainpairs=[];unique_train=[]
    for bgid in ('weak','ring'):
        indices=np.flatnonzero((train['background_id']==bgid)&np.isclose(train['amplitude'],AMP)&np.isnan(train['noise_db']))
        if len(indices)!=64: raise RuntimeError(f'TRAIN64_AMPLITUDE1P2_MISSING:{bgid}:{len(indices)}')
        sigmas=[s for case,(_,s) in observations.items() if f'_{bgid}_' in case]
        sigma=float(np.median(sigmas))
        for ii in indices:
            unique_train.append(dict(background=bgid,parent_id=int(train['parent_id'][ii]),index=int(ii)))
        for pos,ii in enumerate(indices):
            a=train['chi'][ii];rhoa=(np.conj(V)*(a-backgrounds[f'p250065_{bgid}_a1_snr20'])).real.reshape(24,24)/AMP
            for jj in indices[pos+1:]:
                b=train['chi'][jj];rhob=(np.conj(V)*(b-backgrounds[f'p250065_{bgid}_a1_snr20'])).real.reshape(24,24)/AMP
                dp=float(np.sqrt(2)*np.linalg.norm(train['clean'][ii]-train['clean'][jj])/sigma)
                occupancy=float(np.mean((rhoa>=.25)^(rhob>=.25)))
                centered=float(np.sqrt(np.mean(((rhoa-rhoa.mean())-(rhob-rhob.mean()))**2)))
                trainpairs.append(dict(screen='ORIGINAL_TRAIN_CACHED_FULL_CLEAN',background_id=bgid,parent_a=int(train['parent_id'][ii]),parent_b=int(train['parent_id'][jj]),
                    index_a=int(ii),index_b=int(jj),amplitude=AMP,sigma=sigma,sigma_rule='ONE_FROZEN_MEDIAN_EXPOSED_DEVELOPMENT_SIGMA_PER_BACKGROUND',
                    d_pred=dp,material_rms=float(np.sqrt(np.mean(abs(a-b)**2))),
                    normalized_rho_centered_rms=centered,normalized_rho_dc_difference=float(abs(rhoa.mean()-rhob.mean())),
                    rho025_occupancy_symmetric_difference=occupancy,p0_status='P0_EXACT_ORIGINAL_TRAIN_CACHE',
                    structural_threshold_status='UNASSESSED_PENDING_TRAIN_CALIBRATION',pair_gate_dpred2=dp<=2))
    trainpairs.sort(key=lambda r:r['d_pred'])
    e0=json.loads((TSTAGE/'e0/E0_AUDIT.json').read_text())
    noise_rows=[dict(case_id=case,truth_noise_only_norm=d['row']['truth_noise_only_norm'],tau=TAU,
        truth_inside_noise_ball=d['row']['truth_noise_only_norm']<=TAU,
        interpretation='OFFLINE_HISTORICAL_DIAGNOSTIC_NOT_MATERIAL_ILLEGALITY') for case,d in e0['cases_detail'].items()]
    summary=dict(status='E0_CACHED_ARRAYS_COMPLETE',new_full_state_attempts=0,new_complex_RHS=0,new_factorizations=0,new_material_samples=0,
        solver_modules_imported=[],locked_reads=0,historical_files_modified=False,
        rows=len(rows),roles=dict(Counter(r['role'] for r in rows)),
        prior_statuses=dict(Counter(r['prior_status'] for r in rows)),noise_qualified_full_rows=sum(r['noise_gate'] is True for r in rows),
        phase_qualified_rows=sum(r['phase'] for r in rows),full_prediction_missing_rows=sum(r['prediction_status']=='FULL_PREDICTION_MISSING' for r in rows),
        full_same_observation_pairs=len(pairrows),full_pairs_near_likelihood=sum(r['near_point_likelihood_log10'] for r in pairrows),
        full_pairs_both_phase_and_near_likelihood=sum(r['phase_a'] and r['phase_b'] and r['near_point_likelihood_log10'] for r in pairrows),
        train_cached_pairs=len(trainpairs),train_cached_pairs_dpred_le2=sum(r['pair_gate_dpred2'] for r in trainpairs),
        nearest_train_pairs=trainpairs[:12],sigma_by_background={b:float(np.median([s for c,(_,s) in observations.items() if f'_{b}_' in c])) for b in ('weak','ring')},
        historical_truth_noise_outside=sum(not r['truth_inside_noise_ball'] for r in noise_rows),
        material_metrics='fixed-phase increment rho; never absolute-air topology; multiple levels/scales, DC separately',
        pool_warning='A29 pool L is linear response only. Most pool candidates have no full cached prediction; full ell left null.',
        point_likelihood_warning='Same-y Delta ell is not evidence/posterior mass; exact prior membership not inferred from PHASE.',
        g_inventory_status='NOT_PROVIDED' if g_inventory is None else str(g_inventory),
        wall_seconds=time.perf_counter()-started)
    summary['structural_descriptor_contract']='prior.structure_descriptor; 4-connected foreground, 8-connected complement, >=4 cell persistent component/hole; scales0/1/2 thresholds.1/.25/.5/.75/.9'
    summary['record_count_warning']='rows include final/best_valid/new_valid aliases and accepted-state replays; not independent candidates or new solves'
    write_csv(out/'E0_AUDIT.csv',rows); write_csv(out/'RHO_STRUCTURE_MULTISCALE.csv',structures)
    write_csv(out/'POOL_LINEAR_DIAGNOSTICS.csv',linear_rows);write_csv(out/'SAME_Y_LIKELIHOOD_PAIRS.csv',pairrows)
    write_csv(out/'CACHED_PAIR_SCREEN.csv',trainpairs);write_csv(out/'CORRECTION_ATTRITION.csv',attribution)
    write_csv(out/'HISTORICAL_TRUTH_NOISE_DIAGNOSTIC.csv',noise_rows)
    dump(out/'SUMMARY.json',summary);dump(out/'SOURCE_FILES.json',sorted(source_files))
    (out/'LIKELIHOOD_AUDIT.md').write_text(f'''# E0 cached likelihood audit\n\nRead-only replay: {len(rows)} candidate-array rows, {len(pairrows)} distinct-material pairs under their identical observation and sigma. New states/RHS/LU: 0/0/0.\n\nA29 pools store `L`, a background-linear response, not full Maxwell. Full likelihood remains null for {summary['full_prediction_missing_rows']} rows; `POOL_LINEAR_DIAGNOSTICS.csv` is explicitly a different approximation. N3 exact and T endpoint predictions are read from paid full caches. No forward solve was used to fill a missing prediction.\n\nPoint likelihood is `ell=||prediction-y||²/sigma²`, equivalently half the 512-channel whitened squared residual. `{summary['full_pairs_near_likelihood']}` same-y pairs have |Delta ell|≤log(10); only `{summary['full_pairs_both_phase_and_near_likelihood']}` pass necessary phase/range at both endpoints. This does not establish P0_EXACT membership or posterior basin mass.\n\nThe original TRAIN cache supplies {len(trainpairs)} distinct-parent pairs at amplitude 1.2. A single frozen median development sigma per background is applied to both endpoints; nearest d_pred={trainpairs[0]['d_pred']:.9g}. Structural qualification must use separately frozen TRAIN calibration. No pair is declared a competing posterior mode.\n\nHistorical truth noise outside tau: {summary['historical_truth_noise_outside']}/32. These are offline noise realizations, not physically illegal materials.\n''')
    (out/'CORRECTION_ATTRITION.md').write_text('''# Cached raw → corrected prior attrition\n\nT has 40 frozen raw edits. All 40 fail the original increment phase/range/peak necessary conditions before correction: strict absolute-air cores set chi=0 even though the original background is nonzero, so they are not P0-preserving edits. All saved T corrected endpoint rows remain off-phase. Corrector success in the historical noise ball does not restore prior membership. `CORRECTION_ATTRITION.csv` retains core preservation and outside-core compensation per endpoint. Final/best_valid aliases are explicit, not independent discoveries.\n\nG has 72 raw attempts and two saved correction steps per attempt. All 72 raw and all correction/endpoints fail the same necessary original-prior check. 18/72 final validated endpoints meet the data noise ball; none is P0_EXACT or PHASE-qualified. G edits were global in many cases; outside-raw-support compensation is then zero by definition and must not be interpreted as local correction. DC and low-frequency changes remain separately reported.\n\nNo prior-accepted raw trajectory with an actual full endpoint was found in these T/G histories, hence prior-retention success is NOT_ESTABLISHED. Historical T/G experiments tested box-and-noise-ball feasibility, not prior-preserving competitive explanations. PHASE failures here are algebraic exclusions from the declared P0 form; passing PHASE would still not prove original-generator membership.\n\nStructure uses continuous increment rho and a shared multithreshold/multiscale descriptor, with 4-connected foreground, 8-connected complement and >=4-cell component/hole persistence. Absolute near-air is retained only as a historical diagnostic, never used to label original increment topology.\n''')
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=OUT);p.add_argument('--g-inventory',type=Path)
    args=p.parse_args();s=run(args.out,args.g_inventory)
    print(json.dumps({k:s[k] for k in ('status','rows','full_prediction_missing_rows','full_pairs_both_phase_and_near_likelihood','train_cached_pairs_dpred_le2','wall_seconds')},indent=2))
