"""Authorized generation only:24 new public observations with separately sealed truth.

No inverse solver, scoring, GT plotting or parameter selection is invoked.
Original R13 geometry function bodies are extracted as source, not reimplemented.
"""
from __future__ import annotations
import os
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
import argparse
import ast
from datetime import datetime,timezone
import json
from pathlib import Path
import sys
import textwrap
import time
import traceback

FAMILIES=('multi_island_ellipses','continuous_random_field','random_level_patchwork',
          'perforated_slab_closed_holes','open_notched_porous_body','thin_bridge_network',
          'eccentric_multilayers','mixed_scale_disconnected')
VARIANTS=(('full12',.01),('full6',.02),('half12',.01))


def dump(path,value,exclusive=False):
    with Path(path).open('x' if exclusive else 'w',encoding='utf-8') as stream:
        json.dump(value,stream,indent=2,allow_nan=False,ensure_ascii=False);stream.write('\n')


def geometry_source(path):
    source=path.read_text();lines=source.splitlines(keepends=True)
    main=next(n for n in ast.parse(source).body if isinstance(n,ast.FunctionDef) and n.name=='main')
    names=('ellipse','segment','material','average','acquisition')
    functions={n.name:n for n in main.body if isinstance(n,ast.FunctionDef) and n.name in names}
    if set(functions)!=set(names):raise ValueError('Expected exact original geometry interfaces missing')
    pieces=[];locations=[]
    for name in names:
        node=functions[name];code=textwrap.dedent(''.join(lines[node.lineno-1:node.end_lineno]))
        if name=='material':
            old="high=float(rng.uniform(12,18));low=float(rng.uniform(2,5));angle=float(rng.uniform(-np.pi,np.pi));shift=rng.uniform(-.05,.05,2)"
            new="""H=float(rng.uniform(12,18));L=float(rng.uniform(2,5))
    uH=(H-12.)/6.;uL=(L-2.)/3.;regime=int(case['contrast_regime'])
    high,low=((4.+uH,2.+uL) if regime==0 else (15.+3.*uH,10.+2.*uL) if regime==1 else (H,L))
    angle=float(rng.uniform(-np.pi,np.pi));shift=rng.uniform(-.05,.05,2)"""
            if code.count(old)!=1 or code.count('h=low if v==0 else high;')!=1:
                raise ValueError('Original material draw anchors changed; refuse unspecified generator edits')
            code=code.replace(old,new).replace('h=low if v==0 else high;','h=high;')
        pieces.append(code+'\n');locations.append(dict(symbol=name,line=node.lineno,end_line=node.end_lineno))
    body='\n'.join(pieces);compile(body,str(path)+'::localized_geometry','exec')
    return body,locations


def main():
    wall0,cpu0=time.perf_counter(),time.process_time()
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--history-root',type=Path,required=True)
    parser.add_argument('--legacy-repo',type=Path,required=True)
    args=parser.parse_args();root=args.root.resolve();history=args.history_root.resolve();legacy=args.legacy_repo.resolve()
    publicdir=root/'data_public/independent24';sealdir=root/'sealed_truth/independent24';audit=root/'audit'
    audit.mkdir(parents=True,exist_ok=True)
    if any(p.exists() for p in [publicdir,sealdir,audit/'GENERATION24_RECEIPT.json',audit/'GENERATION24_PROTOCOL.json']):
        raise FileExistsError('Independent24 output already exists; no overwrite or retry')
    source=history/'experiments/generate_confirmation.py'
    body,locations=geometry_source(source)
    setup=[];costs=[];privatecases=[];case=None;failure=None
    # The source/protocol record precedes numerical imports, setup and all generated arrays.
    protocol=dict(schema_version='r14_independent24_generation_v1',status='FROZEN_BEFORE_GENERATION',
        created_at_utc=datetime.now(timezone.utc).isoformat(),objects=24,unknowns_per_inverse=576,
        inverse_grid=[24,24],forward_grid_options=[24,36],crossgrid_objects=4,
        quadrature_reference_side=216,material_parameter='real_scalar_chi_epsilon_r_minus1',
        material_box=[.02,20.],host=.1,m_init=.1,prior=dict(smooth=1e-4,tv=1e-4,epsilon=.01),
        k='2*pi',time_convention='exp(-i omega t)',ROI_width_wavelengths=.96,
        geometry_source=str(source),geometry_source_symbols=locations,
        localized_changes='Only original material first two draw remapping and h=high; all geometry branches untouched',
        private_assignment='default_rng(141499).permutation(24) assigns8families*3regimes',
        material_seed_rule='141500+opaque_index',noise_seed_rule='241500+opaque_index',
        crossgrid_rule='First two sorted opaque low-regime IDs,first sorted high ID,first sorted mixed ID; frozen before material calls',
        acquisition_rule='opaque_index%3 => full12/.01,full6/.02,half12/.01',
        regimes=dict(low='high=4+uH,low=2+uL',high='high=15+3uH,low=10+2uL',mixed='original H,L',
                     normalized_draws='uH=(H-12)/6,uL=(L-2)/3; H~U[12,18],L~U[2,5] in original RNG order'),
        noise='independent circular complex Gaussian; E|noise|^2=sigma^2; sigma=relative_RMS*RMS(clean)',
        delta='sigma*sqrt(chi2.ppf(.99,2*Y.size)/2),noise only; no model-mismatch allowance',
        public_family_seed_regime_labels=False,truth_policy='Generator writes separately; no inverse/score/GT plots',
        intended_comparison='Frozen four metric families plus parent-selected strong baselines; adaptive candidate must freeze before its new predictions',
        no_new_backend=True,no_SHA=True)
    dump(audit/'GENERATION24_PROTOCOL.json',protocol,exclusive=True)
    try:
        tick,wall=time.process_time(),time.perf_counter()
        import numpy as np
        from scipy.stats import chi2
        sys.path.insert(0,str(legacy))
        from research.research_fbg_r3_experiments.src.physics import make_grid,MaxwellModel
        setup.append(dict(category='numerical_imports',CPU_seconds=time.process_time()-tick,
                          wall_seconds=time.perf_counter()-wall,LU=0,complex_RHS=0))
        assignments=[(f,regime) for f in range(8) for regime in range(3)]
        permutation=np.random.default_rng(141499).permutation(24)
        mapping={int(permutation[j]):value for j,value in enumerate(assignments)}
        regime_ids={regime:sorted(i for i in range(24) if mapping[i][1]==regime) for regime in range(3)}
        cross=sorted(regime_ids[0][:2]+regime_ids[1][:1]+regime_ids[2][:1])
        for idx in range(24):
            family,regime=mapping[idx];acq,noise=VARIANTS[idx%3]
            privatecases.append(dict(case_id=f'r14_new_{idx:04d}',opaque_index=idx,material_seed=141500+idx,
                noise_seed=241500+idx,family=FAMILIES[family],family_index=family,
                family_variant=regime,contrast_regime=regime,contrast_regime_name=('low','high','mixed')[regime],
                forward_side=36 if idx in cross else 24,inverse_side=24,acquisition=acq,
                relative_noise_rms=noise,quadrature_reference_side=216))
        publicdir.mkdir(parents=True,exist_ok=False);sealdir.mkdir(parents=True,exist_ok=False)
        private=dict(schema_version='r14_independent24_private_v1',status='PREDECLARED_BEFORE_GENERATION',
            objects=24,families=8,objects_per_family=3,assignment_permutation_seed=141499,
            crossgrid_opaque_indices=cross,crossgrid_regime_counts=dict(low=2,high=1,mixed=1),
            expected_LU_by_dimension={'576':24,'1296':4},expected_complex_RHS_by_dimension={
                '576':240,'1296':sum(6 if privatecases[i]['acquisition']=='full6' else 12 for i in cross)},
            truth_semantics='Original216²common quadrature square-cell averages; inherited disk-pulse VIE',cases=privatecases)
        # Freeze resolved private choices before invoking any geometry or Maxwell setup.
        dump(sealdir/'GENERATOR_MANIFEST.json',private,exclusive=True)
        sourceout=root/'sources/generation24';sourceout.mkdir(parents=True,exist_ok=False)
        (sourceout/'localized_original_geometry.py').write_text(body,encoding='utf-8')
        dump(sourceout/'SOURCE_PROVENANCE.json',dict(source=str(source),symbols=locations,
            changes=protocol['localized_changes'],original_source_bytes=source.stat().st_size,no_SHA=True),exclusive=True)
        publicmanifest=dict(schema_version='r14_independent24_public_v1',status='GENERATING',objects=24,
            material_parameter='real_scalar_contrast_chi',unknowns_per_case=576,inverse_grid=[24,24],
            roi_width_wavelengths=.96,k=float(2*np.pi),lower=.02,upper=20.,initial_material=.1,
            noise_model=protocol['noise'],delta_definition=protocol['delta'],
            material_truth_in_public=False,scope='NEW_MATERIALS_FROZEN_BEFORE_PREDICTIONS',case_files=[])
        dump(publicdir/'MANIFEST.json',publicmanifest,exclusive=True)
        grids={}
        for side in [24,36]:
            tick,wall=time.process_time(),time.perf_counter()
            grids[side]=make_grid(side,sources=12,receivers=48,k=2*np.pi)
            setup.append(dict(category='original_operator_setup',side=side,Green_setups=1,LU=0,complex_RHS=0,
                CPU_seconds=time.process_time()-tick,wall_seconds=time.perf_counter()-wall))
        axis=(np.arange(216)+.5)*.96/216-.48;yy,xx=np.meshgrid(axis,axis,indexing='ij')
        namespace=dict(np=np,xx=xx,yy=yy,grids=grids)
        exec(compile(body,str(sourceout/'localized_original_geometry.py'),'exec'),namespace)
        material,average,acquisition=[namespace[k] for k in ('material','average','acquisition')]
        private_diagnostics=[]
        for case in privatecases:
            tick,wall=time.process_time(),time.perf_counter();cid=case['case_id'];models=[]
            try:
                ref,params=material(case);coarse=average(ref,24);fine=average(ref,case['forward_side'])
                forward,si,ri=acquisition(case['forward_side'],case['acquisition'])
                model=MaxwellModel(forward);models.append(model);state=model.state(fine);clean=state['Y'].copy()
                sigma=case['relative_noise_rms']*float(np.linalg.norm(clean)/np.sqrt(clean.size))
                rng=np.random.default_rng(case['noise_seed'])
                noise=sigma/np.sqrt(2)*(rng.normal(size=clean.shape)+1j*rng.normal(size=clean.shape));Y=clean+noise
                delta=float(sigma*np.sqrt(chi2.ppf(.99,2*Y.size)/2));inverse,_,_=acquisition(24,case['acquisition'])
                mismatch=None;coarse_prediction=clean
                if case['forward_side']!=24:
                    coarse_model=MaxwellModel(inverse);models.append(coarse_model)
                    coarse_prediction=coarse_model.state(coarse)['Y']
                    mismatch=float(np.linalg.norm(coarse_prediction-clean)/np.linalg.norm(clean))
                np.savez_compressed(publicdir/(cid+'.npz'),**inverse,Y=Y,sigma=sigma,delta=delta,
                    m_init=np.full(576,.1),components=1,lower=.02,upper=20.,inverse_side=24,
                    forward_side=case['forward_side'],noise_relative_rms=case['relative_noise_rms'])
                np.savez_compressed(sealdir/(cid+'.npz'),material_coarse=coarse,material_forward=fine,
                    Y_clean=clean,Y_coarse_truth=coarse_prediction,noise=noise,forward_side=case['forward_side'])
                params['case_id']=cid;dump(sealdir/(cid+'_GEOMETRY.json'),params,exclusive=True)
                private_diagnostics.append(dict(case_id=cid,crossgrid=case['forward_side']!=24,
                    clean_coarse_forward_relative_mismatch=mismatch,
                    realized_noise_over_delta=float(np.linalg.norm(noise)/delta),
                    generation_forward_relative_residual=float(state['state_relative_residual'])))
                publicmanifest['case_files'].append(dict(case_id=cid,input_file=str((publicdir/(cid+'.npz')).relative_to(root)),
                    inverse_side=24,forward_side=case['forward_side'],sources=len(si),receivers=len(ri),
                    receiver_aperture='half_circle' if case['acquisition']=='half12' else 'full_circle',
                    relative_noise_RMS=case['relative_noise_rms'],sigma=sigma,delta_noise_only=delta,
                    model_mismatch_status='SEPARATE_CROSSGRID_DIAGNOSTIC_NOT_INCLUDED_IN_NOISE_RADIUS' if case['forward_side']!=24 else 'MATCHED_DISCRETE_MODEL',
                    material_truth_in_public=False))
                dump(publicdir/'MANIFEST.json',publicmanifest)
            except BaseException:
                counts={}
                for mod in models:
                    for key,value in mod.cost.items():counts[key]=counts.get(key,0)+value
                costs.append(dict(case_id=cid,status='FAILED',counts=counts,
                    CPU_seconds=time.process_time()-tick,wall_seconds=time.perf_counter()-wall,
                    failure=traceback.format_exc(),dimension_note='Failed and partial model operations retained'))
                dump(audit/f'GENERATION24_CASE_{case["opaque_index"]:04d}.json',costs[-1],exclusive=True)
                raise
            counts={}
            for mod in models:
                for key,value in mod.cost.items():counts[key]=counts.get(key,0)+value
            costs.append(dict(case_id=cid,status='GENERATED',counts=counts,forward_dimension=case['forward_side']**2,
                inverse_dimension=576,CPU_seconds_including_serialization=time.process_time()-tick,
                wall_seconds_including_serialization=time.perf_counter()-wall,solver_run=False))
            dump(audit/f'GENERATION24_CASE_{case["opaque_index"]:04d}.json',costs[-1],exclusive=True)
            print(json.dumps(dict(case_id=cid,status='GENERATED',LU=counts.get('LU',0),complex_RHS=counts.get('actual_complex_RHS',0))),flush=True)
        private['status']='GENERATION_COMPLETE_TRUTH_SEALED';dump(sealdir/'GENERATOR_MANIFEST.json',private)
        dump(sealdir/'GENERATION_DIAGNOSTICS.json',private_diagnostics,exclusive=True)
        publicmanifest['status']='READY_FOR_FROZEN_METHOD_DISPATCH';dump(publicdir/'MANIFEST.json',publicmanifest)
        result='READY'
    except BaseException:
        failure=traceback.format_exc();result='FAILED'
    receipt=dict(status=result,objects_requested=24,objects_generated=sum(c['status']=='GENERATED' for c in costs),
        failed_cases=sum(c['status']=='FAILED' for c in costs),case_costs=costs,setup=setup,
        LU=sum(c['counts'].get('LU',0) for c in costs),complex_RHS=sum(c['counts'].get('actual_complex_RHS',0) for c in costs),
        process_CPU_seconds=time.process_time()-cpu0,process_wall_seconds=time.perf_counter()-wall0,
        public_manifest=str(publicdir/'MANIFEST.json'),private_manifest_path_not_loaded=str(sealdir/'GENERATOR_MANIFEST.json'),
        protocol=str(audit/'GENERATION24_PROTOCOL.json'),source=str(source),legacy_repo=str(legacy),
        new_inverse_runs=0,new_scoring=0,new_GT_plots=0,truth_access='GENERATOR_ONLY_WRITES_SEPARATE_PRIVATE_OUTPUT',
        failure=failure,peak_memory_bytes=None,memory_scope='Not instrumented; do not infer exact peak from array sizes',
        cost_scope='Entry through imports/source prep/setup/forward/NPZ serialization and failure work; final receipt write/print/interpreter exit excluded; no external whole-child CPU captured; no nested CPU sums',
        earlier_interpreter_probe=dict(scope='Prior runtime-audit default Python import NumPy failed before this generation task',
                                      physical_LU=0,physical_RHS=0,CPU_seconds=None,included_in_this_generation_process=False))
    dump(audit/'GENERATION24_RECEIPT.json',receipt,exclusive=True)
    print(json.dumps(dict(status=result,objects=receipt['objects_generated'],LU=receipt['LU'],complex_RHS=receipt['complex_RHS'],
                          receipt=str(audit/'GENERATION24_RECEIPT.json'))),flush=True)
    if failure:print(failure,file=sys.stderr)
    return 0 if result=='READY' else 1


if __name__=='__main__':
    raise SystemExit(main())
