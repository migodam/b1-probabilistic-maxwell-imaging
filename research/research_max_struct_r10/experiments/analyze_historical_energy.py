"""Actually execute sparse saved-state R9 energy diagnostics on one CPU thread."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import argparse,csv,json,sys,time
from pathlib import Path
from collections import defaultdict
import numpy as np
from threadpoolctl import threadpool_limits
HERE=Path(__file__).resolve().parents[1];REPO=HERE.parents[1]
sys.path.insert(0,str(REPO))
from research.research_max_struct_r10.src.energy_diagnostics import Geometry,physical_current,energy,ref
OUT=HERE/'delegated/energy_audit'
R9=REPO.parent/'r9_worktree/research/research_max_struct_r9'

def stamp(p):
    s=p.stat();return dict(path=str(p),bytes=s.st_size,mtime_ns=s.st_mtime_ns)
def load(p):
    with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def csvwrite(p,rows):
    if not rows:return
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=lambda z:z.item() if isinstance(z,np.generic) else str(z))+'\n')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cohort',choices=['development','confirmation','vector','all'],default='development');args=ap.parse_args()
    out=OUT/args.cohort;out.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter();cpu=time.process_time();singles=[];twos=[];rem=[];states=[];aff=[];sources=[];apertures=[];failures=[];geoms=[];physical_cache={};LU=0;RHS=0
    cohorts=['development','confirmation','vector'] if args.cohort=='all' else [args.cohort]
    runs=[]
    for cohort in cohorts:
        if cohort=='vector':
            for p in sorted((R9/'RESULTS/vector_extension').glob('*/*/prediction.npz')):
                receipt=json.loads((p.parent/'RECEIPT.json').read_text())
                if receipt.get('arm') in ('global_dual','a_plus','global_dual_qualified','a_plus_qualified','dual_source'):
                    runs.append((cohort,p.parent,receipt))
            continue
        for p in sorted((R9/'RESULTS'/cohort).glob('*/prediction.npz')):
            if not (p.parent/'RECEIPT.json').exists():continue
            receipt=json.loads((p.parent/'RECEIPT.json').read_text());arm=receipt.get('arm')
            if arm not in ('global_dual','a_plus','cie_dual','a_plus_sourcewise'):continue
            if cohort=='development' and not any(s in p.parent.name for s in ('b20d','40ac')):continue
            runs.append((cohort,p.parent,receipt))
    for cohort,run,receipt in runs:
        try:
            case=receipt['case_id'];arm=receipt['arm'];public={};vector=cohort=='vector'
            components=3 if vector else 1
            if vector:
                from research.research_max_struct_r9.src.vector_extension import StrictVectorInverse
                pp=R9/'data_public/vector_extension'/case/'PUBLIC.json';pm=json.loads(pp.read_text())
                cfg=json.loads((run/'FROZEN_CONFIG.json').read_text())
                model=StrictVectorInverse(8,pm['acquisition']['domain_width'],cfg,auxiliary=False)
                py=pp.with_name('observations.npz')
                public=dict(G=model.G,S=model.S,E0=model.inc,Y=load(py)['Y'],shape=np.array([8,8,8]),
                    grid_shape=np.array([8,8,8]),receivers=model.receivers,sources=model.receivers)
                gtpath=R9/'sealed_truth/vector_extension'/case/'truth.npz';truth=load(gtpath)['chi']
                sources.extend([stamp(pp),stamp(py),stamp(run/'FROZEN_CONFIG.json'),stamp(Path(sys.modules[StrictVectorInverse.__module__].__file__))])
                task=dict(dpcs_rule=receipt['dpcs']['rule'])
            else:
                task=receipt['task']
                for pp in task['public_files']:
                    pp=Path(pp);public.update(load(pp));sources.append(stamp(pp))
                gtpath=R9/'sealed_truth'/cohort/(case+'.npz');truth=load(gtpath)['m']
            sources.append(stamp(gtpath))
            G,S,E,Y=public['G'],public['S'],public['E0'],public['Y']
            apertures.append(dict(task_id=run.name,case_id=case,source_positions=public['sources'].tolist(),receiver_positions=public['receivers'].tolist()))
            match=next((i for i,g in enumerate(geoms) if np.array_equal(g.G,G) and np.array_equal(g.S,S)),None)
            if match is None:
                print(json.dumps(dict(event='spectra_start',case=case,shape=G.shape)),flush=True)
                geoms.append(Geometry(G,S));match=len(geoms)-1
            geom=geoms[match]
            key=(match,truth.tobytes(),E.tobytes())
            if key not in physical_cache:physical_cache[key]=physical_current(G,E,np.repeat(truth,components));LU+=1;RHS+=E.shape[1]
            jgt=physical_cache[key]
            stagepaths=sorted(run.glob('stage_*_joint.npz'))
            stagemeta={p.stem:json.loads(p.with_suffix('.json').read_text()) for p in stagepaths}
            rank=receipt['dpcs']['rank'] if vector else next(iter(stagemeta.values()))['receiver_rank']
            paths=[run/'initial.npz']+stagepaths+[run/('raw.npz' if vector else 'selected_before_corrector.npz'),run/'prediction.npz']
            for p in paths:
                if not p.exists():continue
                z=load(p);sources.append(stamp(p));m=z['m'];k=(match,m.tobytes(),E.tobytes())
                if k not in physical_cache:
                    if 'j' in z or (vector and '_joint' not in p.stem):physical_cache[k]=z.get('j',z.get('J'))
                    else:
                        counterpart=p.with_name(p.name.replace('_joint','_physical'))
                        if counterpart.exists() and 'j' in (zz:=load(counterpart)) and np.array_equal(zz['m'],m):
                            physical_cache[k]=zz['j'];sources.append(stamp(counterpart))
                        else:physical_cache[k]=physical_current(G,E,np.repeat(m,components));LU+=1;RHS+=E.shape[1]
                j=physical_cache[k];dm=m-truth;dj=j-jgt
                base=dict(cohort=cohort,case_id=case,task_id=run.name,arm=arm,rule=task.get('dpcs_rule'),state=p.stem,
                    operator_group=match,truth_access='EXPOSED_OFFLINE_DIAGNOSTIC',sources=E.shape[1],receiver_rows=S.shape[0],
                    receiver_rank=rank,mf=(stagemeta.get(p.stem,{}).get('mf') if vector else stagemeta.get(p.stem,{}).get('basis_info',{}).get('mf')),source_snapshot=str(p),
                    acquisition_metadata='APERTURES.json keyed by task_id and case_id',
                    current_constraint='actual affine FFT-TSOM' if '_joint' in p.stem else 'physical current; not affine-constrained',
                    material_coordinates='isometric repeated-vector lift dm/sqrt(3); coordinate control only' if vector else 'raw scalar grid material')
                row=dict(**base,material_error_energy=energy(dm),true_current_error_energy=energy(dj),
                    internal_error_energy=energy(G@dj),receiver_true_error_energy=energy(S@dj),
                    physical_data_residual=float(np.linalg.norm(S@j-Y)/np.linalg.norm(Y)),
                    linear_physical_defect_energy=energy(j-np.repeat(m,components)[:,None]*(E+G@j)))
                dmlift=np.repeat(dm,components)/np.sqrt(components)
                quantities=[('material' if not vector else 'material_isometric_lift_control',dmlift,'complex'),('material_real' if not vector else 'material_isometric_lift_real_control',dmlift,'real'),('true_current',dj,'complex')]
                if 'J' in z and '_joint' in p.stem:
                    J=z['J'];daux=J-jgt;consistency=J-np.repeat(m,components)[:,None]*(E+G@J);gap=J-j
                    row.update(auxiliary_current_error_energy=energy(daux),physical_consistency_defect_energy=energy(consistency),
                        auxiliary_minus_physical_current_energy=energy(gap),
                        consistency_identity_closure_energy=energy(gap-np.repeat(m,components)[:,None]*(G@gap)-consistency),
                        auxiliary_vs_physical_signed_cross_energy=energy(daux)-energy(dj)-energy(gap))
                    quantities += [('auxiliary_current_error',daux,'complex'),('physical_consistency_defect',consistency,'complex'),('auxiliary_minus_physical_current',gap,'complex')]
                    mf=base['mf'];JD,U,q=geom.affine(Y,rank,tuple(public.get('grid_shape',public.get('shape'))),mf,components)
                    off=jgt-JD;proj=q@(q.conj().T@off);irreducible=off-proj
                    representation=J-JD-q@(q.conj().T@(J-JD))
                    affine_error=J-JD-proj
                    aff.append(dict(**base,fourier_columns=U.shape[1],coefficient_range_rank=q.shape[1],
                        affine_GT_irreducible_energy=energy(irreducible),affine_GT_projectable_energy=energy(proj),
                        affine_auxiliary_in_range_error_energy=energy(affine_error),representation_violation_energy=energy(representation),
                        saved_coefficient_reconstruction_error_energy=energy(J-JD-U@z['gamma' if vector else 'A']),
                        decomposition_closure=energy(daux)-energy(irreducible)-energy(affine_error),
                        JD_energy=energy(JD),data_hardlock_error_energy=energy(geom.receiver_vh[:rank]@(J-JD)),
                        aperture_receivers=repr(public.get('receivers').tolist() if 'receivers' in public else 'NOT_SAVED')))
                    quantities += [('affine_GT_irreducible_current',irreducible,'complex')]
                states.append(row)
                for quantity,x,family in quantities:
                    a,b,c=geom.decompose(x,family)
                    for rows,target in [(a,singles),(b,twos),(c,rem)]:
                        target.extend(dict(**base,quantity=quantity,**r) for r in rows)
            print(json.dumps(dict(event='run_complete',task_id=run.name,states=len(paths),elapsed_seconds=time.perf_counter()-start)),flush=True)
        except Exception as exc:
            import traceback
            failures.append(dict(run=str(run),error=repr(exc),traceback=traceback.format_exc()));print(failures[-1],flush=True)
    for name,rows in [('state_metrics',states),('single_operator_energy',singles),('twofold_energy',twos),('twofold_remainders_cross_energy',rem),('actual_FFT_TSom_affine',aff)]:csvwrite(out/(name+'.csv'),rows)
    # Matched algorithm deltas: absolute ordinary-minus-A+ band energies and proportions.
    index={(r['case_id'],r['rule'],r['state'],r['quantity'],r['tau'],r['operator'],r['band'],r['arm']):r for r in singles}
    paired=[]
    for r in singles:
        if r['arm'] not in ('a_plus','a_plus_sourcewise','a_plus_qualified'):continue
        ordinary={'a_plus':'global_dual','a_plus_qualified':'global_dual_qualified','a_plus_sourcewise':'cie_dual'}[r['arm']]
        key=(r['case_id'],r['rule'],r['state'],r['quantity'],r['tau'],r['operator'],r['band'],ordinary)
        if key not in index:continue
        b=index[key];paired.append(dict(**{k:r[k] for k in ('cohort','case_id','rule','state','quantity','tau','operator','band')},
            ordinary=ordinary,aplus=r['arm'],ordinary_total_energy=b['total_energy'],aplus_total_energy=r['total_energy'],
            ordinary_band_energy=b['component_energy'],aplus_band_energy=r['component_energy'],
            absolute_band_energy_reduction=b['component_energy']-r['component_energy'],
            ordinary_fraction=b['fraction'],aplus_fraction=r['fraction']))
    csvwrite(out/'paired_absolute_energy_reduction.csv',paired)
    dump(out/'OPERATOR_METADATA.json',[dict(group=i,G_shape=g.G.shape,S_shape=g.S.shape,spaces=g.meta) for i,g in enumerate(geoms)])
    dump(out/'SOURCE_PROVENANCE.json',list({r['path']:r for r in sources}.values()))
    dump(out/'APERTURES.json',apertures)
    dump(out/'FAILED_ATTEMPTS.json',failures)
    receipt=dict(status='COMPLETED' if not failures else 'COMPLETED_WITH_FAILURES',cohort=args.cohort,runs=len(runs),saved_states=len(states),operator_groups=len(geoms),
        truth_access='EXPOSED_OFFLINE_DIAGNOSTIC',new_inverse_runs=0,new_physical_LU=LU,new_complex_RHS=RHS,
        wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,GPU=False,threads=1,
        spectra_reused_from_reference=False,spectra_source='fresh SVD of exact historical operators; grouping only by np.array_equal G/S',
        failure_count=len(failures),cost_scope='includes public/GT reads, fresh spectra/joint geometry, exact missing currents, affine SVD, all diagnostic decomposition and serialization',
        reference_linear_algebra_cost=ref.cost)
    dump(out/'EXECUTION_RECEIPT.json',receipt);print(json.dumps(receipt),flush=True)
if __name__=='__main__':
    with threadpool_limits(limits=1):main()
