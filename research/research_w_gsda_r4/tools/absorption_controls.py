#!/usr/bin/env python3
"""E1 offline physical fixtures: full FEM, exact absorption, original GS, reopen.

This is explicitly an offline physics audit. It reads authorized private fixture
parameters and writes their actual fields; it is not an online reconstruction.
"""
from pathlib import Path
import argparse
import json
import os
import sys
import time
sys.dont_write_bytecode=True
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[name]='1'
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from gsda.accounting import BudgetBook,BudgetExceeded,dump
from gsda.fem_model import FEMModel
from gsda.currents import FixedReference,CurrentGreen
from gsda.assimilation import absorb_wall,reopen,relative

CONTROLS=dict(c010='ellipse target, nonuniform passive wall',
              c011='nonconvex target, one exactly air wall cell',
              c012='holed strong target, nonuniform passive wall')


def save_state(path,z,greens,currents,fields,data,**extra):
    arrays=dict(z=z,data=data,frequencies=np.asarray([g.p.freq for g in greens]),
                axis_order='frequency,source,receiver',**extra)
    for ii,(g,j,field) in enumerate(zip(greens,currents,fields)):
        arrays[f'j_f{ii}']=j;arrays[f'total_field_f{ii}']=field
        arrays[f'quadrature_points_f{ii}']=g.p.material_points[g.active]
        arrays[f'quadrature_weights_f{ii}']=g.weights
        arrays[f'absorbed_contrast_f{ii}']=g.absorbed
    np.savez_compressed(path,**arrays)


def finite_state(currents,fields,data):
    return bool(np.isfinite(data).all() and all(np.isfinite(a).all() for a in currents+fields))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',type=Path,default=ROOT/'configs/base.json')
    ap.add_argument('--output',type=Path,default=ROOT/'results/absorption_controls')
    ap.add_argument('--cases',nargs='+',choices=tuple(CONTROLS),default=tuple(CONTROLS))
    ap.add_argument('--job',default='absorption_controls')
    args=ap.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
    if (out/'CONTROL_SUMMARY.json').exists():raise RuntimeError('refusing overwrite of completed control audit')
    if len(set(args.cases))!=len(args.cases):raise ValueError('duplicate controls')
    config=json.loads(args.config.read_text())
    if len(config['frequencies'])!=3:raise ValueError('E1 controls require all three physical frequencies')
    fixtures={}
    for case in args.cases:
        path=ROOT/'inputs/private'/case/'truth.npz'
        with np.load(path,allow_pickle=False) as source:
            z=source['z'].copy();clean=source['clean'].copy() if 'clean' in source.files else None
        if not np.isfinite(z).all():raise ValueError('nonfinite fixture '+case)
        if case=='c011' and not (z[1]==1 and z[7]==0):raise ValueError('c011 must retain its declared air wall cell')
        if case=='c012' and np.max(z[12:61])<6:raise ValueError('c012 must retain the declared strong target')
        fixtures[case]=(z,clean)
    dump(out/'OFFLINE_ACCESS_CONTRACT.json',dict(online_reconstruction=False,
        authorized_private_fixture_access=True,private_paths=[str(ROOT/'inputs/private'/c/'truth.npz') for c in args.cases],
        fixtures={c:CONTROLS[c] for c in args.cases},full_frequencies=config['frequencies'],
        selection='first half of public wall-cell indices, including the exactly air c011 cell, for exact mechanism audit',
        interpretation='same-discrete FEM mechanism control, not inversion recovery or continuum qualification'))
    book=BudgetBook(ROOT,job=args.job);status='FAILED';records=[];fatal=None
    try:
        book.check(new_job=True);model=FEMModel(config,book,phase='R4_E1_OFFLINE_ABSORPTION')
        reference=FixedReference(model)
        for case in args.cases:
            caseout=out/case;caseout.mkdir(parents=True,exist_ok=True)
            before=book.snapshot(pid=os.getpid());cpu=time.process_time();wall=time.perf_counter()
            record=dict(case=case,fixture=CONTROLS[case],truth_access=True,online_reconstruction=False)
            try:
                z,clean=fixtures[case];z=model.validate(z)
                data,js,fields,diag=reference.physical_state(z)
                assembled=model.predict(z)
                native=np.asarray([g.ev.direct+g.apply_GS(j) for g,j in zip(reference.greens,js)])
                save_state(caseout/'full_state.npz',z,reference.greens,js,fields,data,
                           assembled_FEM_data=assembled.data,original_GS0_reconstructed_data=native,
                           current_scope='full physical contrast density chi*E')
                selected=np.arange(int(np.ceil(.5*model.nwall)))
                greens,jr,ea,frozen,absorb=absorb_wall(reference,z,js,fields,fraction=.5,tolerance=1e-8,selected=selected)
                absorb['selection']='fixed first half of public wall-cell indices; c011 includes its exactly air cell'
                absorbed_data=np.asarray([g.ev.direct+g.data_bg+g.apply_GS(j) for g,j in zip(greens,jr)])
                full_reconstructed=[j+g.absorbed[:,None]*field for g,j,field in zip(greens,jr,ea)]
                original_data=np.asarray([g.ev.direct+g.apply_GS(j) for g,j in zip(reference.greens,full_reconstructed)])
                extra=dict(original_GS0_reconstructed_data=original_data,background_scattered_data=np.asarray([g.data_bg for g in greens]),
                           current_scope='residual current after exact wall absorption; full current also retained')
                for ii,j in enumerate(full_reconstructed):extra[f'full_physical_current_f{ii}']=j
                save_state(caseout/'absorbed_state.npz',z,greens,jr,ea,absorbed_data,**extra)
                reopened,ef,reopen_event=reopen(reference,greens,jr,ea,z,tolerance=1e-8)
                reopened_data=np.asarray([g.ev.direct+g.apply_GS(j) for g,j in zip(reference.greens,reopened)])
                save_state(caseout/'reopened_state.npz',z,reference.greens,reopened,ef,reopened_data,
                           current_scope='full current after genuine active-variable REOPEN')
                # A separate all-wall background isolates target feedback. No
                # remaining wall contrast is included in the Born residual.
                allwall=[];born_currents=[];wall_fields=[];wall_data=[];born_induced=[]
                full_target_currents=[]
                for ii,(g0,field) in enumerate(zip(reference.greens,fields)):
                    chiwall=np.where(g0.wall,g0.chi(z),0.)
                    g0.factor=None
                    ga=CurrentGreen(model,ii,absorbed=chiwall,reference='all_wall_born_control')
                    ga.absorbed_mask=g0.wall.copy()
                    target_chi=ga.chi(z)
                    jb=target_chi[:,None]*ga.inc
                    allwall.append(ga);born_currents.append(jb);wall_fields.append(ga.inc)
                    wall_data.append(ga.ev.direct+ga.data_bg);born_induced.append(ga.apply_GS(jb))
                    full_target_currents.append(target_chi[:,None]*field)
                wall_data=np.asarray(wall_data);born_induced=np.asarray(born_induced)
                full_induced=data-wall_data;born_data=wall_data+born_induced
                error=born_data-data
                born_metrics=dict(whitened_L2_error=float(np.linalg.norm(error)/model.sigma),
                    whitened_RMS_independent_real_imag=float(np.linalg.norm(error)/(model.sigma*np.sqrt(2*error.size))),
                    relative_target_induced_error=float(np.linalg.norm(born_induced-full_induced)/max(np.linalg.norm(full_induced),1e-14)),
                    relative_error_denominator='max(norm full target induced, 1e-14)',
                    full_target_induced_norm=float(np.linalg.norm(full_induced)),
                    Born_target_induced_norm=float(np.linalg.norm(born_induced)),
                    approximation='jBorn_target=chi_target*wall_only_total_field',
                    remaining_wall_residual=False,inverse_uses_Born=False,
                    interpretation='measured physical failure boundary, not a pass/fail target-recovery gate')
                extra=dict(full_FEM_data=data,wall_only_data=wall_data,full_target_induced=full_induced,
                    Born_target_induced=born_induced,current_scope='target-only Born current in exact all-wall background',
                    field_scope='all-wall-only total illumination used by Born approximation')
                for ii,(j,field) in enumerate(zip(full_target_currents,fields)):
                    extra[f'full_target_current_f{ii}']=j;extra[f'full_FEM_total_field_f{ii}']=field
                save_state(caseout/'allwall_born_comparison.npz',z,allwall,born_currents,wall_fields,born_data,**extra)
                dump(caseout/'BORN_COMPARISON.json',born_metrics)
                for ga in allwall:ga.factor=None
                checks=dict(full_FEM_assembled_vs_native_current=relative(assembled.data,data),
                    full_FEM_vs_original_GS0=relative(data,native),full_FEM_vs_absorbed=relative(data,absorbed_data),
                    full_FEM_vs_original_GS0_reconstruction=relative(data,original_data),
                    full_FEM_vs_REOPEN=relative(data,reopened_data))
                if clean is not None:checks['offline_generation_clean_vs_FEM']=relative(clean,assembled.data)
                fields_finite=finite_state(js,fields,data) and finite_state(jr,ea,absorbed_data) and finite_state(reopened,ef,reopened_data)
                record.update(status='PASS' if fields_finite and max(checks.values())<1e-8 else 'FAIL',
                    all_actual_arrays_finite=fields_finite,relative_checks=checks,tolerance=1e-8,
                    full_FEM_diagnostics=diag,assembled_FEM_diagnostics=assembled.diagnostics,
                    exact_absorption=absorb,reopen=reopen_event,allwall_Born_comparison=born_metrics,
                    files=['full_state.npz','absorbed_state.npz','reopened_state.npz','allwall_born_comparison.npz','BORN_COMPARISON.json'])
                if record['status']!='PASS':book.failure(model.phase,'E1 fixture mismatch',case,checks=checks)
            except BaseException as error:
                record.update(status='FAILED',error=dict(type=type(error).__name__,message=str(error)))
                book.failure('R4_E1_OFFLINE_ABSORPTION',case,error)
                if isinstance(error,(BudgetExceeded,KeyboardInterrupt,SystemExit)):raise
            finally:
                after=book.snapshot(pid=os.getpid())
                record.update(case_rhs=after['rhs']-before['rhs'],case_full_fe=after['full_fe']-before['full_fe'],
                    case_process_cpu_seconds=time.process_time()-cpu,case_wall_seconds=time.perf_counter()-wall,
                    shared_reference_setup_excluded_from_case_cost=True)
                dump(caseout/'CONTROL_RESULT.json',record);records.append(record)
            dump(out/'CONTROL_PROGRESS.json',dict(records=records))
        status='PASS' if all(r['status']=='PASS' for r in records) else 'FAIL'
    except BaseException as error:
        fatal=dict(type=type(error).__name__,message=str(error));book.failure('R4_E1_OFFLINE_ABSORPTION','controls',error)
        if isinstance(error,(KeyboardInterrupt,SystemExit)):raise
    finally:
        receipt=book.receipt(status=status,output=str(out),offline_fixture_access=True)
        dump(out/'RECEIPT.json',receipt)
        dump(out/'CONTROL_SUMMARY.json',dict(status=status,records=records,error=fatal,
            total_cost_includes_reference_setup=receipt,online_reconstruction=False,
            scope='three declared nondigit same-discrete physical controls; no shape-recovery claim'))
    print(json.dumps(dict(status=status,cases=[dict(case=r['case'],status=r['status']) for r in records],
                         output=str(out),cost=book.snapshot(pid=os.getpid()))),flush=True)
    if status!='PASS':raise SystemExit(1)

if __name__=='__main__':main()
