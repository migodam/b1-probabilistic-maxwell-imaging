#!/usr/bin/env python3
"""Saved-state offline scoring/gallery. Never imports or calls a solver.

The exposed A30 truth archive is read only in this separate process.  The
resulting metrics, figures and material arrays must never initialize searches.
"""
from __future__ import annotations
import argparse
import csv
import html
import json
import os
from pathlib import Path
from time import perf_counter, process_time

START, CPU = perf_counter(), process_time()
import numpy as np
from scipy import ndimage
import scipy.linalg
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
STAGE = ROOT / "research/research_t_20261015"
A30 = ROOT.parents[1] / "A30/a30-maxwell-material-manifold"
AIR_THRESHOLDS = (.075, .1, .125)
F4 = ndimage.generate_binary_structure(2, 1)
B8 = ndimage.generate_binary_structure(2, 2)
REUSE_FIGURES = False


def forbidden(*args, **kwargs):
    raise RuntimeError("OFFLINE_SCORING_FORBIDS_PHYSICS_SOLVES")


for obj, names in ((np.linalg, ("solve", "inv", "pinv", "lstsq", "cholesky")),
                   (scipy.linalg, ("solve", "lu_solve", "lu_factor", "cho_solve", "inv", "lstsq"))):
    for name in names:
        setattr(obj, name, forbidden)


def load(path):
    if any("locked" in part.lower() for part in Path(path).parts):
        raise ValueError("LOCKED_PATH_FORBIDDEN")
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k].copy() for k in z.files}


def plain(value):
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, dict): return {k: plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [plain(v) for v in value]
    return value


def dump(path, value):
    path.write_text(json.dumps(plain(value), indent=2, allow_nan=False)+"\n")


def csvwrite(path, rows):
    if not rows: return
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v) if isinstance(v, (dict,list)) else v for k,v in row.items()})


def iou(a, b):
    union = int(np.count_nonzero(a | b))
    return float(np.count_nonzero(a & b)/union) if union else 1.


def boundaries(mask):
    mask = np.asarray(mask, bool).reshape(24,24)
    return mask & ~ndimage.binary_erosion(mask, F4, border_value=0)


def boundary_metrics(a, b):
    ba, bb = boundaries(a), boundaries(b)
    if not ba.any() or not bb.any():
        return dict(hausdorff_m=None, symmetric_mean_m=None, status="UNDEFINED_EMPTY_BOUNDARY")
    ab = ndimage.distance_transform_edt(~bb)[ba]/24
    b_a = ndimage.distance_transform_edt(~ba)[bb]/24
    return dict(hausdorff_m=float(max(ab.max(),b_a.max())),
                symmetric_mean_m=float((ab.sum()+b_a.sum())/(len(ab)+len(b_a))),status="RASTER_CENTER_DISTANCE")


def topology(material):
    material=np.asarray(material,bool).reshape(24,24)
    _, components=ndimage.label(material,F4)
    labels, aircomponents=ndimage.label(~material,B8)
    exterior=set(map(int,np.r_[labels[0],labels[-1],labels[:,0],labels[:,-1]]))-{0}
    return dict(material_components=int(components),enclosed_air_cavities=int(aircomponents)-len(exterior),
                air_components=int(aircomponents),air_cells=int(np.count_nonzero(~material)))


def validate_saved(state, observation, sigma, tau):
    chi=state['chi'].reshape(-1);prediction=state['prediction']
    residual=float(np.linalg.norm((prediction-observation)*np.sqrt(2)/sigma))
    legal=bool(np.isfinite(chi).all() and np.all(chi.real>=-.5-1e-12) and np.all(chi.real<=4+1e-12)
               and np.all(chi.imag>=-1e-12) and np.all(chi.imag<=.5+1e-12))
    con=None
    if 'field' in state and 'current' in state:
        con=float(np.linalg.norm(state['current']-chi[:,None]*state['field'])/max(float(np.linalg.norm(state['current'])),1e-300))
    back=float(state['backward_residual']) if 'backward_residual' in state else None
    valid=bool(legal and residual<=tau and con is not None and con<=1e-10 and back is not None and back<=1e-10)
    return dict(residual_norm=residual,tau=float(tau),data_compatible=bool(residual<=tau),material_legal=legal,
                constitutive_defect=con,backward_residual=back,
                validation_status="FULL_VALIDATED_NUMERICAL" if valid else "REJECTED",
                output_error_upper_bound=None)


def metrics(chi, truth, background, base=None):
    chi,truth,background=[np.asarray(v).reshape(24,24) for v in (chi,truth,background)]
    error=chi-truth
    threshold=.25*float(np.max(np.abs(truth-background)))
    historical_a=np.abs(chi-background)>threshold
    historical_t=np.abs(truth-background)>threshold
    online_a=np.abs(chi-background)>.18
    online_t=np.abs(truth-background)>.18
    roi=historical_t
    data=dict(full_complex_MAE=float(np.mean(np.abs(error))),full_Re_MAE=float(np.mean(np.abs(error.real))),
        full_Im_MAE=float(np.mean(np.abs(error.imag))),
        complex_permittivity_NMAE=float(np.sum(np.abs(error))/np.sum(np.abs(1+truth))),
        susceptibility_NMAE=float(np.sum(np.abs(error))/max(float(np.sum(np.abs(truth))),1e-300)),
        ROI_complex_MAE=float(np.mean(np.abs(error[roi]))) if roi.any() else None,
        ROI_Re_MAE=float(np.mean(np.abs(error.real[roi]))) if roi.any() else None,
        ROI_Im_MAE=float(np.mean(np.abs(error.imag[roi]))) if roi.any() else None,
        historical_truth_relative_support_IoU=iou(historical_a,historical_t),
        fixed_increment_support_IoU=iou(online_a,online_t),
        increment_boundary=boundary_metrics(online_a,online_t),
        strict_vacuum_cells=int(np.count_nonzero(np.abs(chi)<=1e-10)),
        strict_truth_vacuum_cells=int(np.count_nonzero(np.abs(truth)<=1e-10)),
        strict_truth_vacuum_complex_MAE=float(np.mean(np.abs(error[np.abs(truth)<=1e-10]))) if np.any(np.abs(truth)<=1e-10) else None)
    for threshold in AIR_THRESHOLDS:
        air_a=np.abs(chi)<=threshold;air_t=np.abs(truth)<=threshold
        key=f"abs_air_{threshold:g}"
        data[key]=dict(air_IoU=iou(air_a,air_t),phase_changed_cells_vs_truth=int(np.count_nonzero(air_a^air_t)),
            candidate_topology=topology(~air_a),truth_topology=topology(~air_t),
            absolute_material_boundary=boundary_metrics(~air_a,~air_t),
            air_ROI_complex_MAE=float(np.mean(np.abs(error[air_t]))) if air_t.any() else None,
            air_IoU_empty_both=bool(not air_a.any() and not air_t.any()))
    if base is not None:
        base=np.asarray(base).reshape(24,24);delta=chi-base
        data['vs_qualified_base']=dict(complex_MAE=float(np.mean(np.abs(delta))),
            normalized_material_L2=float(np.sqrt(np.mean(np.abs(delta)**2))/1.2),
            phase_changed_cells=int(np.count_nonzero((np.abs(chi)<=.1)^(np.abs(base)<=.1))),
            increment_support_changed_cells=int(np.count_nonzero((np.abs(chi-background)>.18)^(np.abs(base-background)>.18))))
    return data


def flatten(row):
    m=row['metrics'];a=m['abs_air_0.1']
    return dict(case=row['case'],phase=row['phase'],role=row['role'],arm=row.get('arm'),edit_id=row.get('edit_id'),
        array_ref=row['array_ref'],**row['validation'],
        **{k:v for k,v in m.items() if not isinstance(v,dict)},
        air_0p1_IoU=a['air_IoU'],phase_cells_vs_truth=a['phase_changed_cells_vs_truth'],
        air_cells=a['candidate_topology']['air_cells'],material_components=a['candidate_topology']['material_components'],
        enclosed_air_cavities=a['candidate_topology']['enclosed_air_cavities'],
        increment_boundary_hausdorff_m=m['increment_boundary']['hausdorff_m'],
        increment_boundary_symmetric_mean_m=m['increment_boundary']['symmetric_mean_m'],
        absolute_material_boundary_hausdorff_m=a['absolute_material_boundary']['hausdorff_m'],
        absolute_material_boundary_symmetric_mean_m=a['absolute_material_boundary']['symmetric_mean_m'],
        phase_cells_vs_qualified_base=m.get('vs_qualified_base',{}).get('phase_changed_cells'),
        material_L2_vs_qualified_base=m.get('vs_qualified_base',{}).get('normalized_material_L2'),
        genuine_new_candidate=row.get('genuine_new_candidate'),invalid_to_valid=row.get('invalid_to_valid'))


def material_panel(path, title, columns, *, tau, core=None, guard=None):
    if REUSE_FIGURES and path.exists():return
    fig, axes=plt.subplots(3,len(columns),figsize=(3.4*len(columns),9.1),squeeze=False,constrained_layout=True)
    for j,col in enumerate(columns):
        chi=col['chi'].reshape(24,24)
        ims=((1+chi.real,.5,5,'viridis','Re(epsilon_r)'),(chi.imag,0,.5,'magma','Im(epsilon_r)'),
             ((np.abs(chi)<=.1).astype(float),0,1,'gray_r','Absolute near-air: |chi|<=.1'))
        for i,(value,lo,hi,cmap,label) in enumerate(ims):
            ax=axes[i,j];im=ax.imshow(value.T,origin='lower',extent=(-.5,.5,-.5,.5),vmin=lo,vmax=hi,cmap=cmap,interpolation='nearest')
            if i==0:
                res=col.get('residual');norm='' if res is None else f"\n||W(F-y)||={res:.3f}"
                mae='' if col.get('mae') is None else f"; MAE={col['mae']:.4f}"
                ax.set_title(col['label']+norm+mae,fontsize=10)
            if j==0:ax.set_ylabel(label+'\ny [m]')
            ax.set_xlabel('x [m]')
            if core is not None and j>=2:
                xy=np.linspace(-.5+1/48,.5-1/48,24)
                if np.asarray(core).any():ax.contour(xy,xy,np.asarray(core).reshape(24,24).T,levels=[.5],colors=['#00e5ff'],linewidths=.8)
                if guard is not None and np.asarray(guard).any():ax.contour(xy,xy,np.asarray(guard).reshape(24,24).T,levels=[.5],colors=['#ffcc00'],linewidths=.6)
            if j==len(columns)-1:fig.colorbar(im,ax=axes[i,:].tolist(),shrink=.65,pad=.02)
    fig.suptitle(title+f"\nFrozen tau={tau:.6f}; shared scales; truth is OFFLINE only",fontsize=13)
    fig.savefig(path,dpi=160);plt.close(fig)


def scattering_panel(path, title, observation, sigma, states):
    """All 256 complex channels, source-major; no recomputed prediction."""
    if REUSE_FIGURES and path.exists():return
    fig,axes=plt.subplots(3,1,figsize=(12,8),constrained_layout=True)
    measured=observation.T.reshape(-1);xx=np.arange(len(measured))
    axes[0].plot(xx,measured.real,color='black',lw=1.3,label='observed')
    axes[1].plot(xx,measured.imag,color='black',lw=1.3,label='observed')
    colors=['#176b87','#d97706','#0d9488']
    for color,(name,state) in zip(colors,states):
        pred=state['prediction'].T.reshape(-1)
        axes[0].plot(xx,pred.real,color=color,lw=.8,alpha=.8,label=name)
        axes[1].plot(xx,pred.imag,color=color,lw=.8,alpha=.8,label=name)
        err=np.abs(pred-measured)*np.sqrt(2)/sigma
        axes[2].plot(xx,err,color=color,lw=.8,label=name)
    axes[0].set_ylabel('Re(scattered response)');axes[1].set_ylabel('Im(scattered response)')
    axes[2].set_ylabel('sqrt(2)|F-y| / complex sigma');axes[2].set_xlabel('source-major complex channel (32 receivers x 8 sources)')
    for ax in axes:
        ax.legend(ncol=4,fontsize=9);ax.grid(alpha=.15)
    fig.suptitle(title+'\nSaved full endpoints only; channel errors are not posterior probabilities',fontsize=12)
    fig.savefig(path,dpi=160);plt.close(fig)


def main():
    global REUSE_FIGURES
    parser=argparse.ArgumentParser();parser.add_argument('--stage',type=Path,default=STAGE)
    parser.add_argument('--reuse-figures',action='store_true',help='Reuse already rendered figures only after experiment inputs are frozen; no digest checks')
    parser.add_argument('--frozen',action='store_true',help='Parent explicitly declared all online phases frozen; label final offline snapshot')
    args=parser.parse_args();REUSE_FIGURES=args.reuse_figures
    stage=args.stage.resolve();out=stage/'offline';out.mkdir(parents=True,exist_ok=True);figdir=out/'figures';figdir.mkdir(exist_ok=True)
    trutharchive=load(A30/'inputs/a29_release/results/a29/OFFLINE_TARGETS.npz')
    ids=list(trutharchive['case_id'].astype(str));truthtable={c:trutharchive['truth'][i] for i,c in enumerate(ids)}
    bgtable={c:trutharchive['background'][i] for i,c in enumerate(ids)}
    e0path=stage/'e0/E0_AUDIT.json'
    e0=json.loads(e0path.read_text()).get('cases_detail',{}) if e0path.exists() else {}
    truth_residuals={c:d['row']['truth_noise_only_norm'] for c,d in e0.items()}
    qpath=stage/'qualification/SUMMARY.json';qualification=json.loads(qpath.read_text()) if qpath.exists() else []
    minpath=stage/'qualification/MINCHANGE_SUMMARY.json'
    minentries={r['case']:r for r in json.loads(minpath.read_text())} if minpath.exists() else {}
    rows=[];gallery=[];raw_arrays={};missing=[];residual_rows=[];mechanism_rows=[];observations={}
    metadata=dict(source='cached exposed A30 OFFLINE_TARGETS',truth_usage='separate offline scoring/plots only',
        figures_reused_from_same_frozen_inputs=REUSE_FIGURES,
        metric_definitions=dict(full_MAE='uniform full-grid susceptibility error; epsilon error identical',
            complex_permittivity_NMAE='mean|epsilon_est-epsilon_truth| / mean|epsilon_truth|',
            ROI='offline |truth-background|>.25*max|truth-background|',
            fixed_increment_support='|chi-background|>.18',historical_support='truth-relative 25% increment threshold; offline only',
            absolute_near_air='|chi|<=.1 primary; .075/.125 sensitivity',strict_vacuum='|chi|<=1e-10; different from near-air',
            boundary='4-connected mask boundary; Euclidean cell-centre distances times 1/24m; empty boundary undefined/null',
            air_IoU='both masks empty yields1; reported empty flag prevents treating this as cavity reconstruction success'),
        colors=dict(epsilon_real=[.5,5],epsilon_imag=[0,.5],air=[0,1]),
        new_physics=0,new_RHS=0,solver_imports=[],probability_semantics='candidate_only; posterior mass null')

    def add(case,phase,role,path,observation,sigma,tau,base=None,**extra):
        if not path.exists():
            missing.append(dict(case=case,phase=phase,role=role,path=str(path.relative_to(stage)),status='NOT_SAVED_OR_NOT_RUN'))
            return None
        state=load(path);audit=validate_saved(state,observation,sigma,tau)
        row=dict(case=case,phase=phase,role=role,array_ref=str(path.relative_to(stage)),validation=audit,
                 metrics=metrics(state['chi'],truthtable[case],bgtable[case],base=base),**extra)
        rows.append(row);prefix=f"{phase}__{case}__{role}__{extra.get('edit_id','none')}__{extra.get('arm','none')}"
        raw_arrays[prefix+'__chi']=state['chi'];raw_arrays[prefix+'__prediction']=state['prediction']
        residual_rows.append(dict(case=case,phase=phase,role=role,arm=extra.get('arm'),edit_id=extra.get('edit_id'),
                                  residual_norm=audit['residual_norm'],tau=tau))
        return state,row

    for entry in qualification:
        case=entry['case'];folder=stage/'qualification'/case
        ob=load(folder/'observation.npz');y=ob['observation'];sigma=float(ob['sigma']);observations[case]=(y,sigma)
        tau=float(entry['initial']['tau']);basepath=folder/'qualified_base.npz';base=load(basepath)['chi'] if basepath.exists() else None
        raw_arrays[case+'__OFFLINE_truth']=truthtable[case];raw_arrays[case+'__background']=bgtable[case]
        raw_arrays[case+'__observation']=y
        columns=[dict(label='OFFLINE truth',chi=truthtable[case],residual=truth_residuals.get(case),mae=0.)]
        baseroles=[('original_N3',folder/'original_n3.npz',None),('qualified_base',basepath,entry.get('selected_arm'))]
        if case in minentries:
            baseroles.append(('qualified_minchange_base',folder/'qualified_minchange_base.npz',minentries[case].get('selected_arm')))
        for role,pth,arm in baseroles:
            result=add(case,'qualification',role,pth,y,sigma,tau,base=base,arm=arm)
            if result is not None:
                st,row=result;columns.append(dict(label=role+(f' ({arm})' if arm else ''),chi=st['chi'],residual=row['validation']['residual_norm'],mae=row['metrics']['full_complex_MAE']))
        for arm in entry['arms']:
            armname=arm['arm'];best=folder/armname/'best_valid.npz';final=folder/armname/'final.npz'
            result=add(case,'qualification',f'{armname}_final',final,y,sigma,tau,base=base,arm=armname)
            if best.exists():add(case,'qualification',f'{armname}_actual_best_valid',best,y,sigma,tau,base=base,arm=armname)
            if result is not None:
                st,row=result;columns.append(dict(label=armname,chi=st['chi'],residual=row['validation']['residual_norm'],mae=row['metrics']['full_complex_MAE']))
        fpath=figdir/f'qualification_{case}.png';material_panel(fpath,f'{case}: base qualification',columns,tau=tau)
        gallery.append(dict(kind='qualification',case=case,image=str(fpath.relative_to(out)),qualified=entry['qualified'],selected_arm=entry.get('selected_arm')))

    phasepaths=sorted(stage.glob('mechanism*/SUMMARY.json'))
    for mpath in phasepaths:
        phase=mpath.parent.name;entries=json.loads(mpath.read_text())
        for entry in entries:
            case=entry['case'];edit=entry['edit_id'];folder=mpath.parent/case/f'edit_{edit:02}'
            y,sigma=observations[case];tau=float(entry['base_validation']['tau']);basepath=mpath.parent/case/'base.npz'
            base=load(basepath);geometry=load(folder/'edit_geometry.npz') if (folder/'edit_geometry.npz').exists() else {}
            original_base=load(stage/'qualification'/case/'qualified_base.npz')['chi']
            minchange_path=stage/'qualification'/case/'qualified_minchange_base.npz'
            minchange_base=load(minchange_path)['chi'] if minchange_path.exists() else None
            original_n3=load(stage/'qualification'/case/'original_n3.npz')['chi']
            is_original=bool(np.array_equal(base['chi'],original_base))
            is_minchange=bool(np.array_equal(base['chi'],minchange_base)) if minchange_base is not None else False
            base_policy='SAME_AS_ORIGINAL_AND_MINCHANGE' if is_original and is_minchange else ('MINCHANGE_FROM_N3' if is_minchange else ('MINIMUM_OBSERVED_RESIDUAL' if is_original else 'OTHER_DECLARED_BASE'))
            cols=[dict(label='OFFLINE truth',chi=truthtable[case],residual=truth_residuals.get(case),mae=0.),
                  dict(label='qualified base',chi=base['chi'],residual=entry['base_validation']['residual_norm'],
                       mae=metrics(base['chi'],truthtable[case],bgtable[case])['full_complex_MAE'])]
            for role,pth in [('raw',folder/'raw.npz'),('final',folder/'final.npz'),('actual_best_valid',folder/'best_valid.npz')]:
                result=add(case,phase,role,pth,y,sigma,tau,base=base['chi'],edit_id=edit,
                           genuine_new_candidate=entry.get('genuine_new_candidate'),invalid_to_valid=entry.get('invalid_to_valid'))
                if result is not None:
                    st,row=result
                    if 'core' in geometry:
                        outside=~geometry['core'].reshape(24,24).astype(bool)
                        a=np.abs(st['chi'].reshape(24,24))<=.1;b=np.abs(base['chi'].reshape(24,24))<=.1
                        row['metrics']['phase_changed_cells_outside_core']=int(np.count_nonzero((a^b)&outside))
                        row['metrics']['newly_near_air_cells_outside_core']=int(np.count_nonzero(a&~b&outside))
                        row['metrics']['erased_base_near_air_cells_outside_core']=int(np.count_nonzero(~a&b&outside))
                    if role=='raw' or (role=='actual_best_valid' or (role=='final' and not (folder/'best_valid.npz').exists())):
                        cols.append(dict(label=role,chi=st['chi'],residual=row['validation']['residual_norm'],mae=row['metrics']['full_complex_MAE']))
            fpath=figdir/f'{phase}_{case}_edit{edit:02}.png'
            material_panel(fpath,f"{case}, edit {edit}: raw vs bounded compensation",cols,tau=tau,core=geometry.get('core'),guard=geometry.get('guard'))
            if edit==0 and (folder/'raw.npz').exists():
                finalpath=folder/'best_valid.npz' if (folder/'best_valid.npz').exists() else folder/'final.npz'
                if finalpath.exists():
                    scattering_panel(figdir/f'scattering_fit_{phase}_{case}_edit{edit:02}.png',
                        f'{phase} / {case} / edit{edit}',y,sigma,
                        [('qualified base',base),('raw edit',load(folder/'raw.npz')),('actual corrected endpoint',load(finalpath))])
            gallery.append(dict(kind=phase,case=case,edit_id=edit,image=str(fpath.relative_to(out)),
                genuine_new_candidate=entry.get('genuine_new_candidate'),invalid_to_valid=entry.get('invalid_to_valid')))
            mechanism_rows.append(dict(phase=phase,case=case,edit_id=edit,raw_residual=entry['raw_validation']['residual_norm'],
                base_policy=base_policy,base_matches_original_selected=is_original,base_matches_minchange=is_minchange,
                base_normalized_L2_from_original_N3=float(np.sqrt(np.mean(np.abs(base['chi']-original_n3)**2))/1.2),
                corrected_residual=entry['corrected_validation']['residual_norm'],tau=tau,
                raw_status=entry['raw_validation']['status'],corrected_status=entry['corrected_validation']['status'],
                genuine_new_candidate=entry.get('genuine_new_candidate'),invalid_to_valid=entry.get('invalid_to_valid'),
                cost=entry.get('cost'),material_distance=entry.get('material_distance')))
    if mechanism_rows:
        fig,axes=plt.subplots(1,2,figsize=(13,4.5),constrained_layout=True)
        labels=[f"{r['case'].split('_')[0][1:]}/{r['phase'].replace('mechanism','M')}/e{r['edit_id']}" for r in mechanism_rows]
        xx=np.arange(len(mechanism_rows));raw=np.array([r['raw_residual'] for r in mechanism_rows]);fin=np.array([r['corrected_residual'] for r in mechanism_rows]);tau=mechanism_rows[0]['tau']
        axes[0].plot(xx,raw,'o-',label='raw edit');axes[0].plot(xx,fin,'s-',label='corrected actual endpoint');axes[0].axhline(tau,c='black',ls='--',label='frozen tau')
        axes[0].set_xticks(xx,labels,rotation=75,ha='right',fontsize=7);axes[0].set_ylabel('Whitened residual to OBSERVED data');axes[0].legend()
        cells=[r.get('material_distance',{}).get('phase_changed_cells',0) for r in mechanism_rows]
        axes[1].bar(xx,cells,color=['#0d9488' if r.get('genuine_new_candidate') else '#94a3b8' for r in mechanism_rows]);axes[1].set_ylabel('Absolute near-air/material phase cells changed vs base');axes[1].set_xticks(xx,labels,rotation=75,ha='right',fontsize=7)
        fig.savefig(figdir/'mechanism_residual_and_structure.png',dpi=160);plt.close(fig)
    matched=[]
    phase_names=sorted({r['phase'] for r in mechanism_rows})
    if len(phase_names)>1:
        reference_phase=phase_names[0]
        for later in phase_names[1:]:
            old={(r['case'],r['edit_id']):r for r in mechanism_rows if r['phase']==reference_phase}
            for event in [r for r in mechanism_rows if r['phase']==later]:
                key=(event['case'],event['edit_id'])
                if key not in old:continue
                prior=old[key]
                def endpoint(phase):
                    candidates=[r for r in rows if r['phase']==phase and r['case']==key[0] and r.get('edit_id')==key[1] and r['role'] in ('final','actual_best_valid')]
                    return next((r for r in candidates if r['role']=='actual_best_valid'),next((r for r in candidates if r['role']=='final'),None))
                a,b=endpoint(reference_phase),endpoint(later)
                if a is None or b is None:continue
                matched.append(dict(case=key[0],edit_id=key[1],baseline_phase=reference_phase,revised_phase=later,
                    baseline_base_policy=prior.get('base_policy'),revised_base_policy=event.get('base_policy'),
                    baseline_raw_residual=prior['raw_residual'],revised_raw_residual=event['raw_residual'],
                    raw_residual_difference=event['raw_residual']-prior['raw_residual'],
                    baseline_corrected_residual=a['validation']['residual_norm'],revised_corrected_residual=b['validation']['residual_norm'],
                    baseline_status=a['validation']['validation_status'],revised_status=b['validation']['validation_status'],
                    baseline_complex_MAE=a['metrics']['full_complex_MAE'],revised_complex_MAE=b['metrics']['full_complex_MAE'],
                    baseline_phase_changes_outside_core=a['metrics'].get('phase_changed_cells_outside_core'),
                    revised_phase_changes_outside_core=b['metrics'].get('phase_changed_cells_outside_core'),
                    baseline_new_near_air_outside_core=a['metrics'].get('newly_near_air_cells_outside_core'),
                    revised_new_near_air_outside_core=b['metrics'].get('newly_near_air_cells_outside_core'),
                    baseline_genuine=prior.get('genuine_new_candidate'),revised_genuine=event.get('genuine_new_candidate'),
                    baseline_cost=prior.get('cost'),revised_cost=event.get('cost')))
    csvwrite(out/'MATCHED_ITERATION_COMPARISON.csv',matched)
    dump(out/'MATCHED_ITERATION_COMPARISON.json',matched)
    flat=[flatten(r) for r in rows];csvwrite(out/'MATERIAL_METRICS.csv',flat);csvwrite(out/'RESIDUALS.csv',residual_rows)
    csvwrite(out/'MECHANISM_OFFLINE_CONTEXT.csv',mechanism_rows)
    dump(out/'METRICS.json',dict(metadata=metadata,records=rows,missing=missing))
    np.savez_compressed(out/'OFFLINE_MATERIAL_AND_PREDICTION_ARRAYS.npz',**raw_arrays)
    selected=[r for r in rows if r['phase']=='qualification' and r['role']=='qualified_base']
    minselected=[r for r in rows if r['phase']=='qualification' and r['role']=='qualified_minchange_base']
    n3=[r for r in rows if r['phase']=='qualification' and r['role']=='original_N3']
    group=lambda rr,key:float(np.mean([r['metrics'][key] for r in rr])) if rr else None
    summary=dict(qualification_cases_read=len(qualification),qualification_bases_scored=len(selected),
        qualification_original_N3_mean_complex_MAE=group(n3,'full_complex_MAE'),qualification_selected_mean_complex_MAE=group(selected,'full_complex_MAE'),
        qualification_minchange_mean_complex_MAE=group(minselected,'full_complex_MAE'),
        minchange_bases_scored=len(minselected),
        minchange_policy_status='Post hoc algorithm revision motivated by exposed offline aggregate; each per-case choice uses N3 displacement and data legality only, not truth MAE',
        mechanism_completed_rows=len(mechanism_rows),mechanism_phase_files=[str(p.relative_to(stage)) for p in phasepaths],
        genuine_candidates=sum(bool(r['genuine_new_candidate']) for r in mechanism_rows),
        invalid_to_valid=sum(bool(r['invalid_to_valid']) for r in mechanism_rows),
        mechanism_by_phase={p:dict(completed=sum(r['phase']==p for r in mechanism_rows),
            genuine=sum(bool(r['genuine_new_candidate']) for r in mechanism_rows if r['phase']==p),
            invalid_to_valid=sum(bool(r['invalid_to_valid']) for r in mechanism_rows if r['phase']==p)) for p in phase_names},
        matched_iteration_pairs=len(matched),
        new_full_state_attempts=0,new_RHS=0,probabilities=None,truth_access='OFFLINE_ONLY',
        snapshot_status='FINAL_FROZEN_ONLINE_PHASES_OFFLINE_SCORE' if args.frozen else 'COMPLETED_ROWS_ONLY; rerun after parent freezes all mechanism iterations',
        note='Offline MAE does not select, rank or initialize online candidates; box/data legality remains independent of truth scoring')
    dump(out/'SUMMARY.json',summary);dump(out/'GALLERY_DATA.json',dict(summary=summary,items=gallery))
    write_statistics(out,summary,rows,mechanism_rows,matched)
    build_html(out,summary,gallery,flat)
    receipt=dict(status='COMPLETE_OFFLINE_SNAPSHOT',pid=os.getpid(),new_full_state_attempts=0,new_complex_full_solver_RHS=0,
        solver_modules_imported=[],GPU_seconds=0,old_results_modified=False,truth_usage='existing exposed offline A30 only',
        wall_seconds=perf_counter()-START,process_cpu_seconds=process_time()-CPU,qualification_rows=len(qualification),mechanism_rows=len(mechanism_rows))
    with (out/'EXECUTION_RECEIPTS.jsonl').open('a') as f:f.write(json.dumps(receipt)+'\n')
    dump(out/'EXECUTION_RECEIPT.json',receipt)
    print(json.dumps(summary))


def write_statistics(out,summary,rows,mechanisms,matched):
    mean=lambda values: float(np.mean(values)) if values else None
    qtable=[]
    for role in ('original_N3','full_final','multiscale_haar_final','qualified_base','qualified_minchange_base'):
        rr=[r for r in rows if r['phase']=='qualification' and r['role']==role]
        if not rr:continue
        qtable.append(dict(role=role,cases=len(rr),valid=sum(r['validation']['validation_status']=='FULL_VALIDATED_NUMERICAL' for r in rr),
            observed_residual=mean([r['validation']['residual_norm'] for r in rr]),
            complex_MAE=mean([r['metrics']['full_complex_MAE'] for r in rr]),
            Re_MAE=mean([r['metrics']['full_Re_MAE'] for r in rr]),Im_MAE=mean([r['metrics']['full_Im_MAE'] for r in rr]),
            near_air_phase_cells_vs_truth=mean([r['metrics']['abs_air_0.1']['phase_changed_cells_vs_truth'] for r in rr])))
    mtable=[]
    for phase in sorted({r['phase'] for r in mechanisms}):
        mm=[r for r in mechanisms if r['phase']==phase];endpoints=[]
        for event in mm:
            rr=[r for r in rows if r['phase']==phase and r['case']==event['case'] and r.get('edit_id')==event['edit_id']]
            endpoint=next((r for r in rr if r['role']=='actual_best_valid'),next(r for r in rr if r['role']=='final'))
            endpoints.append(endpoint)
        mtable.append(dict(phase=phase,edits=len(mm),parents=len({r['case'].split('_')[0] for r in mm}),
            new_valid=sum(bool(r['genuine_new_candidate']) for r in mm),invalid_to_valid=sum(bool(r['invalid_to_valid']) for r in mm),
            raw_residual_mean=mean([r['raw_residual'] for r in mm]),
            actual_endpoint_residual_mean=mean([r['validation']['residual_norm'] for r in endpoints]),
            actual_endpoint_complex_MAE_all=mean([r['metrics']['full_complex_MAE'] for r in endpoints]),
            actual_endpoint_Re_MAE_all=mean([r['metrics']['full_Re_MAE'] for r in endpoints]),
            actual_endpoint_Im_MAE_all=mean([r['metrics']['full_Im_MAE'] for r in endpoints]),
            actual_endpoint_complex_MAE_valid_only=mean([r['metrics']['full_complex_MAE'] for r in endpoints if r['validation']['validation_status']=='FULL_VALIDATED_NUMERICAL']),
            outside_core_phase_changes_mean=mean([r['metrics']['phase_changed_cells_outside_core'] for r in endpoints if 'phase_changed_cells_outside_core' in r['metrics']]),
            newly_near_air_outside_core_mean=mean([r['metrics']['newly_near_air_cells_outside_core'] for r in endpoints if 'newly_near_air_cells_outside_core' in r['metrics']]),
            base_policies=sorted({r.get('base_policy','UNKNOWN') for r in mm})))
    dump(out/'AGGREGATE_METRICS.json',dict(qualification=qtable,mechanism=mtable,matched_pairs=len(matched),
        interpretation='Exposed pilot descriptive means, not independent confirmatory validation. Different phases/base policies kept separate.'))
    csvwrite(out/'QUALIFICATION_AGGREGATES.csv',qtable);csvwrite(out/'MECHANISM_AGGREGATES.csv',mtable)
    qlines='\n'.join(f"|{r['role']}|{r['valid']}/{r['cases']}|{r['observed_residual']:.6f}|{r['complex_MAE']:.6f}|{r['Re_MAE']:.6f}|{r['Im_MAE']:.6f}|" for r in qtable)
    mlines='\n'.join(f"|{r['phase']}|{r['new_valid']}/{r['edits']}|{r['invalid_to_valid']}|{r['raw_residual_mean']:.6f}|{r['actual_endpoint_residual_mean']:.6f}|{r['actual_endpoint_complex_MAE_all']:.6f}|{r['outside_core_phase_changes_mean']:.2f}|" for r in mtable)
    text=f'''# Research T：最终离线材料评分

状态：{summary['snapshot_status']}。读取原已暴露A30真值，未读取locked truth，未运行新forward/adjoint/tangent；新物理/RHS均为0。所有MAE只在独立离线进程产生，不进入候选选择或在线优化。字段、数组与图像可复现于本目录。

## 基础解：数据相容不等于材料准确

|规则/端点|数据与数值合格|平均白化残差|complex MAE|Re MAE|Im MAE|
|---|---:|---:|---:|---:|---:|
{qlines}

原最低观测残差选择改善相容性，但明显增加full-grid材料误差。随后最小N3位移规则是受这轮exposed离线汇总启发的**post-hoc算法修订**，不是预注册盲验证。该规则在每个case仅用已合格端点与N3间的材料距离，不用真值MAE选择。它减少了材料偏离，仍不能宣称比原N3更接近真值。

## 有限编辑与补偿：逐轮保留

|轮次|显著合格编辑|不合格→合格|raw平均残差|实际选中端点平均残差|全部端点complex MAE|core外平均相变化细胞|
|---|---:|---:|---:|---:|---:|---:|
{mlines}

“实际端点”优先使用保存的best_valid；没有合法端点时保留final失败端点。全部端点均值没有悄悄去除失败；合法子集均值另在AGGREGATE_METRICS.json。raw→corrected通过数是物理相容与保护结构的机制见证，**不是**重建真值、posterior mass、独立材料族数、或搜索方法优越性。

v1与v2采用同一最低残差base和相同raw编辑；v3换成minchange base，只比较固定编辑indices0/2，必须与原base分开。MATCHED_ITERATION_COMPARISON保留24个配对，包含base身份及raw residual差。跨轮32个通过记录不能直接当32个已去重材料世界。

## 图、掩膜和边界

所有epsilon实部共用[.5,5]、损耗共用[0,.5]。原增量support与绝对近空气掩膜分开：primary近空气为|chi|≤.1，另列.075/.125敏感性，严格真空为|chi|≤1e-10，不能把这些定义混用。空气IoU在两空集时取1，并附empty_both标志，不能因此宣称发现空腔。边界距离是24×24像素中心栅格量，不是连续几何或网格收敛证书；任一边界空集时为null。

图库index.html含8case基础反演、40条机制编辑、实际best_valid、全8源散射拟合、保护core/guard与材料差别；原始full-material/prediction保存在OFFLINE_MATERIAL_AND_PREDICTION_ARRAYS.npz。真值图仅离线展示。

## 费用与失败

物理求解费用由父阶段唯一ledger计费，离线脚本不写该ledger。独立零物理后处理费用存EXECUTION_RECEIPTS.jsonl。首次最终刷新在图已生成后出现重复字段TypeError，FAILED_REFRESH_RECEIPT.json保留异常；其wall未被外层计时，明示null/UNRECORDED，不冒写0。修复后复用同一冻结输入已生成的图，重算全部指标，无新物理。所有概率字段为null/candidate_only。
'''
    (out/'OFFLINE_SCORING_REPORT_ZH.md').write_text(text,encoding='utf-8')


def build_html(out,summary,gallery,metricsrows):
    cards=[]
    for g in gallery:
        title=g['case']+' / '+g['kind']+(f" / edit {g['edit_id']}" if 'edit_id' in g else '')
        note=('Qualified='+str(g['qualified'])+'; selected='+str(g['selected_arm'])) if g['kind']=='qualification' else ('new significant candidate='+str(g['genuine_new_candidate'])+'; invalid→valid='+str(g['invalid_to_valid']))
        cards.append(f'<article><h2>{html.escape(title)}</h2><p>{html.escape(note)}</p><a href="{g["image"]}"><img loading="lazy" src="{g["image"]}" alt="Offline material scoring"></a></article>')
    table=[]
    for r in metricsrows:
        table.append('<tr>'+''.join('<td>'+html.escape(str(x))+'</td>' for x in (r['case'],r['phase'],r['role'],r.get('edit_id'),f"{r['residual_norm']:.4f}",r['validation_status'],f"{r['full_complex_MAE']:.5f}",f"{r['full_Re_MAE']:.5f}",f"{r['full_Im_MAE']:.5f}",r['phase_cells_vs_truth']))+'</tr>')
    body=f'''<!doctype html><html lang="zh"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Research T · 离线材料评分</title><style>body{{font:16px system-ui,sans-serif;margin:30px auto;max-width:1500px;padding:0 18px;color:#182230;background:#f6f8fb}}h1,h2{{color:#16364c}}article{{background:white;padding:18px;margin:22px 0;border-radius:14px}}img{{width:100%;height:auto}}.note{{padding:18px;background:#e5f1f3;border-left:4px solid #0d9488}}.scroll{{overflow:auto}}table{{border-collapse:collapse;background:white;font-size:12px}}td,th{{padding:7px;border-bottom:1px solid #ddd;white-space:nowrap}}a{{color:#176b87}}</style><h1>Research T：材料相、基础解与补偿端点</h1><div class="note">严格离线评分：真值没有进入proposal、corrector或在线排序。所有图共用epsilon实部[0.5,5]、损耗[0,0.5]色标。近空气阈值|chi|≤0.1与严格真空χ=0不同；原暗增量图也不等于空气。只展示摘要已保存完成的记录，未完成阶段不记失败。</div><p>基础解评分数：{summary['qualification_bases_scored']}；机制已完成编辑记录：{summary['mechanism_completed_rows']}；显著新候选记录：{summary['genuine_candidates']}；补偿前失败而补偿后合格记录：{summary['invalid_to_valid']}。跨迭代记录不宣称去重材料族数。</p><p><a href="METRICS.json">完整指标JSON</a> · <a href="MATERIAL_METRICS.csv">逐材料CSV</a> · <a href="RESIDUALS.csv">真实观测残差CSV</a> · <a href="MATCHED_ITERATION_COMPARISON.json">相同编辑跨轮比较</a> · <a href="OFFLINE_MATERIAL_AND_PREDICTION_ARRAYS.npz">原始材料/预测数组</a></p><p>下图空腔core青色、guard黄色；它们来自冻结truth-free编辑，不来自真值。MAE只是离线诊断。局部白化残差变小不能替代冻结相容阈值与结构合法性。</p>'''
    if summary['mechanism_completed_rows']:body+='<article><h2>补偿残差与结构变化</h2><img src="figures/mechanism_residual_and_structure.png"></article>'
    for fpath in sorted((out/'figures').glob('scattering_fit_*.png')):
        rel=str(fpath.relative_to(out));body+=f'<article><h2>全部8源散射拟合：{html.escape(fpath.stem)}</h2><img loading="lazy" src="{rel}"></article>'
    body+=''.join(cards)+'<h2>逐端点指标</h2><div class="scroll"><table><tr>'+''.join('<th>'+s+'</th>' for s in ['case','phase','role','edit','whitened residual','validation','complex MAE','Re MAE','Im MAE','phase cells vs truth'])+'</tr>'+''.join(table)+'</table></div><p>candidate_only，概率为null。相容材料不等于连续分支、posterior mode或完整覆盖。</p></html>'
    (out/'index.html').write_text(body,encoding='utf-8')


if __name__=='__main__':main()
