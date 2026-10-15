"""Bounded independent physical-coordinate integration of frozen histories.

This diagnostic wraps the existing Reference.adaptive_physical 5/9 rule. It
does not change policies, sampling, primary statistics, legacy audit statuses,
the underlying continuous proxy, prior measure or persistent gain model.
"""
from __future__ import annotations

import csv
import gc
import json
import os
import resource
import sys
import time
from pathlib import Path

for variable in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[variable]='1'
import numpy as np

REPO=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(REPO))
from src.b1_joint.response import Response
from src.b1_joint.reference import Reference
from src.b1_joint.material import RMAX
from cost import Cost,ROOT,dump,BudgetExceeded

WIDTHS=np.array([*RMAX,.6,.12])


class DiagnosticDeadline(RuntimeError):
    pass


def _peak_bytes():
    peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(peak if sys.platform=='darwin' else peak*1024)


def _compress(blocks):
    """Sufficient repeat means preserve a single gain, not per-block gains."""
    groups={}
    for block in blocks:
        groups.setdefault(block['frequency'],[]).append(block)
    compressed=[];constant=0.
    for frequency,items in groups.items():
        y=np.array([item['y'] for item in items],float)
        sd=np.array([item['sigma'] for item in items],float)
        precision=1/(sd*sd)
        mean=np.sum(precision*y,axis=0)/precision.sum(axis=0)
        effective_sd=1/np.sqrt(precision.sum(axis=0))
        # Product of full block densities = compressed density * exp(constant).
        # The factor is independent of material AND the persistent gain.
        constant+=float(-.5*np.sum(((y-mean)/sd)**2)-np.log(sd).sum()
            +np.log(effective_sd).sum()-.5*(len(items)-1)*y.shape[1]*np.log(2*np.pi))
        compressed.append(dict(frequency=frequency,y=mean,sigma=effective_sd,
            original_actions=[item['action'] for item in items],repeated_blocks=len(items)))
    return compressed,constant


class HistoryResponse:
    """Reference's F1 is the entire actually selected compressed history."""
    def __init__(self,response,blocks,cost,deadline):
        self.response,self.blocks=response,blocks
        self.cost,self.deadline=cost,deadline
        self.bounds=response.bounds
        self.mesh_calls=0;self.proxy_material_frequency_actions=0

    def check(self):
        self.cost.check_time()
        if time.perf_counter()>=self.deadline:
            raise DiagnosticDeadline('Frozen diagnostic computation deadline reached; reserve time for export')

    def mesh(self,axes,k=1.6):
        self.check()
        count=int(np.prod([len(axis) for axis in axes]))
        if float(k)==2.4:
            # Reference.cell computes a second cache even for condition C. It
            # is unused; keep the number of rows and zero columns explicit.
            return np.empty((count,0),dtype=float)
        if float(k)!=1.6:
            raise ValueError('HistoryResponse only exposes full joint F1 to condition C')
        self.mesh_calls+=1
        fields=[]
        for block in self.blocks:
            self.check()
            frequency=block['frequency']
            if frequency is None:
                fields.append(np.ones((count,1)))
            else:
                fields.append(self.response.mesh(axes,float(frequency)))
                self.proxy_material_frequency_actions+=count
        return np.concatenate(fields,axis=1)


def _difference(estimate,target):
    mean=np.asarray(estimate['mean']);sd=np.asarray(estimate['SD'])
    other_mean=np.asarray(target['mean']);other_sd=np.asarray(target['SD'])
    return dict(event_mass_max=float(np.max(abs(np.asarray(estimate['event_mass'])-target['event_mass']))),
        mean_fraction_prior_width=float(np.max(abs(mean-other_mean)/WIDTHS)),
        mean_abs=abs(mean-other_mean),
        SD_relative=float(np.max(abs(sd-other_sd)/np.maximum(other_sd,WIDTHS*1e-6))),
        log_evidence_abs=float(abs(estimate['log_evidence']-target['log_evidence'])))


def _write(results,report,final=False):
    dump(results/('HISTORY_ADAPTIVE_AUDIT.json' if final else 'HISTORY_ADAPTIVE_AUDIT_PARTIAL.json'),report)
    rows=[]
    for row in report['rows']:
        adaptive=row.get('adaptive',{})
        rows.append(dict(history=row['history'],status=row['status'],legacy_13_17_status=row['legacy_13_17_status'],
            converged=adaptive.get('converged'),cells=adaptive.get('cells'),heuristic_relative_error=adaptive.get('estimated_relative_error'),
            event_difference_vs13=row.get('difference_vs13',{}).get('event_mass_max'),
            event_difference_vs17=row.get('difference_vs17',{}).get('event_mass_max'),
            event_difference_vs_finite2048=row.get('difference_vs_finite2048',{}).get('event_mass_max'),
            mean_fraction_difference_vs13=row.get('difference_vs13',{}).get('mean_fraction_prior_width'),
            mean_fraction_difference_vs17=row.get('difference_vs17',{}).get('mean_fraction_prior_width'),
            mean_fraction_difference_vs_finite2048=row.get('difference_vs_finite2048',{}).get('mean_fraction_prior_width'),
            wall_seconds=row.get('wall_seconds')))
    if rows:
        with (results/'HISTORY_ADAPTIVE_AUDIT.csv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def run(cost):
    root=ROOT;results=root/'results'
    if (results/'HISTORY_ADAPTIVE_AUDIT.json').exists():
        raise ValueError('Existing independent audit is preserved; no automatic repeat')
    cfg=json.loads((root/'configs/HISTORY_ADAPTIVE_AUDIT_FROZEN.json').read_text())
    if cfg['relative']!=2e-6 or cfg['max_cells']!=512 or cfg['total_additional_wall_seconds_cap']!=300:
        raise ValueError('Independent diagnostic must retain its frozen scope/tolerances')
    original=json.loads((results/'HISTORY_REFERENCE_AUDIT.json').read_text())
    if set(cfg['fixed_unique_histories'])!=set(original['unique_history_details']):
        raise ValueError('Frozen history set differs from original fixed16 audit')
    response=Response(results/'E1_RESPONSE_BANK.npz')
    deadline=cost.start+float(cfg['computation_deadline_seconds'])
    report=dict(schema='research-b-independent-history-adaptive-v1',frozen_config=cfg,
        new_Maxwell=0,GPU=0,old_results_modified=False,old_UNRESOLVED_statuses_preserved=True,
        independent_parent_count=4,unique_histories=16,policy_labels_in_original=32,
        target='same original continuous coefficient prior, same saved-response proxy, one persistent proper truncated gain',
        method='Reference.adaptive_physical existing 5/9 physical-coordinate rule with unchanged prior Jacobian and atoms',
        certified=False,relative=2e-6,max_cells=512,
        probabilities_use_compressed_history_with_same_gain=True,
        original_log_evidence_constant_restored=True,
        risk_not_integrated='existing adaptive feature vector evaluates events/moments, not selected-material L1 risk',rows=[])
    stopped=False
    for key in cfg['fixed_unique_histories']:
        source=original['unique_history_details'][key]
        row=dict(history=key,legacy_13_17_status=source['status'],legacy_status_unchanged=True,
            status='UNRESOLVED_NOT_RUN_WALL_LIMIT')
        if stopped or time.perf_counter()>=deadline:
            stopped=True;report['rows'].append(row);continue
        begin=time.perf_counter();wrapper=None
        try:
            blocks,constant=_compress(source['selected_observation_history'])
            wrapper=HistoryResponse(response,blocks,cost,deadline)
            reference=Reference(wrapper)
            y=np.concatenate([block['y'] for block in blocks]);sd=np.concatenate([block['sigma'] for block in blocks])
            baseline=reference.solve(y,sd,order=13,condition='C')
            baseline_summary=dict(event_mass=baseline['event_mass'],mean=baseline['mean'],SD=baseline['SD'],
                                  log_evidence=float(baseline['log_evidence']+constant))
            row['compression_identity_difference_vs_original13']=_difference(baseline_summary,source['order13'])
            row['compressed_history']=blocks;row['theta_independent_log_constant']=constant
            light={name:baseline[name] for name in ('known_eps','condition','loglik','context','y','sigma')}
            reference.cache.clear();del baseline;gc.collect()
            target_difference=row['compression_identity_difference_vs_original13']
            if target_difference['event_mass_max']>1e-8 or target_difference['log_evidence_abs']>1e-8:
                row['status']='UNRESOLVED_HISTORY_TARGET_MISMATCH'
            else:
                adaptive=reference.adaptive_physical(light,max_cells=cfg['max_cells'],
                    relative=cfg['relative'],prior_measure='coefficient_uniform')
                adaptive['compressed_log_evidence']=adaptive['log_evidence']
                adaptive['log_evidence']=float(adaptive['log_evidence']+constant)
                row['adaptive']=adaptive
                row['status']='EMPIRICAL_ADAPTIVE_CONVERGENCE_NO_CERTIFICATE' if adaptive['converged'] else 'UNRESOLVED_MAX_CELLS'
                row['difference_vs13']=_difference(adaptive,source['order13'])
                row['difference_vs17']=_difference(adaptive,source['order17'])
                row['difference_vs_finite2048']=_difference(adaptive,source['finite2048']['summary'])
            del reference,light
        except (DiagnosticDeadline,BudgetExceeded) as error:
            row['status']='UNRESOLVED_WALL_INTERRUPTED'
            row['error']=type(error).__name__+': '+str(error)
            stopped=True
        except Exception as error:
            row['status']='UNRESOLVED_NUMERICAL_EXCEPTION'
            row['error']=type(error).__name__+': '+str(error)
        row['wall_seconds']=time.perf_counter()-begin
        if wrapper is not None:
            row['mesh_calls']=wrapper.mesh_calls
            row['proxy_material_frequency_actions']=wrapper.proxy_material_frequency_actions
        report['rows'].append(row)
        report['peak_RSS_bytes']=_peak_bytes()
        _write(results,report)
        print(json.dumps(dict(history=key,status=row['status'],
            cells=row.get('adaptive',{}).get('cells'),heuristic=row.get('adaptive',{}).get('estimated_relative_error'),
            elapsed=time.perf_counter()-cost.start,new_Maxwell=0)),flush=True)
        gc.collect()
    report['elapsed_before_final_export']=time.perf_counter()-cost.start
    return finalize_saved(results,report)


def finalize_saved(results,report):
    """Export saved numerical evidence; never runs an integration or solver."""
    original=json.loads((results/'HISTORY_REFERENCE_AUDIT.json').read_text())
    rows=report['rows']
    complete=[row for row in rows if all(name in row for name in
        ('adaptive','difference_vs13','difference_vs17','difference_vs_finite2048'))]
    stopped=any('WALL' in row['status'] for row in rows)
    report['execution_status']='BUDGET_STOPPED_WITH_UNFINISHED_ITEMS' if stopped else 'BOUNDED_DIAGNOSTIC_COMPLETED'
    report['adaptive_converged_histories']=[row['history'] for row in rows if row.get('adaptive',{}).get('converged')]
    report['unresolved_or_unfinished_histories']=[row['history'] for row in rows if row['status'].startswith('UNRESOLVED')]
    report['legacy_UNRESOLVED_histories_unchanged']=original['unresolved_unique_histories']
    report['summary']={}
    for target in ('13','17','finite2048'):
        for metric in ('event_mass_max','mean_fraction_prior_width','SD_relative'):
            suffix='_finite2048' if target=='finite2048' else target
            values=[row[f'difference_vs{suffix}'][metric] for row in complete]
            report['summary'][f'maximum_{metric}_difference_vs{target}']=max(values) if values else None
    report['completed_adaptive_evaluations']=len(complete)
    report['peak_RSS_bytes']=max(int(report.get('peak_RSS_bytes',0)),_peak_bytes())
    report['limits']='Heuristic stopping only; no new full-Maxwell calls, order21, CDF, all64-case campaign, algorithm update, primary-table update, or old status promotion.'
    _write(results,report,final=True)
    short=['# 独立physical-coordinate history审计','',
        f'固定16唯一history；完成adaptive估计{len(complete)}，满足原heuristic relative准则{len(report["adaptive_converged_histories"])}，未稳/未完成{len(report["unresolved_or_unfinished_histories"])}。',
        '原13/17审计、全部旧UNRESOLVED、算法与主MAE/配对表均保持原样。','',
        f'总阶段状态：{report["execution_status"]}。max_cells512、relative2e-6未放宽。',
        f'新增Maxwell0、GPU0；最终计费以COMPUTE_LEDGER中本job receipt为准。','']
    for key,value in report['summary'].items():short.append(f'- {key}: {value}')
    short+=['','这些差异只属于同一连续proxy target的独立数值诊断。即使heuristic收敛，也不是certificate或full posterior calibration；未重选policy或更新reference-based风险结论。','']
    (results/'HISTORY_ADAPTIVE_AUDIT_SHORT.md').write_text('\n'.join(short))
    return report


if __name__=='__main__':
    with Cost().job('e2-history-adaptive-physical-audit-v1','E2') as cost:
        run(cost)
