"""Prepare public-only independent11arm24case and exposed self-only16 manifests.

No solver execution, array/truth reads, scoring, or case-specific tuning.
"""
from __future__ import annotations
import argparse
from copy import deepcopy
from datetime import datetime,timezone
import json
from pathlib import Path

SELF_INDICES=(2,5,7,16,17,24,29,37,38,42,44,45,8,20,32,46)
PAIRS=('00','01','10','11')
BASELINES=('cie','dual_cie','g19','dr','global_dual','direct')


def write_new(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x',encoding='utf8') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')


def common(root):
    if not (root/'INDEPENDENT_CONTRACT.md').is_file():
        raise FileNotFoundError('Write frozen INDEPENDENT_CONTRACT.md before preparation')
    base=json.loads((root/'configs/dualmetric48_v1.json').read_text())
    protocol=deepcopy(base['protocol'])
    if protocol['CPU_cap']!=120 or protocol['proposal_fraction']!=.6:
        raise ValueError('Expected120totalCPU/72proposal protocol')
    protocol.update(lane_policy='All arms per case sequential on the same pinned P core; complete every case/arm despite other failures',
        truth_policy='NONE online; freeze all predictions before independent offline scoring',
        parameter_policy='Frozen complete cohort; no outcome-based exclusions, tuning, or budget extensions',
        corrector_policy='SHARED_ABSOLUTE_TOTAL120_UNUSED_PROPOSAL_REALLOCATED_ALL_ARMS',
        corrector='Original numerical L-BFGS-B; total absoluteCPU120 after<=72proposal; unused proposal reallocated for all arms')
    return base,protocol


def prepare_self_only(root):
    base,protocol=common(root);campaign='self_only16_v1';tasks=[]
    by={(t['case_id'],t['pair']):t for t in base['tasks']}
    for index in SELF_INDICES:
        cid=f'r13_fresh_{index:04d}'
        for pair in ('01','11'):
            t=deepcopy(by[cid,pair]);t.update(campaign=campaign,task_id=cid+'__pair'+pair+'__self_only',onehop_self_only=True,
                cohort='EXISTING_EXPOSED16_FIXED_SELF_ONLY_CONTROL')
            if not (root/t['input_file']).is_file():raise FileNotFoundError(t['input_file'])
            tasks.append(t)
        t=deepcopy(by[cid,'01'])
        t.update(campaign=campaign,task_id=cid+'__global_dual',method='global_dual',algorithm='global_dual',
                 cohort='EXISTING_EXPOSED16_FIXED_SELF_ONLY_CONTROL',
                 current_metric='legacy_scalar_beta_current',material_metric='legacy_scalar_beta_material')
        t.pop('pair',None)
        tasks.append(t)
    protocol.update(scope='EXPOSED16_SELF_ONLY_MECHANISM_NOT_INDEPENDENT',
        source_scope='OriginalexposedR13publicinputs;R14pairs01/11withactualonehopGoffzeroed')
    output=root/'configs'/f'{campaign}.json'
    record=dict(campaign=campaign,status='FROZEN_BEFORE_SELF_ONLY_RUN',
        scope='EXPOSED16_FIXED_SELF_ONLY_NOT_NEW_CONFIRMATION',protocol=protocol,
        fixed_case_indices=list(SELF_INDICES),methods=['pair01','pair11','global_dual'],onehop_self_only='pair arms only; global_dual unchanged',
        frozen_readout_options=deepcopy(base['frozen_readout_options']),tasks=tasks)
    write_new(output,record)
    write_new(root/'audit/SELF_ONLY_PREPARATION.json',dict(status='PREPARED_NOT_EXECUTED',created_utc=datetime.now(timezone.utc).isoformat(),
        manifest=str(output),cases=16,tasks=48,truth_access='NONE',scoring='NOT_RUN',optimization='NOT_RUN'))
    return output


def prepare_independent(root,public_manifest):
    # Parse only the generator's explicit public manifest schema. No guessed
    # generator labels, ranges, input names or material arrays are permitted.
    base,protocol=common(root)
    source=Path(public_manifest).resolve()
    public_dir=(root/'data_public/independent24').resolve()
    if not source.is_relative_to(public_dir):raise ValueError('Require manifest inside declared publicindependent24directory')
    data=json.loads(source.read_text())
    if data.get('schema_version')!='r14_independent24_public_v1' or data.get('material_truth_in_public') is not False:
        raise ValueError('Require the actual declaredpublicgenerator schema withouttruth')
    records=data['case_files']
    if len(records)!=24:raise ValueError('Exactly24publiccase recordsrequired')
    ids=[r['case_id'] for r in records]
    expected=[f'r14_new_{i:04d}' for i in range(24)]
    if len(set(ids))!=24 or sorted(ids)!=expected:raise ValueError('PubliccaseIDs do not match declared24cohort')
    # Require positive public evidence for generator/material-regime provenance;
    # the schema binding below is completed only against the actual generator.
    campaign='independent24_v1';tasks=[]
    for rec in records:
        cid=rec['case_id'];input_file=rec['input_file']
        path=(root/input_file).resolve()
        if not path.is_relative_to(public_dir) or not path.is_file():raise FileNotFoundError(input_file)
        contract=deepcopy(rec)
        for method in tuple('pair'+x for x in PAIRS)+('pair01_to11',)+BASELINES:
            t=dict(task_id=cid+'__'+method,case_id=cid,method=method,algorithm=method,
                campaign=campaign,cohort='INDEPENDENT24_FROZEN_GENERATOR',physics_kind='scalar_tm_2d',
                input_file=input_file,seed=13013,CPU_cap=120,observation_contract=contract)
            if method.startswith('pair'):
                pair='01' if method=='pair01_to11' else method[4:]
                if method=='pair01_to11':t['metric_schedule']=['01','01','01','11']
                t.update(pair=pair,current_metric='M'+pair[0],material_metric='M'+pair[1],
                    sourcewise_cg=False,initialize_physical_current=False,
                    readout_version='R7_FINAL_READOUT' if pair[1]=='1' else 'ORDINARY_SCALAR_B0',
                    readout_options=deepcopy(base['frozen_readout_options']) if pair[1]=='1' else {})
            else:
                t.update(runner_family='HISTORICAL_R13_FROZEN_C2',
                    initialization_scope='Original method-specific full auxiliary initialization retained')
            tasks.append(t)
    protocol.update(scope='INDEPENDENT24_NEW_OBJECTS_FROZEN_BEFORE_PREDICTIONS',
        source_scope='Generator public manifest frozen before predictions; four R14 pairs plus stage candidate plus five C2 numerical baselines with shared-budget extension',
        sourcewise_cg='R14pairsfalse;historicalC2baseline-specificunchanged',
        selected_logical_CPUs=[0,2,4,6,8,10,12,14,16,18,20,22],
        core_type_policy='Hostresourceauditlabels;unspecifiedUNVERIFIED',
        no_absolute_CPU_crosscohort_comparison=True)
    output=root/'configs'/f'{campaign}.json'
    record=dict(campaign=campaign,status='FROZEN_BEFORE_INDEPENDENT_RUN',
        scope='INDEPENDENT24_8_GEOMETRY_FAMILIES_X_3_MATERIAL_REGIMES',protocol=protocol,
        methods=[*('pair'+x for x in PAIRS),'pair01_to11',*BASELINES],generator_public_manifest=str(source.relative_to(root)),
        generator_public_provenance=deepcopy(data),frozen_readout_options=deepcopy(base['frozen_readout_options']),tasks=tasks)
    write_new(output,record)
    write_new(root/'audit/INDEPENDENT_PREPARATION.json',dict(status='PREPARED_NOT_EXECUTED',created_utc=datetime.now(timezone.utc).isoformat(),
        manifest=str(output),generator_public_manifest=str(source),cases=24,tasks=264,truth_access='NONE',scoring='NOT_RUN',optimization='NOT_RUN'))
    policy=root/'configs/public_selection_independent24_v1.json'
    write_new(policy,dict(status='FROZEN_BEFORE_INDEPENDENT_RESULTS',source_rule='Pro scripts/replay_selector_audit.py::public_replay feasibility_gated_saved_uniform',
        original_BASE='a_plus_arc_r7',r14_BASE='pair01',order=['pair01','g19','dr'],
        public_columns=['case_id','solver','Phi','residual_over_delta','noise_feasible','noise_certificate_eligible','external_child_CPU_seconds'],
        stop_exact='if bool(incumbent.noise_feasible): break',replace_exact='if float(trial.Phi) < float(incumbent.Phi): replace incumbent',
        cost_exact='Start with BASE external child CPU; add every invoked child external CPU before comparison, including failed improvements',
        child_initialization='Independent uniform-start stored results, never terminal warmrestart',
        scope='Prospectively frozen saved-result replay benchmark, not measured conditional deployment; no reoptimization',
        source_scope='Only original source public-selection conditions inspected; no archived result thresholds retuned',
        excludes=['NRMSE','IoU','relativeL2','family','truth_arrays','pair01_to11'],
        proposed_CPU_baseline='Sum recorded invoked childCPU; no equal120s totalbudgetclaim; controller/transfer overhead not measured'))
    return output


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--public-manifest',type=Path)
    p.add_argument('--only-self-only',action='store_true')
    p.add_argument('--only-independent',action='store_true')
    a=p.parse_args();root=a.root.resolve()
    if a.only_self_only and a.only_independent:p.error('Choose atmostoneonlyflag')
    if not a.only_self_only and a.public_manifest is None:p.error('--public-manifestrequiredforindependentpreparation')
    outputs=[]
    if not a.only_independent:outputs.append(prepare_self_only(root))
    if not a.only_self_only:outputs.append(prepare_independent(root,a.public_manifest))
    print(json.dumps(dict(status='PREPARED_NOT_EXECUTED',manifests=[str(x) for x in outputs],truth_access='NONE')))


if __name__=='__main__':main()
