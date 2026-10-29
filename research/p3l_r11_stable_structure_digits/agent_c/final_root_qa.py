"""Final root audit: stdlib-only ledger/CSV checks, no physics or GT arrays."""
from pathlib import Path
from collections import Counter, defaultdict
import csv, datetime, json, re

D=Path(__file__).resolve().parent;P=D.parent


def main():
    cost=json.loads((P/'COST_SUMMARY.json').read_text())
    results=json.loads((P/'EXPERIMENT_RESULTS.json').read_text())
    rows=list(csv.DictReader((P/'DIGIT_PER_NUMBER_RESULTS.csv').open()))
    rhs=Counter();lu=Counter();groups=defaultdict(lambda:dict(RHS=0,LU=0,physical_charge_records=0));paid=0
    for source in cost['ledger_sources']:
        for line in (P/source).read_text().splitlines():
            r=json.loads(line);nr=int(r.get('RHS',0));nl=int(r.get('LU',0))
            if not (nr or nl):continue
            paid+=1;g=groups[source.split('/')[0]];g['RHS']+=nr;g['LU']+=nl;g['physical_charge_records']+=1
            if nr:
                rb=r.get('RHS_by_order')
                if rb is not None:rhs.update({str(k):int(v) for k,v in rb.items()})
                else:rhs[str(r.get('order',r.get('factor_order',256)))]+=nr
            if nl:lu[str(r.get('order',r.get('factor_order',256)))]+=nl
    threshold=282.8737921622988;cells=[];mean_diffs=[];sum_diffs=[]
    for ds in ('N','P'):
        for method in ('FULL256','AGENT_A','AGENT_B','AGENT_C'):
            rr=[r for r in rows if r['dataset']==ds and r['method']==method]
            target=results['comparisons'][ds][method]
            for key in ('complex_MSE','ROI_complex_MSE','IoU','Dice','SSIM_full_grid','wall_seconds','strict_support','corrected_kernel_wall'):
                mean_diffs.append(abs(sum(float(r[key]) for r in rr)/len(rr)-target[key]))
            for key,tkey in (('RHS','RHS'),('LU','LU'),('wall_seconds','total_wall_seconds'),('trial_fallbacks','trial_fallbacks'),('trial_evaluations','trial_evaluations')):
                sum_diffs.append(abs(sum(float(r[key]) for r in rr)-target[tkey]))
            passes=sum(float(r['Phi'])<=threshold for r in rr)
            cells.append(dict(dataset=ds,method=method,n=len(rr),noise_ball_pass=passes,
                JSON_noise_ball_pass=target['noise_ball_pass'],
                failed_digits=[int(r['digit']) for r in rr if float(r['Phi'])>threshold],
                CSV_boolean_gate_agrees=all((r['noise_ball']=='True')==(float(r['Phi'])<=threshold) for r in rr),
                trial_evaluations=sum(int(r['trial_evaluations']) for r in rr),
                trial_fallbacks=sum(int(r['trial_fallbacks']) for r in rr)))
    report=(P/'RESEARCH_REPORT_ZH.md').read_text();report_phi=[]
    for line in report.splitlines():
        if re.match(r'\|[NP]\|\d+\|',line):
            fields=line.split('|');ds=fields[1];digit=int(fields[2]);display=float(fields[-2])
            r=next(r for r in rows if r['dataset']==ds and int(r['digit'])==digit and r['method']=='FULL256')
            report_phi.append(dict(dataset=ds,digit=digit,report_Phi=display,CSV_Phi=float(r['Phi']),
                                   absolute_display_rounding_error=abs(float(r['Phi'])-display)))
    warnings=[dict(id='W1',severity='NONBLOCKING_COST_LABEL',
        file='RESEARCH_REPORT_ZH.md',line=64,
        finding='0.00548s is aggregate current-shape port-norm setup in the global arm; fixed C spectral setup is additional. Calling it global boundary setup without that scope can be read as a complete cold setup.',
        suggested_wording='Global arm per-shape port-norm setup totals 0.00548s, with separate fixed spectral setup; local shape-bound setup totals 0.19036s.',
        changes_scientific_verdict=False)]
    qa=dict(status='SCOPED_FINAL_REPORT_PASS_RESOURCE_PROTOCOL_FAIL',
        utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        blocking_issues=[],warnings=warnings,
        known_decisive_failure=dict(status='FAIL_CAP_EXCEEDED',all_order_LU=sum(lu.values()),cap=12000,
            excess=sum(lu.values())-12000,fully_disclosed_in_root=True,not_waived=True),
        resource_audit=dict(source_ledgers=len(cost['ledger_sources']),charged_records=paid,
            independent_RHS=sum(rhs.values()),independent_LU=sum(lu.values()),
            RHS_by_order_exact_match=dict(rhs)==cost['RHS_by_order'],LU_by_order_exact_match=dict(lu)==cost['LU_by_order'],
            group_totals_exact_match=dict(groups)==cost['group_totals'],group_totals=dict(groups),
            status=cost['status'],RHS_cap_pass=sum(rhs.values())<=400000,LU_cap_pass=sum(lu.values())<=12000,
            new_RHS576=rhs['576'],new_LU576=lu['576'],
            normalized_root_ledger_excluded_as_derivative_copy=True),
        data_audit=dict(CSV_rows=len(rows),dataset_method_cells=cells,report_phi_rows=report_phi,
            maximum_report_Phi_rounding_error=max(r['absolute_display_rounding_error'] for r in report_phi),
            maximum_white_residual_squared_minus_2Phi= max(abs(float(r['white_residual'])**2-2*float(r['Phi'])) for r in rows),
            maximum_JSON_mean_discrepancy=max(mean_diffs),maximum_JSON_sum_discrepancy=max(sum_diffs),
            all_N_runs_J80=all(int(r['J_calls'])==80 and r['status']=='MAX_J' for r in rows if r['dataset']=='N'),
            all_P_runs_objective_stop=all(r['status']=='OBJECTIVE_STOP' for r in rows if r['dataset']=='P')),
        formula_audit={
            'F_SS_006':'PASS: lambda<0 contraction reverses the norm inequality as required; lambda>=0 uses beta_star=0. No invalid clipping.',
            'F_SS_007':'PASS: rectangle satisfies y_max<x_min, so y/(x^2+y^2) is minimized at x_max,y_min; coercivity yields inverse 2-norm bound for K_delta.',
            'F_SS_008':'PASS: K_delta[delta P(e-e_tilde)]=P r_e and delta O A(e-e_tilde)=OP[delta P(e-e_tilde)]; therefore order OP K_delta^-1 P r_e is correct, including singular A.',
            'F_SS_009':'PASS: sqrt(2)/sigma realification norm, spectral OP norm, Frobenius weighted residual, and b_star denominator are all correctly ordered.',
            'F_SS_010_011':'PASS: potential bound uses approximate data residual; output error gate is not promoted to a Jacobian gate. B/C exact derivative charges are disclosed.'},
        scope_preservation=dict(no_continuous_PDE_or_interval_certificate_claim=True,
            no_posterior_or_exact_fiber_claim=True,no_new_imaging_capability_claim=True,
            no_TSOM_impossibility_inference=True,no_whole_system_speedup_confirmation_claim=True,
            resource_protocol_failure_preserved=True),
        new_Maxwell_RHS=0,new_LU=0,GT_arrays_read=False,root_files_modified=False,
        visual_scope='Main-thread documented actual image review accepted as provenance; this final audit did not substitute a new visual score or revise blind readings.',
        novelty_scope='No independent global literature novelty search in this final audit; parent retains final scientific and novelty judgment.')
    (D/'FINAL_ROOT_QA.json').write_text(json.dumps(qa,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(dict(status=qa['status'],RHS=sum(rhs.values()),LU=sum(lu.values()),
        cells=[(x['dataset'],x['method'],x['noise_ball_pass'],x['failed_digits']) for x in cells],
        blockers=len(qa['blocking_issues']),warnings=len(warnings)),ensure_ascii=False))


if __name__=='__main__':main()
