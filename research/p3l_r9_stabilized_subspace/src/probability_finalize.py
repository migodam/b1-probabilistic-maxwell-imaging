from probability_transfer import *
from probability_integrate import global_audit
import shutil

def run():
    if (RES/'audit_detail_replay_receipt.json').exists():raise RuntimeError('already executed')
    b=Book('audit_detail_replay');status='FAILED';error=None
    try:
        ph=Physics(b);ph.load_setup();integration=json.loads((RES/'integration_results.json').read_text());final={}
        for name,rr in integration.items():
            # Original sparse audit retained; detailed replay explicitly charged.
            shutil.copyfile(RES/f'{name}_whole_support_audit.json',RES/f'{name}_whole_support_audit_first_run.json')
            mode=json.loads((RES/f'{name}_fresh_mode.json').read_text());mode['center']=np.array(mode['center']);mode['C']=np.array(mode['C'])
            audit=global_audit(ph,name,mode)
            sources=rr['source_coordinate_integrals'];ref=rr['reference'];older=json.loads((RES/f'{name}_reduced_GH21.json').read_text())
            means=[np.asarray(x['mean_increment']) for x in [sources[-1],sources[-2],ref,older]]
            meanerr=max(np.max(abs(means[0]-means[1])),np.max(abs(means[0]-means[2])),np.max(abs(means[2]-means[3])))
            keys=['source_layer_mass','latent_phase_gt_half','material_real_increment_gt_002','material_nonhost_event','ROI256_hard_event']
            eventerr=max(float(np.max(abs(np.asarray(a[k])-np.asarray(bb[k])))) for k in keys for a,bb in [(sources[-1],sources[-2]),(sources[-1],ref),(ref,older)])
            # Prior-integrated excluded union has mass <=1, not one per leaf.
            Znumeric=min(sources[-1]['scaled_evidence'],ref['scaled_evidence'])*.999
            omitted=audit['numerical_scaled_omitted_evidence_upper']/Znumeric
            meanerr+=RHI*omitted;eventerr+=omitted
            q=dict(status='NUMERICAL_REFERENCE_QUALIFIED' if meanerr<=2e-6 and meanerr/ref['material_ROI_SD']<=.002 and eventerr<=.0004 else 'REFERENCE_UNQUALIFIED',mean_absolute_numeric_uncertainty_proxy=meanerr,mean_uncertainty_over_SD=meanerr/ref['material_ROI_SD'],event_or_layer_numeric_uncertainty_proxy=eventerr,scaled_evidence_numerical_lower_estimate=Znumeric,numerical_omitted_posterior_upper=omitted,relative_evidence_difference_source_fine_vs_GH31=abs(sources[-1]['scaled_evidence']/ref['scaled_evidence']-1),relative_evidence_source_refinement_discrepancy=abs(sources[-1]['scaled_evidence']/sources[-2]['scaled_evidence']-1),strict_reference_interval=None,qualification='numerical convergence and whole-support float audit; uncertainty proxies are not mathematically rigorous error bounds',first_run_qualification_correction='v2 takes maximum event discrepancy across BOTH source refinement and cross-coordinate pairs; original qualification reported only final cross-coordinate event difference')
            comparisons=[]
            for method in rr['baselines']:
                merr=float(np.max(abs(np.array(method['mean_increment'])-means[2])));eerr=max(float(np.max(abs(np.asarray(method[k])-np.asarray(ref[k])))) for k in keys)
                comparisons.append(dict(method=method['method'],physics=method['physics'],order=method['order'],numerical_mean_error=merr,numerical_mean_error_over_SD=merr/ref['material_ROI_SD'],reference_uncertainty_proxy=meanerr,mean_abs_gate_with_proxy=merr+meanerr<=1e-5,mean_SD_gate_with_proxy=(merr+meanerr)/ref['material_ROI_SD']<=.01,event_or_layer_error=eerr,event_gate_with_proxy=eerr+eventerr<=.002,wall=method['wall_seconds'],cost_counts_delta=method['cost_counts_delta'],task_status='NUMERICAL_TASK_PASS' if merr+meanerr<=1e-5 and (merr+meanerr)/ref['material_ROI_SD']<=.01 and eerr+eventerr<=.002 else 'FAIL'))
            byorder={}
            setup=json.loads((RES/'probe_receipt.json').read_text());setupwall=sum(setup['timings'][k] for k in ['setup_full_host_wall','Schur256_and_transforms_wall','setup_diagnostics_wall'])
            for order in (9,11):
                dr=next(x for x in comparisons if x['physics']=='direct' and x['order']==order);red=next(x for x in comparisons if x['physics']=='reduced' and x['order']==order)
                byorder[order]=dict(direct_integral_wall=dr['wall'],shared_integral_wall=red['wall'],shared_full_setup_wall=setupwall,shared_cold_numerical_path_wall=red['wall']+setupwall,direct_to_shared_integral_wall_ratio=dr['wall']/red['wall'],direct_to_shared_setup_included_wall_ratio=dr['wall']/(red['wall']+setupwall),interpretation='physical reuse factor only; same probability rule, same fresh reduced mode geometry shared. Excludes common mode/reference/validation costs, listed separately. Single-run timing, not statistical speed estimate.')
            # Historical failure replay is algebra on saved old summary, no new Maxwell calls.
            hist=json.loads((R6/'track_a/RESULTS.json').read_text());hrows=[]
            for row in hist['rows']:
                if row['observation']==name and row['method'].startswith('S3'):
                    er=np.linalg.norm(np.asarray(row['mean_increment'])-means[2])/ref['material_ROI_SD']
                    hrows.append(dict(method=row['method'],prior_mass_quadrature=row['prior_mass_quadrature'],mean_error_vector_norm_SD=float(er),old_source_layer_mass=row['source_layer_mass'],qualification='historical cached output replay; no new physical discovery, not timing comparator'))
            final[name]=dict(reference=ref,reference_qualification=q,numerical_old_gate_comparisons=comparisons,two_factor_physical_reuse=byorder,historical_global_grid_failure_replay=hrows,flat_material_prior_wrong_shift_SD=rr['flat_prior_wrong_mean_shift_SD'],whole_support_audit=audit)
            write(RES/f'{name}_reference_qualification_v2.json',q)
        write(RES/'final_numerical_comparison.json',final)
        status='COMPLETE'
    except BaseException as e:error=repr(e);write(AUD/'finalize_failure.json',dict(error=error));raise
    finally:b.close(status,error)
if __name__=='__main__':run()
