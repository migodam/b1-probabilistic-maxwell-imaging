#!/usr/bin/env python3
"""Two auxiliary negative controls, no new FEM, no experiment on import.

Neither control establishes the main new 2D unknown-wall nonrecoverable pair.
"""
import argparse
import ast
import json
from pathlib import Path
import time
import numpy as np
from wslam.accounting import BudgetBook,json_default
from wslam.finite_bayes import PRIOR_ID,update_logweights,conditional_target_only,summarize

ROOT=Path(__file__).resolve().parents[1]
SIGMA=.01


def write_json(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,default=json_default,allow_nan=False)+'\n')
    temporary.replace(path)


def save_npz(path,**arrays):
    temporary=path.with_suffix('.tmp.npz');np.savez_compressed(temporary,**arrays);temporary.replace(path)


def loglikelihood(predictions,observation):
    d=(np.asarray(predictions)-np.asarray(observation))/SIGMA
    axes=tuple(range(1,d.ndim))
    return -.5*np.sum(d.real*d.real+d.imag*d.imag,axis=axes)


def load_audited_functions(path):
    """Load actual definitions unchanged, omitting legacy plotting/file effects."""
    text=path.read_text();tree=ast.parse(text,filename=str(path))
    names=('scattering','full_target_reflection');nodes=[];locations={}
    for name in names:
        found=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==name]
        if len(found)!=1:raise RuntimeError('audited analytic interface missing '+name)
        node=found[0]
        expected=['layers','freq'] if name=='scattering' else ['layers','freq','position','rho']
        if [a.arg for a in node.args.args]!=expected or node.decorator_list:
            raise RuntimeError('audited analytic interface differs: '+name)
        nodes.append(node);locations[name]=dict(start=node.lineno,end=node.end_lineno)
    namespace={'np':np,'__builtins__':__builtins__}
    module=ast.Module(body=nodes,type_ignores=[])
    exec(compile(module,str(path),'exec'),namespace)
    return namespace['full_target_reflection'],locations


def analytic_control(out,source):
    out.mkdir();online=out/'online';offline=out/'offline';online.mkdir();offline.mkdir()
    provenance=dict(source_path=str(source),formula='Audited source scattering and full_target_reflection unchanged',
        scope='AUXILIARY_EXISTING_1D_MAXWELL_SLAB_PEC_SINGLE_FREQUENCY',
        main_new_2D_unknown_wall_control=False,frequency=1.,reference_plane=0.,rho=-1.,
        walls=[[[2.,.25]],[[3.,1/6]]],target_positions=[1.5,1.5-1/12],
        proper_prior=[.5,.5],sigma_each_real_imag=SIGMA,
        parents=12,noise_repeats=2,independent_factors=6,world_seed=2026102509,noise_seed=2026102510)
    try:
        function,locations=load_audited_functions(source)
        old_results=source.parent/'analytic_results/results.json'
        original=json.loads(old_results.read_text()) if old_results.exists() else None
        if original is None or original.get('status')!='executed_analytic_1d_Maxwell_check_only':
            raise RuntimeError('original analytic audit receipt unavailable')
        means=np.array([function([(2.,.25)],[1.],1.5)[0],
                        function([(3.,1/6)],[1.],1.5-1/12)[0]],complex)
        provenance.update(function_source_lines=locations,original_audit_path=str(old_results),
            original_center_difference=original.get('tests',{}).get('joint_center_data_difference_abs'),
            independent_joint_amplitudes=means,center_difference_abs=float(abs(means[0]-means[1])),
            equality_reason='Existing half-wave slab/PEC compensation at fixed reference and center frequency; not equality at other frequencies or a 2D wall result.')
    except Exception as exc:
        provenance.update(status='NOT_RUN',reason='AUDITED_ANALYTIC_INTERFACE_OR_AUDIT_UNAVAILABLE',error=repr(exc))
        write_json(out/'RESULTS.json',provenance);return provenance
    rng=np.random.default_rng(2026102509);truth=rng.choice(2,size=12,p=[.5,.5])
    noise_rng=np.random.default_rng(2026102510)
    noise=SIGMA*(noise_rng.normal(size=(12,2,6))+1j*noise_rng.normal(size=(12,2,6)))
    observations=means[truth,None,None]+noise
    prior=np.array([.5,.5]);posterior=np.empty((12,2,6,2));refresh=np.empty_like(posterior)
    factor_ids=np.empty((12,2,6),dtype='U80')
    for parent in range(12):
        for repeat in range(2):
            factors=[]
            for step in range(6):
                identifier=f'AUX1D-parent-{parent:02d}-repeat-{repeat}-independent-{step}'
                factor_ids[parent,repeat,step]=identifier
                ll=loglikelihood(means[:,None],observations[parent,repeat,step:step+1])
                factors.append((identifier,ll))
                posterior[parent,repeat,step]=np.exp(update_logweights(prior,factors))
                # Reopening refreshes the same six unique factors from the proper prior.
                refresh[parent,repeat,step]=np.exp(update_logweights(prior,list(factors)))
    save_npz(online/'PROBABILITIES_AND_OBSERVATIONS.npz',prior=prior,independently_evaluated_means=means,
        observations=observations,factor_id=factor_ids,posterior_by_step=posterior,reopened_by_step=refresh)
    save_npz(offline/'TRUTHS_AND_NOISE.npz',truth_atom=truth,independent_noise=noise)
    provenance.update(status='AUXILIARY_CONTROL_EXECUTED',maximum_posterior_deviation_from_half=float(np.max(abs(posterior-.5))),
        maximum_reopen_difference=float(np.max(abs(refresh-posterior))),
        interpretation='Distinct wall/PEC-position atoms share center-frequency joint amplitude up to floating point; fresh repetitions and exact refresh do not resolve this inherited analytical equality.',
        claim_limit='Existing analytic 1D two-atom control only; not the new 2D finite material prior or general wall result.')
    write_json(out/'RESULTS.json',provenance)
    return provenance


def known_air_control(out,bank_root):
    out.mkdir();online=out/'online';offline=out/'offline';online.mkdir();offline.mkdir()
    with np.load(bank_root/'online/FINITE_GRID.npz',allow_pickle=False) as archive:
        grid={name:archive[name] for name in archive.files}
    with np.load(bank_root/'online/FINITE_BANK.npz',allow_pickle=False) as archive:
        if int(archive['completed'])!=2950 or str(archive['prior_id'])!=PRIOR_ID:
            raise RuntimeError('complete explicit finite bank required for known-air control')
        bank=np.array(archive['data'])
    indices=np.flatnonzero(~grid['support'].any(axis=1))
    if len(indices)!=25:raise RuntimeError('known-air conditional support must have exactly25 physical poses')
    prior=grid['prior_weight'][indices];prior=prior/prior.sum()
    predictions=bank[indices,:2,:]
    conditional={}
    for name in ('support_id','support','epsinf','ell','pose','wall_group_id'):
        conditional[name]=grid[name][indices]
    conditional.update(prior_weight=prior,legal_supports=grid['legal_supports'],label='AUXILIARY_KNOWN_AIR_CONDITIONAL_25_POSE')
    groups=np.zeros(25,dtype=int)  # The one known constant wall, not a truth-selected freeze.
    rng=np.random.default_rng(2026102511);parent_pose_index=rng.choice(25,size=12,p=prior)
    noise_rng=np.random.default_rng(2026102512)
    noise=SIGMA*(noise_rng.normal(size=(12,2,2,16))+1j*noise_rng.normal(size=(12,2,2,16)))
    observations=predictions[parent_pose_index,None,:,:]+noise
    exact=np.empty((12,2,3,25));stale=np.empty_like(exact)
    ids=np.empty((12,2,3),dtype='U80');summaries=[]
    remaining=np.array([r for r in range(16) if r!=8])
    for parent in range(12):
        for repeat in range(2):
            factors=[];current=prior.copy()
            # Three unique factors collectively exhaust the two base sources:
            # t1 TX0/RX8; t2 TX1/RX8; t3 bothTX/remaining15RX.
            for step in range(3):
                identifier=f'AUXAIR-parent-{parent:02d}-repeat-{repeat}-t{step+1}'
                ids[parent,repeat,step]=identifier
                if step==0:
                    p=predictions[:,0,8:9];y=observations[parent,repeat,0,8:9]
                elif step==1:
                    p=predictions[:,1,8:9];y=observations[parent,repeat,1,8:9]
                else:
                    p=predictions[:,:,remaining];y=observations[parent,repeat][:,remaining]
                ll=loglikelihood(p,y);factors.append((identifier,ll))
                exact[parent,repeat,step]=np.exp(update_logweights(prior,factors))
                current=conditional_target_only(current,ll,groups)['weights']
                stale[parent,repeat,step]=current
                summaries.append(dict(parent=parent,repeat=repeat,step=step+1,
                    summary=summarize(exact[parent,repeat,step],conditional),
                    conditional_vs_batch_max_error=float(np.max(abs(current-exact[parent,repeat,step])))))
    save_npz(online/'CONDITIONAL_GRID.npz',original_bank_index=indices,prior_weight=prior,
        pose=conditional['pose'],known_wall_support=np.zeros((2,2),int),predictions=predictions)
    save_npz(online/'PROBABILITIES_AND_OBSERVATIONS.npz',observations=observations,prior=prior,
        factor_id=ids,exact_by_step=exact,conditional_target_only_by_step=stale,
        known_wall_group=groups,t3_receiver_indices=remaining)
    save_npz(offline/'TRUTHS_AND_NOISE.npz',conditional_pose_index=parent_pose_index,
        original_bank_index=indices[parent_pose_index],true_pose=conditional['pose'][parent_pose_index],independent_noise=noise)
    write_json(online/'SUMMARIES.json',summaries)
    result=dict(status='AUXILIARY_CONTROL_EXECUTED',scope='KNOWN_AIR_WALL_CONDITIONAL_FINITE_POSE_CONTROL',
        proper_conditional_prior=True,pose_states=25,parents=12,noise_repeats=2,
        world_seed=2026102511,noise_seed=2026102512,wall_known_before_data=True,
        no_ground_truth_wall_selection=True,unknown_target_pose=True,
        maximum_conditional_vs_batch_error=float(np.max(abs(exact-stale))),
        interpretation='A single known wall group has no uncertain wall marginal to revisit. Conditional target updates equal exact finite updates; pose remains inferred.',
        channel_rule='t1 TX0/RX8, t2 TX1/RX8, t3 both baseTX/remaining15RX; no repeated channels/factors',
        new_FEM_rhs=0,new_FEM_factorizations=0,
        claim_limit='Known constant air wall only; not the main unknown-wall inference contract or general no-revisit result.')
    write_json(out/'RESULTS.json',result);return result


def run(args):
    out=Path(args.output);out.mkdir(parents=True,exist_ok=False)
    book=BudgetBook(args.ledger);before=book.totals('E2E3')
    report=dict(status='RUNNING',scope='AUXILIARY_CONTROLS_ONLY',main_E3_exact_2D_negative='NOT_ESTABLISHED',
        replaces_main_contract=False,continuous_probability_certificate=False,scientific_acceptance='PARENT_REVIEW_REQUIRED',
        cost_before=before,cpu_accounting='Parent run_job wrapper measures aggregate process CPU; this tool never charges CPU twice',
        statistical_unit='12 parent worlds, two noise repeats per world; not24 independent parents')
    write_json(out/'RESULTS.json',report)
    started=time.perf_counter();cpu=time.process_time()
    try:
        report['analytic_1D']=analytic_control(out/'analytic_1D',Path(args.analytic_source))
        try:report['known_air_2D']=known_air_control(out/'known_air_2D',Path(args.bank))
        except Exception as exc:report['known_air_2D']=dict(status='NOT_RUN',error=repr(exc),main_unknown_wall_contract=False)
        report.update(status='AUXILIARY_CONTROLS_COMPLETE',new_rhs=book.totals('E2E3')['rhs']-before['rhs'],
            new_factorizations=book.totals('E2E3')['factorizations']-before['factorizations'])
        if report['new_rhs'] or report['new_factorizations']:
            raise RuntimeError('unexpected physical charge during algebra-only auxiliary controls')
    except BaseException as exc:
        report.update(status='PARTIAL_FAILURE',error=repr(exc));raise
    finally:
        report.update(cost_after=book.totals('E2E3'),tool_process_cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-started)
        write_json(out/'RESULTS.json',report)
        book.event('E2E3','auxiliary_negative_controls',status=report['status'],new_rhs=report.get('new_rhs',0),
            new_factorizations=report.get('new_factorizations',0),output=str(out),
            diagnostic_tool_cpu_seconds=report['tool_process_cpu_seconds'],main_E3_exact_2D_negative='NOT_ESTABLISHED')


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bank',default=str(ROOT/'results/e2e3/finite_bank'))
    parser.add_argument('--output',default=str(ROOT/'results/e2e3/negative_controls'))
    parser.add_argument('--analytic-source',default=str(ROOT/'legacy_inputs/p0a/analytic_witness.py'))
    parser.add_argument('--ledger',default=None)
    run(parser.parse_args(argv))


if __name__=='__main__':main()
