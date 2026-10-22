from pathlib import Path
import sys,json,gc,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,Target,write,save,plain
from p3l.conditional import make_frame,conditional_reference,importance_correct,evaluate_slice

def main():
    parents=json.loads((OUT/'A/DEEP_REFERENCE_PARENTS.json').read_text())['parents'][:4]
    protocol=dict(parents=parents,Gaussian_weak_dimensions=[2,4,8,16,24],frame='posterior_material_risk',outer_weak_contexts=['zero','fixed_seed_half_posterior_scale'],properties='native anchor properties, explicitly conditioned; no box Gaussian surrogate',strong_reference='GN mode at most6 full evaluations',importance_sizes=[8,32,128],same_inner_draw_prefix=True,target_family='FULL conditional in b at declared a,u; outer distribution NOT inferred here',GT_used=False,seed=2026102260,primary_question='volume variation, decoder mean-mode gap, exact-weight degeneracy versus strong dimension',new_budget=False)
    path=OUT/'configs/B_CONDITIONAL_PILOT.json'
    if path.exists():raise RuntimeError('conditional pilot already registered; do not replay')
    write(path,protocol)
    b=Book('B','strong_conditional_pilot',max_wall=3600);status='COMPLETE';error=None;rows=[]
    try:
        T=Target(b)
        for ip,parent in enumerate(parents):
            y=T.observation(parent);anchor=T.load_cache(OUT/f'A/anchors/{parent}.npz',owner='P3L A paid optimized anchor')
            J=T.jacobian(anchor)
            for d in protocol['Gaussian_weak_dimensions']:
                frame=make_frame(anchor.q,J,d,kind='posterior_material_risk',material_derivative=T.prior.derivative(anchor.q))
                save(OUT/f'B/conditionals/{parent}_d{d}_frame.npz',**{k:v for k,v in frame.items() if isinstance(v,np.ndarray)})
                rng=np.random.default_rng(np.random.SeedSequence([2026102260,ip,d]))
                for context,a in [('zero',np.zeros(d)),('shift',.5*rng.normal(size=d))]:
                    label=f'{parent}_d{d}_{context}';dest=OUT/'B/conditionals'/label;dest.mkdir(parents=True,exist_ok=True)
                    events=(dest/'EVENTS.jsonl').open('a')
                    def callback(row):events.write(json.dumps(plain(row),allow_nan=False)+'\n');events.flush()
                    sliced=evaluate_slice(T,y,frame,a,anchor.q[64:],callback=callback)
                    cond=conditional_reference(T,y,frame,a,anchor.q[64:],max_evaluations=6,callback=callback)
                    ref=cond['reference'];mom=cond['material_moments'];P=cond['point']
                    save(dest/'CONDITIONAL_REFERENCE.npz',a=a,u=anchor.q[64:],mode_q=P.q,mode_chi=P.state.chi,mean_b=ref['mean_b'],precision_b=ref['precision_b'],covariance_b=ref['covariance_b'],gaussian_material_mean=mom['mean'],gaussian_material_variance=mom['variance'],field_mean=mom['field_mean'],field_sd=mom['field_sd'])
                    state=rng.bit_generator.state;imp=importance_correct(T,y,frame,a,anchor.q[64:],ref,samples=128,rng=rng,callback=callback)
                    save(dest/'IMPORTANCE.npz',strong_samples=imp['strong_samples'],logweights=imp['logweights'],weights=imp['weights'],material_mean=imp['material_mean'],material_variance=imp['material_variance'],endpoint_keys=np.array([p.key for p in imp['points']]))
                    image=np.asarray([p.state.chi for p in imp['points']]);convergence=[]
                    from scipy.special import logsumexp
                    for k in [8,32,128]:
                        lw=imp['logweights'][:k];lz=logsumexp(lw)-np.log(k);ww=np.exp(lw-logsumexp(lw));mean=np.einsum('i,ij->j',ww,image[:k]);convergence.append(dict(K=k,logZ_estimate=float(lz),importance_ESS=float(1/(ww@ww)),max_weight=float(max(ww)),mean_difference_to_K128_RMS=float(np.linalg.norm(mean-imp['material_mean'])/24)))
                    row=dict(parent=parent,weak_dimension=d,strong_dimension=64-d,context=context,target_family='FULL_CONDITIONAL_DIAGNOSTIC_NOT_WEAK_MARGINAL',properties=anchor.q[64:].tolist(),slice_log_density=sliced['log_gamma'],profile_log_density=cond['profile']['log_profile'],GN_logZ=ref['logZ_GN'],importance_logZ=imp['logZ_estimate'],logdet_precision=ref['logdet_precision'],strong_gradient_inf=cond['strong_gradient_inf'],mode_converged=cond['conditional_mode_converged'],GN_material_mean_minus_mode_RMS=float(np.linalg.norm(mom['mean']-P.state.chi)/24),importance_mean_minus_GN_RMS=float(np.linalg.norm(imp['material_mean']-mom['mean'])/24),importance_ESS=imp['importance_ess'],max_weight=imp['maximum_normalized_weight'],importance_prefix_convergence=convergence,probability_status='FINITE_K_IMPORTANCE_NOT_QUALIFIED_REFERENCE',SNIS_moment_bias=True,all_new_physics_paid=True,RNG_before_importance=state)
                    rows.append(row);write(dest/'SUMMARY.json',row);write(OUT/'B/CONDITIONAL_PILOT_RESULTS.json',dict(protocol=protocol,rows=rows));events.close()
                    print(json.dumps(dict(parent=parent,d=d,context=context,ESS=imp['importance_ess'],mean_mode_gap=row['GN_material_mean_minus_mode_RMS'])),flush=True)
                    del imp,cond,image;gc.collect()
    except BaseException as e:status='ERROR';error=repr(e);raise
    finally:b.close(status,error)

if __name__=='__main__':main()
