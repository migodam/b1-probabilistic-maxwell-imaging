"""Synthetic array-only implementation checks; source AST isolation, no main/I/O dependencies."""
import ast
import json
import math
from pathlib import Path
import numpy as np
from scipy.special import logsumexp

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
SRC=ROOT/'src/wr'
def isolate(name,names,extra=None):
    path=SRC/(name+'.py');tree=ast.parse(path.read_text())
    selected=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names]
    assert {n.name for n in selected}==set(names)
    env=dict(np=np,math=math,logsumexp=logsumexp,OBSERVATION_SHAPE=(5,4,12),AMPLITUDE_BOUNDS=(.99,1.01),PHASE_BOUNDS=(-.01,.01),SIGMA=.01,CENTERS=np.array([[-.9,2.15],[-.8,2.15]]),MAX_ARRAY_BYTES=1_200_000_000,AMPLITUDE=(.99,1.01),PHASE=(-.01,.01))
    env.update(extra or {})
    exec(compile(ast.Module(body=selected,type_ignores=[]),str(path),'exec'),env)
    return env

pm=isolate('posterior_math',['_shape','_observations','_blocks','sufficient_statistics','nuisance_quadrature','_logsumexp','log_marginal_likelihood'])
pi=isolate('progressive_inference',['material_grid','gamma_rule','source_statistics','_likelihood_panel','joint_finite_target','profile_shared_gamma','_particle_summary'])
sp=isolate('sequential_particles',['beta_schedule','_at_beta','_neighbors','_log_ratio'])
off=isolate('progressive_offline',['_completion_record'])
hand=isolate('prob_handoff',['dependent_prediction'])
checks={}
material=pi['material_grid'](65,5);gamma=pi['gamma_rule'](21)
prior=material['weights'];mean=float(prior@material['x'])
checks['uniform_GL']=dict(material_mass=float(prior.sum()),branch_mass=[float(prior[material['branches']==b].sum())for b in (0,1)],offset_mean=float(prior@material['u']),contrast_mean=float(prior@material['contrast']),gamma_mass=float(gamma['weights'].sum()),gamma_amplitude_mean=float(gamma['weights']@gamma['amplitude']),gamma_phase_mean=float(gamma['weights']@gamma['phase']),prior_position_mean=mean,prior_position_squared_risk=float(prior@((material['x']-mean)**2)),analytic_prior_risk=.0025+.0001/3)

rng=np.random.default_rng(231023)
F=.05*(rng.normal(size=(2,5,4,12))+1j*rng.normal(size=(2,5,4,12)))
y=.05*(rng.normal(size=(5,4,12))+1j*rng.normal(size=(5,4,12)))
stats=pi['source_statistics'](F,y);g=pi['gamma_rule'](3);panel=pi['_likelihood_panel'](stats,g,0,2)
constant=-120*np.log(.01*np.sqrt(2*np.pi));direct=np.empty(panel.shape)
for m in range(2):
 for k,(a,p) in enumerate(zip(g['amplitude'],g['phase'])):
  for s in range(4):direct[m,k,s]=constant-.5*np.sum(abs(a*np.exp(1j*p)*F[m,:,s,:]-y[:,s,:])**2)/.01**2
joint=logsumexp(panel.sum(2)+np.log(g['weights']),axis=1)
wrong=sum(logsumexp(panel[:,:,s]+np.log(g['weights']),axis=1)for s in range(4))
whole=pm['log_marginal_likelihood'](F,y,n=3);blocks=pm['log_marginal_likelihood'](F,y,n=3,source_blocks=True)
checks['shared_gamma']=dict(max_panel_vs_direct_error=float(np.max(abs(panel-direct))),max_algebra_module_vs_progressive_error=float(np.max(abs(np.asarray(whole['logZ'],float)-joint))),max_joint_vs_block_statistics_error=float(np.max(abs(whole['logZ']-blocks['logZ']))),independent_per_source_integral_changes_target_max=float(np.max(abs(wrong-joint))))

shape=(2,3,1,2,2);n=math.prod(shape);prior=np.arange(1,n+1,dtype=float);prior/=prior.sum();llblocks=rng.normal(size=(n,4))*.3
proposal=np.zeros((n,n))
for i in range(n):
 for axis in range(5):
  for sign in (-1,1):
   j=sp['_neighbors'](np.array([i]),np.array([axis]),np.array([sign]),shape)[0];proposal[i,j]+=.1
balance=rowsum=increment_error=0.
for name in ('direct','progressive'):
 beta=sp['beta_schedule'](name);d=np.diff(np.vstack((np.zeros(4),beta)),axis=0)
 increment_error=max(increment_error,float(np.max(abs(sum(sp['_at_beta'](llblocks,x)for x in d)-llblocks.sum(1)))))
 for stage in beta:
  ll=sp['_at_beta'](llblocks,stage);target=prior*np.exp(ll);target/=target.sum();K=np.zeros((n,n))
  for i in range(n):
   for j in range(n):
    ag=min(1.,math.exp(ll[j]-ll[i]));al=min(1.,math.exp(ll[j]-ll[i])*prior[j]/prior[i]);move=.25*prior[j]*ag+.75*proposal[i,j]*al
    K[i,j]+=move;K[i,i]+=.25*prior[j]*(1-ag)+.75*proposal[i,j]*(1-al)
  balance=max(balance,float(np.max(abs(target[:,None]*K-(target[:,None]*K).T))));rowsum=max(rowsum,float(np.max(abs(K.sum(1)-1))))
checks['SMC_kernel']=dict(nodes=n,stages_checked=8,local_proposal_symmetry_error=float(np.max(abs(proposal-proposal.T))),maximum_detailed_balance_defect=balance,transition_row_sum_error=rowsum,increment_sum_vs_full_error=increment_error,actual_SMC_runs=0)

result=dict(particle_indexes=np.array([0,1]),weights=np.array([1.,0.]),log_evidence_estimate=0.,surviving_initial_ancestors=2)
finite=dict(branch_mass=[.5,.5],logZ=0.,CDF={'-0.85':.5},x_mean=-.85,x_sd=.05,qualified=True)
continuous=dict(CDF={'-0.85':{'order9':.5}})
summary=pi['_particle_summary'](result,dict(x=np.array([-.9,-.8]),branches=np.array([0,1])),dict(n=1),finite,continuous,[-.85])
checks['zero_weight_branch_counterexample']=dict(actual_branch_mass=summary['branch_mass'],eligible=summary['branch_miss_eligible'],reported_branch_missing=summary['branch_missing'],reported_branch_MISS=summary['branch_MISS'],expected_positive_weight_branch_missing=[False,True],unweighted_occupancy_semantics='different from effective posterior support')
markers=[]
for marker in ({},{'status':'INCOMPLETE_NO_CONVERGENCE_CLAIM','observations':32}):
 try:off['_completion_record'](marker,'synthetic',32);accepted=True
 except Exception:accepted=False
 markers.append(dict(marker=marker,accepted_as_frozen=accepted))
checks['completion_marker_counterexamples']=markers

w=np.array([[.4,.1],[.1,.4]]);response=np.zeros((2,5,4,12),complex);response[0]=1;response[1]=2j
prediction=hand['dependent_prediction'](w,response,np.array([.99,1.01]),np.array([0.]))
gamma2=np.array([.99,1.01]);explicit=sum(w[m,k]*gamma2[k]*response[m]for m in range(2)for k in range(2))
checks['handoff_prediction']=dict(max_joint_vs_explicit_error=float(np.max(abs(prediction['posterior_predictive_approx']-explicit))),factorized_difference_norm=prediction['factorized_predictive_difference_norm'])
checks['execution']=dict(synthetic_arrays_only=True,PDE=0,RHS=0,LU=0,Book_instances=0,new_observations_read=0,GT_read=0,main_functions_called=0,SMC_runs=0,source_modules_imported=[])
(OUT/'ARRAY_CHECKS.json').write_text(json.dumps(checks,indent=2,allow_nan=False))
print(json.dumps(checks,indent=2))
