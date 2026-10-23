"""Synthetic readback of parent repairs and public quadrature/coordinate functions.
No source-module imports, filesystem data inputs, Family construction, main or PDE.
"""
import ast,json,types
from pathlib import Path
import numpy as np
def logsumexp(values):
    values=np.asarray(values,float)
    if not values.size:return -np.inf
    anchor=np.max(values)
    if anchor==-np.inf:return -np.inf
    return float(anchor+np.log(np.exp(values-anchor).sum()))
OUT=Path(__file__).resolve().parent
SRC=OUT.parents[1]/'src/wr'
def isolate(filename,names,extra=None):
    path=SRC/(filename+'.py');tree=ast.parse(path.read_text());selected=[]
    for node in tree.body:
        if isinstance(node,ast.FunctionDef) and node.name in names:selected.append(node)
        if isinstance(node,ast.ClassDef):
            selected.extend(n for n in node.body if isinstance(n,ast.FunctionDef) and n.name in names)
    assert {n.name for n in selected}==set(names)
    env=dict(np=np,logsumexp=logsumexp,CENTERS=((-0.9,2.15),(-0.8,2.15)),U_BOUNDS=(-.01,.01),CONTRAST_BOUNDS=(2.94,3.06),RADIUS=.12,TARGET_CONDUCTIVITY=.02)
    env.update(extra or {});exec(compile(ast.Module(body=selected,type_ignores=[]),str(path),'exec'),env);return env
pi=isolate('progressive_inference',['_particle_summary'])
off=isolate('progressive_offline',['_completion_record','_cost_fields'])
ref=isolate('prob_reference',['posterior'])
bank=isolate('prob_bank',['node_rule'])
family=isolate('probability_family',['_branch','_bounded','coordinates'])
finite=dict(branch_mass=[.5,.5],logZ=0.,CDF={'-0.85':.5},x_mean=-.85,x_sd=.05,qualified=True)
continuous=dict(CDF={'-0.85':{'order9':.5}})
material=dict(x=np.array([-.9,-.8]),branches=np.array([0,1]))
result=dict(particle_indexes=np.array([0,1]),weights=np.array([1.,0.]),log_evidence_estimate=0.,surviving_initial_ancestors=2)
s=pi['_particle_summary'](result,material,dict(n=1),finite,continuous,[-.85])
assert s['branch_MISS'] and s['branch_particle_slot_present']==[True,True] and s['branch_positive_weight_present']==[True,False]
branch={k:s[k] for k in ('branch_mass','branch_MISS','branch_missing','branch_particle_slot_present','branch_positive_weight_present','branch_effective_particle_count')}
result['weights']=np.array([1.,1e-300]);tiny=pi['_particle_summary'](result,material,dict(n=1),finite,continuous,[-.85])
assert tiny['branch_effective_particle_count']==[1.,1.]
branch['tiny_positive_weight_effective_counts']=tiny['branch_effective_particle_count']
records=[{},dict(status='INCOMPLETE_NO_CONVERGENCE_CLAIM',observations=32),dict(status='FULL_NODE_INFERENCE_COMPLETE',inference_frozen=True,truth_access=False,observations=32),dict(status='FULL_NODE_INFERENCE_COMPLETE',inference_frozen=True,truth_access=False,observations=31),dict(status='FULL_NODE_INFERENCE_COMPLETE',inference_frozen=True,truth_access=True,observations=32)]
markers=[]
for record in records:
    try:off['_completion_record'](record,'inference freeze',32);accepted=True
    except (ValueError,RuntimeError):accepted=False
    markers.append(dict(record=record,accepted=accepted))
assert [r['accepted'] for r in markers]==[False,False,True,False,False]
validation=dict(status='FULL_POINT_VALIDATION_COMPLETE',inference_frozen=True,truth_access=False,observations=32)
off['_completion_record'](validation,'point validation',32)
costs=off['_cost_fields'](dict(SMC_runs={'direct':[{'cpu_seconds':3.}]},methods={'point':{'wall_seconds':2.}},shared_proposal_library_cost=[{'cpu_seconds':1.}]))
assert costs['SMC_runs.direct[0].cpu_seconds']==3 and costs['methods.point.wall_seconds']==2
rule=dict(weights=np.array([.5,.5]),x=np.array([-.9,-.8]),contrast=np.array([3.,3.]),branches=np.array([0,1]))
odds=ref['posterior'](rule,np.array([0.,1000.]));assert odds['posterior_log_odds']==1000. and odds['posterior_odds'] is None
cut_checks=[]
for order in (3,5,9):
    for cut,expected in [(None,1.),(-.92,0.),(-.909,.025),(-.9,.25),(-.85,.5),(-.8,.75),(-.791,.975),(-.78,1.)]:
        nodes=bank['node_rule'](order,cut);mass=sum(n['prior_quadrature_mass'] for n in nodes)
        assert abs(mass-expected)<1e-13
        cut_checks.append(dict(order=order,cut=cut,mass=mass,analytic_mass=expected))
self=types.SimpleNamespace()
coordinates=[]
for b,u,c in [(0,-.01,2.94),(1,.01,3.06),(1,np.nextafter(.01,np.inf),3.)]:
    coordinates.append(family['coordinates'](self,b,u,c))
rejects=[]
for b,u,c in [(True,0.,3.),(2,0.,3.),(0,np.nextafter(np.nextafter(.01,np.inf),np.inf),3.),(0,0.,3.061)]:
    try:family['coordinates'](self,b,u,c);rejected=False
    except ValueError:rejected=True
    assert rejected;rejects.append(dict(branch=str(b),u=u,contrast=c,rejected=rejected))
output=dict(parent_repairs=dict(branch_diagnostic=branch,completion_markers=markers,point_validation_complete_accepted=True,nested_costs=costs,odds={k:odds[k] for k in ('branch1','posterior_log_odds','posterior_odds','posterior_odds_status')}),public_truncated_rules=cut_checks,public_coordinates=dict(accepted=coordinates,rejected=rejects),execution=dict(synthetic_arrays_only=True,PDE=0,RHS=0,LU=0,Book_instances=0,Family_instances=0,Family_point_calls=0,material_decoder_calls=0,new_observations_read=0,GT_read=0,main_functions_called=0,SMC_runs=0,source_modules_imported=[]))
(OUT/'REPAIR_READBACK_CHECKS.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
print(json.dumps(output,indent=2))
