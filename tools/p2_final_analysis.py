"""Actual final comparisons, prior refinement, and complete resource receipts."""
import os
from pathlib import Path
import sys,json,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.b1_joint.accounting import Ledger,dump
from src.b1_joint.reference import Reference
from src.b1_joint.response import Response
from src.b1_joint.package import update_context
def main():
    out=ROOT/'research/joint_p2/results';ledger=Ledger(out);C1=Path(os.environ.get('C1_P2',ROOT.parent.parent/'C1/c1-joint-p2'));methods=[]
    with ledger.job('p2-final-comparison-schema-resume'):
        sensitivity=json.loads((out/'PRIOR_NUISANCE_SENSITIVITY.json').read_text());ref=Reference(Response(out/'COMMON_RESPONSE_BANK.npz'));items=json.loads((out/'DEVELOPMENT_CASES.json').read_text())['cases']
        for item in ([] if (out/'PRIOR_NUISANCE_SENSITIVITY_REFINED.json').exists() else items):
            r=ref.solve(np.array(item['y'])[:32],np.array(item['sigma'])[:32],order=17);a=ref.adaptive_physical(r,prior_measure='radius_uniform_with_atoms')
            if not a['converged']:raise RuntimeError('separate radius-prior sensitivity unresolved')
            old=sensitivity['cases'][item['id']]['radius_uniform_gain015'];old['old_coefficient_node_reweighting']=old['posterior_event_mass'];old['posterior_event_mass']=a['event_mass'];old['independent_coordinate_adaptive']=a;old['numerical_status']='independent smooth-coordinate adaptive integral for distinct radius-uniform prior';old['refinement_difference']=float(abs(np.array(old['old_coefficient_node_reweighting'])-a['event_mass']).max())
        if not (out/'PRIOR_NUISANCE_SENSITIVITY_REFINED.json').exists():dump(out/'PRIOR_NUISANCE_SENSITIVITY_REFINED.json',sensitivity)
        summary=json.loads((C1/'p2/results/summary/METHOD_SUMMARY.json').read_text());train=json.loads((C1/'p2/results/models/TRAINING_SUMMARY.json').read_text());corr=json.loads((out/'CORRECTED_PROPOSALS.json').read_text());smc=json.loads((out/'GENERIC_SMC.json').read_text());p=np.array(json.loads((out/'FULL_MIE_CORRECTED_REFERENCE.json').read_text())['cases']['large_cavity']['event_mass'])
        for row in summary['seed_separated']:
            methods.append(dict(method=row['method']+':'+str(row['seed']),event_mass_MAE=row['parent_equal_weight_means']['event_mass_MAE'],wall_seconds=None,target='fresh common proxy target',metrics=row['parent_equal_weight_means'],scope='raw conditional density; training plus paidteacher extra; cannot compare timing with fullMie SMC as if equal accuracy'))
        for name,a in corr['methods'].items():methods.append(dict(method=name+' full-IS64',event_mass_MAE=a['comparison_to_full_Mie_reference']['corrected_event_MAE'],wall_seconds=a['physical_cost']['full_forward_wall_seconds'],full_calls=a['physical_cost']['new_material_frequency_evaluations'],scope='preregistered development large_cavity; finite importance sampling',target='full-Mie numerical C target',ESS=a['ESS']))
        for a in smc['runs']:methods.append(dict(method='generic SMC64:'+str(a['seed']),event_mass_MAE=float(abs(np.array(a['event_mass'])-p).mean()),wall_seconds=None,full_calls=a['actual_full_evaluations'],target='full-Mie numerical C target',scope='single development observation, finite population seed variation; no per-seed wall receipt, combined job wall in ledger; not model-wide superiority trial'))
        dump(out/'METHOD_COMPARISON.json',dict(methods=methods,C1_fresh_summary=summary,training_receipt=train,matched_corrected_proposals=corr,SMC=smc,posterior_target_separation=True))
        sample=np.load(out/'packages/solid_samples.npz');r=update_context(sample['latent'],sample['weights'],lambda x:np.zeros(len(x)));dump(out/'HANDOFF_IDENTITY.json',dict(no_context_max_weight_change=float(abs(r['weights']-sample['weights']).max()),whole_modeled_mass=float(r['weights'].sum()),event_partition_sum=float(r['event_mass'].sum()),outside_model_support_mass=None,physics=0))
    dump(out/'RESOURCE_SUMMARY.json',ledger.summary())
if __name__=='__main__':main()
