from pathlib import Path
import json
root=Path(__file__).resolve().parents[1]
obj=lambda props, req=None: {'type':'object','properties':props,'required':req or list(props),'additionalProperties':False}
string={'type':'string','minLength':1}
number={'type':'number'}
nonnegative={'type':'number','minimum':0}
prob={'type':'number','minimum':0,'maximum':1}
nullable=lambda s:{'anyOf':[s,{'type':'null'}]}
array=lambda s:{'type':'array','items':s}
ref=obj({'uri':string,'format':string,'shape':array({'type':'integer','minimum':0}),'dtype':string,'units':string,'provenance_id':string})
uncertainty=obj({'status':{'enum':['certified_bound','finite_empirical_audit','mc_only','unquantified']},'value':nullable(nonnegative),'metric':string,'scope':string,'confidence':nullable(prob),'source':string})
weight=obj({'kind':{'enum':['none','design_prior','proposal','importance_corrected','independent_reference']},'normalized_weights':nullable(ref),'log_prior':nullable(ref),'log_proposal':nullable(ref),'log_likelihood':nullable(ref),'base_measure':string,'data_history_id':string,'target_id':string,'normalization_status':string,'error':uncertainty,'coverage_status':string})
noise=obj({'packing':string,'covariance':ref,'cross_history_covariance':nullable(ref),'real_channel_variance_convention':string,'dwell_scaling':string,'reference_id':string})
cost=obj({'measurement_seconds':nullable(nonnegative),'switching_seconds':nullable(nonnegative),'energy_joules':nullable(nonnegative),'financial_cost':nullable(nonnegative),'currency':nullable(string),'compute_seconds':nullable(nonnegative),'maxwell_attempts':{'type':'integer','minimum':0},'factorizations':{'type':'integer','minimum':0},'rhs':{'type':'integer','minimum':0},'setup_amortization_count':{'type':'integer','minimum':1},'cost_status':{'enum':['measured','estimated','not_available']}})
action=obj({'id':string,'kind':{'enum':['frequency','illumination','polarization','receiver','near_field','reference_calibration','source_combination','external_context','repeat','stop','bundle']},'frequency_kR':nullable(number),'frequency_hz':nullable(nonnegative),'source_geometry':nullable(ref),'receiver_geometry':nullable(ref),'polarization_or_source_weights':nullable(ref),'constituent_ids':array(string),'feasibility':{'enum':['verified_hardware','simulation_only','unverified','infeasible']},'feasibility_evidence':string,'noise':noise,'nuisance_group_ids':array(string),'new_independent_noise':{'type':'boolean'},'cost':cost})
nuisance=obj({'parameter_names':array(string),'prior_definition':string,'sharing_graph':ref,'posterior_representation':string,'joint_particles_or_conditional_model':ref,'drift_transition':nullable(string),'bounds':nullable(ref)})
rom=obj({'used':{'type':'boolean'},'model_id':nullable(string),'valid_material_domain':string,'valid_action_ids':array(string),'response_error':uncertainty,'refresh_rule':string,'setup_and_refresh_cost':cost})
candidate=obj({'id':string,'material':ref,'structure_description':string,'event_ids':array({'type':'integer','minimum':0}),'base_solution':{'type':'boolean'},'full_maxwell_predictions':array(ref),'prediction_history_id':string,'material_distance_to_base':nonnegative,'whitened_distance_to_data':nonnegative,'distance_definitions':string,'compensation_constraints':string,'physical_feasibility':string,'rom':rom})
event=obj({'id':{'type':'integer','minimum':0},'definition':string,'mass':nullable(prob),'mass_error':uncertainty,'conditional_particles':nullable(ref),'conditional_weights':nullable(ref),'representative_candidate_ids':array(string),'is_separated_mode_proven':{'type':'boolean'}})
package=obj({'schema_version':{'const':'research-b-mhp-1.0'},'package_id':string,'created_at':string,'producer':string,'source_revision':string,'target_id':string,'data_history_id':string,'model':obj({'physical_equations':string,'solver_id':string,'discretization_id':string,'material_domain':string,'material_metric':string,'dispersion_model':string,'time_convention':string,'reference_frequency':number}),'data':obj({'observations':ref,'acquisition_history':array(action),'noise':noise}),'nuisance':nuisance,'candidates':array(candidate),'weights':weight,'partition':obj({'version':string,'events':array(event),'partition_covers_modeled_support':{'type':'boolean'},'displayed_event_ids':array({'type':'integer','minimum':0}),'unshown_modeled_mass':nullable(prob),'unshown_modeled_representation':nullable(ref),'outside_model_mass':nullable(prob),'outside_model_mass_status':string,'unexplored_modes':array(string)}),'available_actions':array(action),'resource_ledger':ref,'limitations':array(string)})
package.update({'$schema':'https://json-schema.org/draft/2020-12/schema','$id':'https://example.invalid/research-b/MaterialHypothesisPackage.schema.json','title':'MaterialHypothesisPackage'})
interval=obj({'estimate':nullable(number),'lower':nullable(number),'upper':nullable(number),'confidence':nullable(prob),'error_scope':string})
recommendation=obj({'schema_version':{'const':'research-b-acquisition-1.0'},'input_package_id':string,'input_data_history_id':string,'status':{'enum':['RECOMMEND_MEASURE','STOP_TASK_SUFFICIENT','STOP_NOT_WORTH_COST','PRESERVE_AMBIGUITY','ABSTAIN_POSTERIOR_UNRELIABLE','EXPAND_SUPPORT_FIRST','ABSTAIN_UTILITY_RANK_UNRESOLVED','BUDGET_EXHAUSTED','SIMULATION_DESIGN_ONLY']},'selected_action_ids':array(string),'criterion':string,'criterion_probability_status':string,'planning_horizon':{'type':'integer','minimum':0},'loss_definition':string,'current_risk':interval,'expected_future_risk':interval,'risk_reduction':interval,'event_eig':interval,'posterior_event_quantiles':nullable(ref),'posterior_quantile_definition':string,'resolved_ambiguities':array(string),'remaining_ambiguities':array(string),'calibration_risks':array(string),'missing_support_status':string,'cost':cost,'utility_error_components':array(uncertainty),'action_ranking_scope':string,'stop_reason':nullable(string),'update':obj({'joint_likelihood_rule':string,'persistent_nuisance_rule':string,'candidate_refresh_rule':string,'normalization_and_remainder_rule':string,'reuse_or_new_data_ids':array(string),'failure_action':string}),'limitations':array(string)})
recommendation.update({'$schema':'https://json-schema.org/draft/2020-12/schema','$id':'https://example.invalid/research-b/AcquisitionRecommendation.schema.json','title':'AcquisitionRecommendation'})
for name,data in [('MaterialHypothesisPackage',package),('AcquisitionRecommendation',recommendation)]:
    (root/'schemas'/f'{name}.schema.json').write_text(json.dumps(data,indent=2),encoding='utf-8')

protocol={
 'protocol_id':'RESEARCH-B-DESIGN-20261014-v1','status':'SPECIFICATION_NOT_EXECUTED',
 'E0':{'maxwell_attempt_cap':0,'training_enabled':False,'allowed_kR':[1.6,2.4],
       'rules':['Only existing responses; cache misses fail closed.',
                'Existing observed channels are not new independent measurements.',
                'A retrospective subset test must recompute its conditioning posterior.',
                'Every new-noise simulation is labeled simulated, not observed hardware data.']},
 'E1':{'execution_authorized':False,'proposed_attempt_cap':4096,
       'proposed_kR_menu':[1.2,1.6,2.0,2.4,2.8,3.2],
       'already_cached_kR':[1.6,2.4],'new_kR':[1.2,2.0,2.8,3.2],
       'allocation':{'new_tensor_bank_4x625':2500,'independent_field_audit_64x4':256,
                     'new_truth_material_frequency_32x4':128,'selected_full_field_checks':512,
                     'failure_and_reference_refinement_reserve':700},
       'fresh_parents':32,'noise_replicates_per_parent':2,
       'shared_gain':'one persistent proper truncated real gain per episode; never per frequency',
       'noise_convention':'P2 real variance sigma^2, complex total variance 2*sigma^2',
       'frequency_dispersion':'constant epsilon radial calibration benchmark only; separate causal-model stress target',
       'polarization':'existing two transverse channels; span/redundancy control only on radial isotropic support',
       'nonconcentric_vector_branch':'separate approval and budget; not auto-started',
       'historical_joint_p2_balance_reused':False},
 'E2':{'execution_authorized':False,'max_new_measurement_blocks':3,'lookahead_depth':2,
       'require_independent_budget':True,'keep_unresolved_worlds':True},
 'E3':{'execution_authorized':False,'arms':['B1_reference','A_candidates','A_plus_limited_full_correction','TSOM_base_only'],
       'primary_endpoint':'true full material MAE at matched total cost',
       'require_known_proposal_for_importance_weights':True,'require_independent_budget':True},
 'policies':['stop_no_new_data','preset_kR_2.4','random','uniform_sweep','raw_candidate_separation',
             'classical_Fisher_D','joint_nuisance_projected_Fisher','event_EIG','expected_material_risk',
             'compensation_aware_robust'],
 'statistics':{'unit':'parent object','paired_bootstrap':2000,'primary_material_metric':'volume-weighted componentwise complex-permittivity MAE',
               'mae_normalization':'fixed prespecified physical scale, never truth-dependent',
               'model_mismatch_panels':'reported separately, not pooled with in-model performance'},
 'go_no_go_proposals':{'E1':'same cost >=10% relative MAE improvement OR >=20% acquisition-cost saving at matched MAE, paired interval excludes zero; no meaningful coverage deterioration',
                      'ROM':'same accepted utility/decision quality, >=20% inclusive end-to-end savings including refresh, or disable ROM',
                      'graph':'retain only if measurable quality/cost benefit over direct risk computation',
                      'notes':'Engineering relevance thresholds proposed before experiments, not physical laws or already passed gates.'}
}
assert sum(protocol['E1']['allocation'].values())==protocol['E1']['proposed_attempt_cap']
(root/'configs'/'EXPERIMENT_PROTOCOL.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
print('Wrote two schemas and the proposed independent protocol.')
