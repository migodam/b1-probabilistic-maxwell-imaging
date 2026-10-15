"""Necessary contract checks only. All fixtures are SYNTHETIC, never discoveries."""
from __future__ import annotations
import copy,json,sys,tempfile,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'src'))
from research_a.contracts import (write_packages,validate_native,validate_exchange,
    validate_schema_subset,evaluate_new_action,zero_cost,error_record,CONTRACT_DIR)


def fixture():
    context=dict(package_id='synthetic-contract-fixture-NOT-A-REAL-DISCOVERY',target_id='synthetic-target',
      solver_id='synthetic-no-physics',backend_commit='synthetic-fixture',material_law_id='synthetic-law',
      time_convention='exp(-i omega t)',acquisition=dict(id='synthetic-acq',frequencies=[1.6],frequency_units='kR',n_sources=2),
      observation=np.zeros(4),noise_covariance=np.eye(4),whitener=np.eye(4),noise_id='synthetic-unit-noise',
      material_support={'scope':'synthetic'},nuisance_support={'gain':[.7,1.3]},prior_manifest=None,
      grid={'shape':[2,2],'scope':'SYNTHETIC'},fit_threshold=1.,total_cost=zero_cost(),
      distance_definition_id='synthetic-l2',model_discrepancy_policy='unmodeled, synthetic contract only',
      self_term_convention='synthetic',polarizability_map='synthetic',boundary_conditions='synthetic',
      dispersion_model='synthetic nondispersive',candidate_ordering_rule='creation order, no truth',
      unexplored_regions=['all real materials'],search_restrictions=['SYNTHETIC ONLY'],
      upstream_manifest={'scope':'synthetic contract fixture'},failed_attempts=[],
      new_action_interface=dict(evaluator_entrypoint='caller-supplied explicit callable, no automatic physics',environment={'python':'test'},
        allowed_actions=[{'id':'test-action','frequency':2.4}],cache_manifest={}))
    row=dict(id='synthetic-world-0',material=np.full((2,2),2.5+.1j),prediction=np.full(4,.1),
      whitened_residual_l2=.2,validated_full_maxwell=True,
      validation=dict(constitutive_relative_defect=1e-14,solver_relative_residual=1e-14, material_constraints_passed=True,all_illuminations_same_material=True,
        numerical_error=error_record(value=1e-14,status='solver_tolerance_only',method='SYNTHETIC fixture, not physical')),
      structure=dict(support_volume=1.,cavity_count=0,component_count=1,feature_scale=.1,
        descriptor_definition_id='synthetic',cluster_ids=[],event_ids=[],topology_thresholds=[0.],
        topology_persistence_rule='synthetic not a claim'),
      compensation=dict(edit_family='synthetic',edit_parameters={},material_compensation_parameters={},
        nuisance_compensation_parameters={'gain':1.},coordinates=np.array([0.,0.,2.5,.1])),
      material_parameterization_id='synthetic',material_distance_to_base=0.,cost=zero_cost(),
      nuisance_parameters={'gain':1.},history=[{'status':'SYNTHETIC fixture only'}])
    return context,row


class ContractChecks(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.ctx,self.row=fixture()
        self.result=write_packages(self.root,self.ctx,[self.row])
    def tearDown(self): self.tmp.cleanup()
    def native(self): return copy.deepcopy(self.result['manifest'])
    def exchange(self): return copy.deepcopy(self.result['exchange'])
    def test_positive_native_and_common(self):
        self.assertEqual(validate_native(self.native(),self.root)['semantic_checks'],'PASS')
        self.assertEqual(validate_exchange(self.exchange(),self.root)['physics_evaluations'],0)
    def test_original_digest_requirement_disclosed(self):
        original=json.loads((CONTRACT_DIR/'original/MaterialHypothesisPackage.schema.json').read_text())
        with self.assertRaises(ValueError): validate_schema_subset(self.native(),original)
        self.assertIsNone(self.native()['observation']['sha256'])
    def test_candidate_only_probability_rejected(self):
        d=self.native();d['candidates'][0]['posterior_weight']=.5
        with self.assertRaises(ValueError):validate_native(d,self.root)
    def test_target_mismatch_rejected(self):
        d=self.native();d['candidates'][0]['validation']['physics_target_hash']='wrong'
        with self.assertRaises(ValueError):validate_native(d,self.root)
    def test_failed_fit_rejected(self):
        d=self.native();d['candidates'][0]['fit'].update(whitened_residual_l2=2.,squared_residual=4.)
        with self.assertRaises(ValueError):validate_native(d,self.root)
    def test_wrong_k_rejected(self):
        d=self.native();d['valid_k']=2
        with self.assertRaises(ValueError):validate_native(d,self.root)
    def test_truth_leak_rejected(self):
        d=self.native();d['candidates'][0]['compensation']['edit_parameters']['truth_mae']=.1
        with self.assertRaises(ValueError):validate_native(d,self.root)
    def test_fake_information_gain_rejected(self):
        d=self.native();d['next_measurement_interface']['posterior_information_gain_allowed']=True
        with self.assertRaises(ValueError):validate_native(d,self.root)
    def test_prediction_replay_rejected(self):
        d=self.exchange();d['material_worlds'][0]['fit']['whitened_residual_l2']=.1
        d['material_worlds'][0]['fit']['squared_residual']=.01
        with self.assertRaises(ValueError):validate_exchange(d,self.root)
    def test_asset_path_escape_rejected(self):
        d=self.exchange();d['material_worlds'][0]['full_material_array']['path']='../elsewhere.npy'
        with self.assertRaises(ValueError):validate_exchange(d,self.root)
    def test_no_candidates_is_not_real_handoff(self):
        with self.assertRaises(ValueError):write_packages(self.root,self.ctx,[])
    def test_missing_physical_defect_is_not_fabricated(self):
        row=copy.deepcopy(self.row);row['validation'].pop('solver_relative_residual')
        with self.assertRaises(ValueError):write_packages(self.root,self.ctx,[row])
    def test_analytic_mie_unexposed_constitutive_is_null(self):
        row=copy.deepcopy(self.row)
        row['validation'].update(constitutive_relative_defect=None,constitutive_relative_defect_status='NOT_EXPOSED_BY_ANALYTIC_MIE_INTERFACE')
        result=write_packages(self.root,self.ctx,[row])
        self.assertIsNone(result['manifest']['candidates'][0]['validation']['constitutive_relative_defect'])
    def test_fabricated_constitutive_zero_rejected(self):
        row=copy.deepcopy(self.row)
        row['validation'].update(constitutive_relative_defect=0.,constitutive_relative_defect_status='NOT_EXPOSED_BY_ANALYTIC_MIE_INTERFACE')
        with self.assertRaises(ValueError):write_packages(self.root,self.ctx,[row])
    def test_no_new_action_without_budget(self):
        called=[]
        def callback(*args):called.append(True)
        with self.assertRaises(ValueError):evaluate_new_action(self.exchange(),self.root,self.row['id'],'test-action',{},callback)
        self.assertEqual(called,[])
    def test_explicit_callback_contract_without_physics(self):
        called=[]
        def callback(material,action,budget):
            called.append(material.shape)
            return dict(target_identity=self.ctx['target_id'],action_id=action['id'],
                realification_order='all-real-then-all-imaginary',prediction=[0.,0.],cost_receipt={'physics_attempts':0,'scope':'synthetic callback'},validation_status='SYNTHETIC-NO-PHYSICS')
        token=dict(authorized=True,remaining_attempts=1,ledger_id='synthetic-contract-only',target_identity=self.ctx['target_id'])
        out=evaluate_new_action(self.exchange(),self.root,self.row['id'],'test-action',token,callback)
        self.assertEqual(called,[(2,2)]);self.assertEqual(out['cost_receipt']['physics_attempts'],0)


if __name__=='__main__':
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(ContractChecks)
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    report=dict(scope='SYNTHETIC interface/semantic checks only; no actual candidate discoveries',
        passed=result.testsRun-len(result.errors)-len(result.failures),tests_run=result.testsRun,
        failures=[str(x[0])+': '+x[1] for x in result.failures],errors=[str(x[0])+': '+x[1] for x in result.errors],
        new_physics_evaluations=0,digest_computation='NOT_RUN_BY_USER_POLICY',
        original_validator_dependency='jsonschema unavailable in existing runtime; original files preserved unchanged; equivalent required semantics checked by dependency-free schema-subset validator')
    (CONTRACT_DIR/'ZERO_PHYSICS_CONTRACT_CHECKS.json').write_text(json.dumps(report,indent=2)+'\n')
    raise SystemExit(0 if result.wasSuccessful() else 1)
