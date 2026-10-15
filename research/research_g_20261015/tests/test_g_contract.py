"""Sidecar software contracts; fixtures are not physical imaging evidence."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from g_contract import validate_episode, ContractError
ROOT=Path(__file__).resolve().parents[1]
class SidecarContract(unittest.TestCase):
    def setUp(self):self.example=json.loads((ROOT/'examples/GSearchEpisode.minimal.json').read_text())
    def reject(self,change):
        item=copy.deepcopy(self.example);change(item)
        with self.assertRaises(ContractError):validate_episode(item)
    def test_minimal_and_schema_required_fields(self):
        self.assertTrue(validate_episode(self.example))
        schema=json.loads((ROOT/'GSearchEpisode.schema.json').read_text())
        self.assertEqual(set(schema['required']),set(self.example))
        self.assertEqual(schema['properties']['version']['const'],self.example['version'])
        for key in self.example:self.reject(lambda x,k=key:x.pop(k))
    def test_wrong_target_packing_gain_and_episode_versions(self):
        self.reject(lambda x:x.update(version='GSearchEpisode/v2'))
        for key in ['id','packing_version','gain_version','gain']:
            self.reject(lambda x,k=key:x['target'].update({k:'wrong'}))
        self.reject(lambda x:x['target'].update(gain=True))
    def test_string_paths_and_array_references(self):
        for key in self.example['candidate']:
            if key.endswith('_path'):self.example['candidate'][key]='arrays/field.npy'
        self.assertTrue(validate_episode(self.example))
        for value in ['/absolute.npy','../escape.npy','arrays/../../escape','C:/file.npy',{'path':'field.npz','array_key':''}]:
            self.reject(lambda x,v=value:x['candidate'].update(raw_field_path=v))
    def test_missing_diagnostics_need_reason_and_no_nan_bool(self):
        for key in ['solver_diagnostics','constitutive_diagnostics','full_residual']:
            for value in [{'value':None,'reason':None},{'value':1.,'reason':'claimed'},{'value':True,'reason':None},{'value':float('nan'),'reason':None}]:
                self.reject(lambda x,k=key,v=value:x['validation'].update({k:v}))
            self.example['validation'][key]={'value':1e-12,'reason':None}
        self.assertTrue(validate_episode(self.example))
    def test_validation_required(self):
        for key in ['full_residual','material_legal','full_validation_receipt']:
            self.reject(lambda x,k=key:x['validation'].pop(k))
    def test_mass_is_not_invented(self):
        for key in ['posterior_weight','represented_mass','unmodeled_mass']:
            self.reject(lambda x,k=key:x['probability'].update({k:0.5}))
        self.example['probability']['proposal_density']={'value':.2,'interpretation':'PROPOSAL_ONLY'}
        self.assertTrue(validate_episode(self.example))
        self.reject(lambda x:x['probability']['proposal_density'].update(interpretation='POSTERIOR'))
    def test_forbidden_online_inputs_recursive(self):
        for key in ['truth','oracle','hole_count','topology_label','true_truth_path','holeCount']:
            self.reject(lambda x,k=key:x['history'].update(nested=[{k:'forbidden'}]))
    def test_duplicate_observation_id(self):
        self.reject(lambda x:x['observations'].append(copy.deepcopy(x['observations'][0])))
    def test_asset_presence_and_npz_keys(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            refs=self.example['native_package_refs']+[o['data_path'] for o in self.example['observations']]+[v for k,v in self.example['candidate'].items() if k.endswith('_path')]+[self.example['validation']['full_validation_receipt']]
            for ref in refs:
                path=root/(ref if isinstance(ref,str) else ref['path']);path.parent.mkdir(parents=True,exist_ok=True)
                if isinstance(ref,str):path.write_text('{}')
            np.savez(root/'arrays/example.npz',**{ref['array_key']:np.zeros(1) for ref in refs if isinstance(ref,dict)})
            self.assertTrue(validate_episode(self.example,root,True))
            bad=copy.deepcopy(self.example);bad['candidate']['raw_field_path']['array_key']='missing'
            with self.assertRaisesRegex(ContractError,'ARRAY_KEY_MISSING'):validate_episode(bad,root,True)
            (root/'native/package.json').unlink()
            with self.assertRaisesRegex(ContractError,'ASSET_MISSING'):validate_episode(self.example,root,True)
if __name__=='__main__':unittest.main()
