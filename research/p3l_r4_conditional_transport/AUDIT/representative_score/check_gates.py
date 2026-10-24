"""Pure schema/gate checks with synthetic files; no actualjob claim."""
from pathlib import Path
import json,time,os,numpy as np
from p3r4.representative_eval import check_representative_gates,action_properties
base=Path('research/p3l_r4_conditional_transport/AUDIT/representative_score/SCHEMA_ONLY');base.mkdir(exist_ok=True)
physical=base/'fake_job';physical.mkdir(exist_ok=True);book=base/'II/jobs/fake_job';book.mkdir(parents=True,exist_ok=True)
mu=base/'mock_estimator.npz';np.savez(mu,parent=np.array(['R4DEV_00000']),prediction_raw=np.full((1,576),.1+.04j));s=mu.stat()
now=time.time_ns();prep=base/'FROZEN_REQUEST.json'
prepared=dict(status='FROZEN_READY',GT_read=False,request=dict(parents=['R4DEV_00000'],estimator_id='mock_estimator',frozen_estimator=True,no_GT=True),source_predictions=dict(path=str(mu.resolve()),size=s.st_size,mtime_ns=s.st_mtime_ns),qualified=[],skips=[])
source_gate=base/'ONLINE_COMPLETE.json';source_gate.write_text(json.dumps(dict(status='ONLINE_COMPLETE',truth_read=False,pid=os.getpid()+100,completed_time_ns=now,freeze=dict(freeze_time_ns=now-10000),parents=['R4DEV_00000'])))
prep.write_text(json.dumps(prepared));os.utime(prep,ns=(now+10000,now+10000));(physical/'FROZEN_REQUEST.json').write_text(json.dumps(prepared));os.utime(physical/'FROZEN_REQUEST.json',ns=(now+20000,now+20000))
(book/'START.json').write_text(json.dumps(dict(kind='physics',job='fake_job',stage='II',pid=os.getpid()+200,start_unix=(now+15000)/1e9)))
(physical/'RESULTS.json').write_text(json.dumps(dict(status='COMPLETE',GT_read=False)));os.utime(physical/'RESULTS.json',ns=(now+30000,now+30000))
receipt=book/'CLOSED.json';receipt.write_text(json.dumps(dict(status='COMPLETE')));os.utime(receipt,ns=(now+40000,now+40000))
check_representative_gates(physical,prep,source_gate,receipt)
receipt.write_text(json.dumps(dict(status='STOPPED')))
try:check_representative_gates(physical,prep,source_gate,receipt)
except ValueError:pass
else:raise AssertionError('unclosedjob accepted')
cfg=json.loads(Path('research/p3l_r4_conditional_transport/configs/TARGET.json').read_text());pts=json.loads(Path('research/p1_material_target_20261018/inputs/geometry.json').read_text())['points'];roi=np.flatnonzero(np.linalg.norm(pts,axis=1)<=cfg['roi_radius'])
a=np.full(576,complex(*cfg['host_chi']));prop=action_properties(a,cfg,roi);assert prop['decoder_realizability']=='NOT_PROVEN' and prop['passive_imag_nonnegative']
a[roi[0]]=-1j;assert not action_properties(a,cfg,roi)['passive_imag_nonnegative']
Path('research/p3l_r4_conditional_transport/AUDIT/representative_score/SCHEMA_CHECK.json').write_text(json.dumps(dict(status='PASS',scope='synthetic chronologyfiles only, no actual representative outputs/scientific claims',closedcomplete_required=True,actualinference_link_required=True,rawbounds_do_not_prove_decoder=True,physics_RHS=0,training_updates=0),indent=2)+'\n')
print('representative gate/materialschema PASS')
