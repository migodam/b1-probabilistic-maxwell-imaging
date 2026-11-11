"""Exact copies of four exposed public starts. No new observations or GT."""
from pathlib import Path
import json,shutil,time

ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT.parent/'research_max_struct_r12/supplement_paired'
cases=['A__r10_twod_0004','B__r12_confirm24_0010','C__r8_test_0001','A__r10_twod_0001']
methods=['dr','ap','ap_dr','event_integer','event_fractional','random_event','normal_double','anderson','reset_fractional']


def main():
    start=time.process_time();out=ROOT/'data_public/development';out.mkdir(exist_ok=True)
    for c in cases:shutil.copy2(OLD/'inputs'/(c+'.npz'),out/(c+'.npz'))
    cfg=dict(CPU_cap=120,proposal_fraction=.6,steps=3000,audit_every=50,AP_preparation=128,
        scale=6,omega=1,corrector='inherited original physical L-BFGS-B',selection='minimum paid physical Phi; secondary minR within original CN99',
        material_box=[.02,20],prior=dict(smooth=1e-4,tv=1e-4,epsilon=.01),
        random_seed=13013,primary='all paid physical candidates including original',prefixes=[15,30,60,90,120],
        scope='four already exposed posthoc development cases; not blind confirmation')
    tasks=[dict(task_id=c+'__'+m,case_id=c,method=m,campaign='development_v1',input_file='data_public/development/'+c+'.npz',seed=13013) for c in cases for m in methods]
    data=dict(campaign='development_v1',status='FROZEN_BEFORE_NEW_RUNS',protocol=cfg,tasks=tasks)
    (ROOT/'configs/development_v1.json').write_text(json.dumps(data,indent=2)+'\n')
    (ROOT/'audit/DEVELOPMENT_INPUT_PREPARATION.json').write_text(json.dumps(dict(cases=cases,truth_access='NONE',new_physics=0,CPU_seconds=time.process_time()-start),indent=2)+'\n')
    print(json.dumps(dict(tasks=len(tasks),truth_access='NONE')))


if __name__=='__main__':main()
