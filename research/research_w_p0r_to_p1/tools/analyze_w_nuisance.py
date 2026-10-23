"""Declared shared-gain stress on fixed, already paid full-Maxwell pairs.

The reference gain stays1. One bounded complex gain may change on world1.
This is profiled identifiability only, never Bayesian evidence or new physics.
"""
import csv,json,time
from pathlib import Path
import numpy as np

root=Path(__file__).resolve().parents[1]

def profile(f,y):
    f=np.asarray(f).ravel();y=np.asarray(y).ravel()
    cross=np.vdot(y,f);norm=float(np.vdot(f,f).real)
    phases=[-.01,.01,float(np.clip(-np.angle(cross),-.01,.01))]
    answers=[]
    for phase in phases:
        a=float(np.clip(np.real(np.exp(1j*phase)*cross)/norm,.99,1.01))
        answers.append((float(np.linalg.norm(a*np.exp(1j*phase)*f-y)/.01),a,phase))
    return min(answers)

def main():
    start=time.process_time();cfg=json.loads((root/'configs/ORIGINAL_P0.json').read_text())
    assert cfg['nuisance']['stress_gain_amplitude']==[.99,1.01]
    assert cfg['nuisance']['stress_shared_phase_rad']==[-.01,.01]
    assert cfg['nuisance']['per_frequency_free_gain'] is False
    records=json.loads((root/'results/C/phase_diagram/PHASE_TABLE.json').read_text());rows=[]
    for r in records:
        f0=np.load(r['response0'])['data'];f1=np.load(r['response1'])['data']
        distance,a,p=profile(f1,f0)
        rows.append(dict(case=r['case'],band=r['band'],wall_loss_scale=r['common_wall_loss_scale'],
            primary_gain1_distance=float(np.linalg.norm(f1-f0)/.01),profiled_one_shared_gain_distance=distance,
            amplitude=a,phase=p,reference_gain=1,at_gain_boundary=bool(a in [.99,1.01] or p in [-.01,.01]),
            interpretation='FIXED_PAIR_BOUNDED_NUISANCE_PROFILE_NOT_POSTERIOR_EVIDENCE',source0=r['response0'],source1=r['response1']))
    out=root/'results/C/nuisance_stress';out.mkdir(exist_ok=True)
    (out/'RESULTS.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
    with (out/'TABLE.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (out/'COST.json').write_text(json.dumps(dict(new_rhs=0,new_LU=0,pure_analysis_cpu_seconds=time.process_time()-start,
        original_saved_full_endpoint_cost_not_free=True,prior_odds_not_computed=True),indent=2)+'\n')
    print(json.dumps([r for r in rows if r['case']=='joint_contrast3' and r['wall_loss_scale']==1]))

if __name__=='__main__':main()
