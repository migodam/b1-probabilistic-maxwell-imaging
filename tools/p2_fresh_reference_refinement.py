"""All fresh observations' moment/event reference refinement, zero solves."""
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.b1_joint.response import Response
from src.b1_joint.reference import Reference
from src.b1_joint.material import RMAX
from src.b1_joint.accounting import dump
def main():
    out=ROOT/'research/joint_p2/results';d=np.load(out/'FRESH_FINAL_OBSERVATIONS.npz');ref=Reference(Response(out/'COMMON_RESPONSE_BANK.npz'));width=np.array([*RMAX,.6,.12]);rows=[]
    for i,y in enumerate(d['y']):
        a=ref.solve(y[:32],d['sigma'][:32],order=13);b=ref.solve(y[:32],d['sigma'][:32],order=17)
        error=dict(event=float(abs(a['event_mass']-b['event_mass']).max()),mean_fraction=float((abs(a['mean']-b['mean'])/width).max()),SD_relative=float((abs(a['SD']-b['SD'])/np.maximum(b['SD'],width*1e-6)).max()))
        passed=error['event']<=.01 and error['mean_fraction']<=.001 and error['SD_relative']<=.01
        row=dict(observation=i,parent=int(d['parent'][i]),original13_refinement=error,original13_pass=passed,order17_event_mass=b['event_mass'],order17_mean=b['mean'],order17_SD=b['SD'])
        if not passed:
            c=ref.solve(y[:32],d['sigma'][:32],order=21);e=dict(event=float(abs(b['event_mass']-c['event_mass']).max()),mean_fraction=float((abs(b['mean']-c['mean'])/width).max()),SD_relative=float((abs(b['SD']-c['SD'])/np.maximum(c['SD'],width*1e-6)).max()));row.update(order17_21_refinement=e,order17_21_pass=e['event']<=.01 and e['mean_fraction']<=.001 and e['SD_relative']<=.01,order21_event_mass=c['event_mass'])
        rows.append(row)
    dump(out/'ALL_FRESH_REFERENCE_REFINEMENT.json',dict(observations=len(rows),rows=rows,failed_original13_observations=[a['observation'] for a in rows if not a['original13_pass']],unresolved_after21=[a['observation'] for a in rows if not a.get('order17_21_pass',True)],target_unchanged=True,old_metrics_retained=True,scope='event/moments refinement only; CDF partition and independent direct quantiles separate',new_physics=0))
if __name__=='__main__':main()
