"""OFFLINE ONLY scoring. Never imported by the online inverse search."""
from pathlib import Path
import csv,json
import numpy as np
from scipy.ndimage import label,binary_fill_holes
from wq.experiments import continuous
from wq.inverse_support import PixelGrid,pixel_material
root=Path(__file__).resolve().parents[1];g=PixelGrid();r=json.loads((root/'results/p0b/search-corrected/RESULTS.json').read_text());t=json.loads((root/'results/p0b/task-endpoints/RESULTS.json').read_text());a=np.load(root/'results/p0b/search-corrected/ARRAYS.npz');b=np.load(root/'results/p0b/task-endpoints/ARRAYS.npz');ref=np.load(root/'results/wall-pml/ARRAYS.npz')
xx,yy=np.meshgrid(-2+(np.arange(768)+.5)*4/768,(np.arange(576)+.5)*3/576);pts=np.c_[xx.ravel(),yy.ravel()];rows=[]
for row in r:
    if row['case']!='calibration_candidate':continue
    seed,band,method=row['seed'],row['band'],row['method'];key=f'{seed}_{band}_{method}';theta=a[key];estimated=pixel_material(theta,1.,g)(pts);truth=continuous(seed,1.)(pts)
    mask=(estimated.real>1.1).reshape(yy.shape);tm=(truth.real>1.1).reshape(yy.shape)
    ncomp=label(mask,np.ones((3,3)))[1];holes=label(binary_fill_holes(mask)&~mask)[1]
    targetrow=next(x for x in t if x['seed']==seed and x['band']==band and x['method']==method)
    dt=b[key+'_task_data'];dc=a[key+'_data'];rt=ref[f'seed{seed}_{band}_task_volume'];rc=ref[f'seed{seed}_{band}_cal_volume']
    rows.append(dict(**{k:row[k] for k in ('seed','band','method','d_cal','method_rhs','method_wall_seconds','accepted_updates','occupied_area')},fixed_probe8_total_distance=targetrow['task8_response_distance'],fixed_probe8_induced_distance=float(np.linalg.norm((dt-dc)-(rt-rc))/.01),material_L2_error=float(np.linalg.norm(estimated-truth)*np.sqrt(12/len(pts))),support_IoU=float(np.sum(mask&tm)/max(np.sum(mask|tm),1)),components=int(ncomp),holes=int(holes),truth_scoring='OFFLINE_ONLY_AFTER_SAVING_ONLINE_OUTPUT',independent_parent_count=2,task_class_coverage='1/24',joint_tested=False,standalone_proposal_rhs=row['method_rhs']+(160 if method=='TASK_AWARE' else 80)))
out=root/'results/p0b'
with (out/'P0_B_CALIBRATION_RESULTS.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
(out/'SCORED_RESULTS.json').write_text(json.dumps(rows,indent=2))
with (out/'P0_B_JOINT_RESULTS.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=['status','reason']);w.writeheader();w.writerow(dict(status='NOT_RUN_NO_CALIBRATION_WITNESS',reason='No observed calibration distance <=1; task8 uses same target and is not joint compensation'))
print(json.dumps(dict(scored=len(rows),minimum_dcal=min(x['d_cal'] for x in rows),oracle_selection_used_online=False)))
