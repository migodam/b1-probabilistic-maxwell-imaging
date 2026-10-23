"""Read saved arrays/source only; never imports wq or performs any PDE action."""
from pathlib import Path
import csv
import json
import subprocess
import numpy as np
from scipy.ndimage import label

OLD = Path("/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/Research_W/research-w-continuation")
REPO = Path("/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-w-p0r-to-p1")
OUT = Path(__file__).resolve().parent
COMMIT = '483d260c'
rel = 'research/research_w_continuation/'
provenance = {}
for name in ['src/wq/search.py', 'src/wq/inverse_support.py', 'src/wq/fem.py', 'src/wq/experiments.py']:
    committed = subprocess.run(['git', '-C', str(REPO), 'show', f'{COMMIT}:{rel}{name}'], check=True, stdout=subprocess.PIPE).stdout
    provenance[name] = dict(original_equals_fixed_commit=(OLD/name).read_bytes() == committed,
                            worktree_equals_fixed_commit=(REPO/rel/name).read_bytes() == committed)
rows = json.loads((OLD/'results/p0b/search-corrected/RESULTS.json').read_text())
n = 48*32
scale = np.r_[np.full(n,8.),np.full(n,.15)]
x = np.linspace(-2,2,49); y=np.linspace(0,3,33)
wx=np.maximum(0,np.minimum(x[1:],1.2)-np.maximum(x[:-1],-1.2))
wy=np.maximum(0,np.minimum(y[1:],2.8)-np.maximum(y[:-1],1.6))
areas=12/n-wy[:,None]*wx[None,:]
stats=[]
with np.load(OLD/'results/p0b/search-corrected/ARRAYS.npz') as z:
    for r in rows:
        if r['case']!='calibration_candidate':continue
        key=f"{r['seed']}_{r['band']}_{r['method']}"
        theta=z[key]; common=z[f"{r['seed']}_{r['band']}_COMMON"]
        active=((theta[:n]!=0)|(theta[n:]!=0)).reshape(32,48)
        base=((common[:n]!=0)|(common[n:]!=0)).reshape(32,48)
        cc,count=label(active)
        sizes=np.bincount(cc.ravel())[1:]
        with np.load(OLD/'results/wall-pml/ARRAYS.npz') as truth:
            dcal=np.linalg.norm((z[key+'_data'].ravel()-truth[f"seed{r['seed']}_{r['band']}_cal_volume"].ravel()))/.01
        stats.append(dict(seed=r['seed'],band=r['band'],method=r['method'],d_cal=r['d_cal'],
                          replay_d_cal=float(dcal),accepted_updates=r['accepted_updates'],method_rhs=r['method_rhs'],
                          active_pixels=int(active.sum()),active_area=float(areas[active].sum()),
                          exact_same_as_common=bool(np.array_equal(theta,common)),
                          support_added=int((active&~base).sum()),support_removed=int((base&~active).sum()),
                          components_4_neighbor=int(count),singleton_components=int((sizes==1).sum()),
                          normalized_parameter_linf=float(np.abs(theta/scale).max()),
                          contrast_min_active=float(theta[:n][active.ravel()].min()),
                          contrast_max=float(theta[:n].max())))
with (OUT/'saved_candidate_diagnostics.csv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=list(stats[0]));writer.writeheader();writer.writerows(stats)
with np.load(OLD/'results/p0b/preflight/ARRAYS.npz') as p:
    keys=list(p.files)
    preflight=dict(keys=keys,shapes={k:list(p[k].shape) for k in keys})
    if all(k in p for k in ['theta','direction','Jdirection','FD']):
        t=p['theta'];v=p['direction']
        preflight.update(relative_error_replayed=float(np.linalg.norm(p['FD']-p['Jdirection'])/np.linalg.norm(p['Jdirection'])),
                         nonzero_parameter_count=int(np.count_nonzero(t)),
                         direction_outside_active_parameter_count=int(np.count_nonzero(v[t==0])))
ledger=[json.loads(s) for s in (OLD/'results/p0b/PHYSICS_LEDGER.jsonl').read_text().splitlines()]
reserves=[r for r in ledger if r.get('event')=='reserve']
pid_stats={}
for r in ledger:
    pid_stats.setdefault(str(r.get('pid')), {'events':0,'reserved_rhs':0,'backend_states':0,'ndofs':set()})
    q=pid_stats[str(r.get('pid'))];q['events']+=1
    if r.get('event')=='reserve':q['reserved_rhs']+=r.get('rhs',0)
    if r.get('event')=='backend_state':q['backend_states']+=1
    if 'ndofs' in r:q['ndofs'].add(r['ndofs'])
for q in pid_stats.values():q['ndofs']=sorted(q['ndofs'])
summary=dict(physics_actions_in_this_audit=0, historical_commit=COMMIT, provenance=provenance,
             source_root=str(OLD),preflight_saved_array_check=preflight,
             historical_ledger_by_pid=pid_stats,
             candidates=stats,
             note='All arithmetic uses saved outputs; no original result file is written. Component counts are descriptive, not a minimum-feature certification.')
(OUT/'saved_only_audit.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(dict(physics_actions=0,source_match=all(v['original_equals_fixed_commit'] for v in provenance.values()),
                     candidates=len(stats),max_saved_dcal_replay_error=max(abs(r['d_cal']-r['replay_d_cal']) for r in stats),
                     preflight=preflight,output=str(OUT)),indent=2))
