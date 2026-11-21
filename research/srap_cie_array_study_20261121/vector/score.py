"""All 48 arm predictions must be sealed before opening the first material truth."""
import csv
from scipy.ndimage import label,binary_fill_holes
from physics import *
LABELS=['eccentric_shell','split_frame','separated_bodies','crossed_bars_layer']
def measures(m,truth):
    norm=np.linalg.norm(truth);diff=m-truth;half=.5*np.max(truth);a=(m>=half).reshape((12,)*3);b=(truth>=half).reshape(a.shape);inside=truth>0
    return dict(rel_L2=float(np.linalg.norm(diff)/norm),half_IoU=float((a&b).sum()/max((a|b).sum(),1)),inside_L2=float(np.linalg.norm(diff[inside])/norm),background_L2=float(np.linalg.norm(diff[~inside])/norm),foreground_mean_ratio=float(m[b.flatten()].mean()/truth[b.flatten()].mean()),peak_chi=float(m.max()),true_peak_chi=float(truth.max()),components=label(a)[1],true_components=label(b)[1],closed_cavities=label(binary_fill_holes(a)&~a)[1],true_closed_cavities=label(binary_fill_holes(b)&~b)[1])
def score():
    start=time.perf_counter();cpu=time.process_time()
    paths={(i,g,m):ROOT/'predictions'/f'{g}_{i:02d}_{m}' for i in range(12) for g in ['plane_grid','hemi'] for m in ['Dual','SRAP']}
    for p in paths.values():
        if not (p/'SEALED.json').is_file():raise RuntimeError(f'Missing {p.name} seal')
        with np.load(p/'materials.npz') as z:
            assert set(z.files)=={'raw_60','raw_180','light_60','light_180'} and all(np.isfinite(z[k]).all() and z[k].size==1728 for k in z.files)
    write(ROOT/'COHORT_SEAL.json',dict(arms=[p.name for p in paths.values()],unique_arrays=48,selected_rows=192,checked_at_unix=time.time(),scientific_raw_and_corrector_only=True))
    mods={g:Model(12,g) for g in ['plane_grid','hemi']};rows=[]
    for (i,g,arm),p in paths.items():
        truth=np.load(ROOT/'sealed_truth'/f'{i:02d}.npz')['chi'];mod=mods[g];Y=t(np.load(ROOT/'data'/f'{g}_{i:02d}.npz')['Y']);info=json.loads((p/'SEALED.json').read_text())
        with np.load(p/'materials.npz') as z:
            for key in z.files:
                m=z[key];view,num=key.split('_');budget=int(num);chosen=info['selected'][key];seed=chosen if view=='raw' else chosen['seed'];ph,res=mod.physical(t(m).real,Y)
                corr=next(a for a in info['correctors'] if a['budget']==budget)
                rows.append(dict(object=i,family=LABELS[i//3],array=g,method=arm,view=view,budget=budget,**measures(m,truth),data_rel=res,physical_phi=ph,paid_phi=chosen['phi'],seed_step=seed['step'],seed_cycles=seed['cycles'],seed_elapsed=seed['time'],corrector_accepted=corr['num_accepted'] if view=='light' else 0,corrector_elapsed=corr['seconds'] if view=='light' else 0,raw_trajectory_seconds=info['raw_s'],all_arm_seconds=info['total_s'],peak_gpu_bytes=info['peak_gpu_bytes'],route_rank=info['route'].get('route_rank',0),route_setup_s=info['route'].get('setup_s',0)))
    with (ROOT/'objects.csv').open('w',newline='') as fd:
        w=csv.DictWriter(fd,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    write(ROOT/'scoring_cost.json',dict(wall=time.perf_counter()-start,cpu=time.process_time()-cpu,by_array={g:dict(m.count) for g,m in mods.items()},max_phi_discrepancy=max(abs(r['physical_phi']-r['paid_phi']) for r in rows)))
if __name__=='__main__':score()
