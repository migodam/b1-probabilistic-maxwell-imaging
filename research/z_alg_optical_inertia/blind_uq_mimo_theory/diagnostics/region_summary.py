"""Post-freeze score all stable ROIs. Labels appear only in the scoring block."""
import csv,json
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parents[1]; res=P/'results'
def read(f):return json.loads(f.read_text())
def csvfile(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def summary():
    bounds={v['case']:v for v in read(P/'bounds/results/summary.json')['stationary']}
    allrows=[];subjects=[]; catalog=[]
    for case,br in sorted(bounds.items()):
        base=P/'imaging/results'/case/'stationary';a=np.load(base/'baseline.npz')
        rois=a['rois'].reshape(len(a['rois']),-1);C=len(rois)
        candidates=[];continuous=[None]*C
        def material(x,src,ref,upper):
            if upper is None or not 0<=upper<=1 or min(x)<0 or max(x)>4.5:return
            peaks=[float(max(x[r])) for r in rois]
            candidates.append((src,ref,float(upper),peaks));catalog.append(dict(case=case,file=src,replay=ref,upper_over_r=upper,peaks=peaks))
        for kind in ('stable_profiles','path_profiles','gn_profiles','gap_profiles'):
            for f in sorted((res/kind/case).glob('roi_??_inner.json')):
                b=read(f);i=b['roi']-1;ff=f.with_name(f.name.replace('_inner.json','.npz'))
                with np.load(ff) as z:x=z['material'].ravel();mask=z['roi'].ravel()
                if not np.array_equal(mask,rois[i]):raise ValueError(f'Changed final mask {ff}')
                material(x,str(ff.relative_to(P)),str(f.relative_to(P)),b['base_residual_ratio_upper'])
                sel=b['selected']
                if b['strict_feasible_base'] and sel:
                    high=x.copy();high[sel['pixel']]=sel['cell_to']
                    material(high,str(ff.relative_to(P))+' + single cell '+str(sel['pixel'])+'='+str(sel['cell_to']),str(f.relative_to(P)),sel['end_ratio_upper'])
                    if sel['full_segment_strict'] and (continuous[i] is None or sel['peak_to']-sel['peak_from']>continuous[i][1]-continuous[i][0]):continuous[i]=(sel['peak_from'],sel['peak_to'])
            for f in (res/kind/case).glob('roi_??_pin.json'):
                j=read(f)
                if j['verified']:
                    src=f.with_suffix('.npz')
                    with np.load(src) as z:material(z['material'].ravel(),str(src.relative_to(P)),str(f.relative_to(P)),j['original_residual_ratio_upper'])
        ext=res/'extended'/(case+'.json')
        if ext.exists():
            e=read(ext)
            with np.load(ext.with_suffix('.npz')) as z:material(z['material'].ravel(),str(ext.with_suffix('.npz').relative_to(P)),str(ext.relative_to(P)),e['upper_residual_over_r'])
        for rr in sorted((res/'native_replay'/case).glob('*.json')):
            j=read(rr)
            if j['residual_upper_over_r']<=1:
                with np.load(P/j['file']) as z:material(z[j['key']].ravel(),j['file']+'#'+j['key'],str(rr.relative_to(P)),j['residual_upper_over_r'])
        # Opening material labels and posterior scores only after constructing every physical endpoint.
        with np.load(P/'ports/private'/(case+'.npz')) as z:truth=z['material'].ravel()
        gg=np.load(res/'gaussian'/(case+'.npz'))
        assert np.array_equal(gg['roi_masks'].reshape(C,-1),rois)
        for k,rm in enumerate(rois):
            c=continuous[k];bw=br['profile_windows'][f'roi_{k+1:02}'];outer_lo=bw.get('largest_strict_positive_peak_floor') or 0.
            lo=min(candidates,key=lambda v:v[3][k]) if candidates else None;hi=max(candidates,key=lambda v:v[3][k]) if candidates else None
            gl,_,gh=gg['quantiles'][k];gt=float(truth[rm].max())
            allrows.append(dict(case=case,roi=k+1,cells=int(rm.sum()),outer_lower=outer_lo,outer_upper=4.5,
                 known_attained_min=lo[3][k] if lo else None,known_attained_max=hi[3][k] if hi else None,
                 min_receipt=lo[1] if lo else '',max_receipt=hi[1] if hi else '',min_map=lo[0] if lo else '',max_map=hi[0] if hi else '',
                 continuous_from=c[0] if c else None,continuous_to=c[1] if c else None,
                 width_exact=bool(lo and lo[3][k]==0 and hi[3][k]==4.5),full_continuous=bool(c and c==(0,4.5)),
                 zero_witness=bool(lo and lo[3][k]==0),true_peak=gt,true_mean=float(truth[rm].mean()),true_cells=int((truth[rm]>.25).sum()),
                 covers=bool(outer_lo<=gt<=4.5),gaussian_lower=float(gl),gaussian_upper=float(gh),gaussian_covers=bool(gl<=gt<=gh),
                 compressed_neg_numeric=bw['numeric_smallest_negative_count']))
        b=read(base/'baseline.json');subjects.append(dict(case=case,grid=m if (m:=int(case.split('_')[0][1:])) else None,
            ROI_count=C,certified_maps=len(candidates),baseline_updates=b['updates'],baseline_converged=b['stable'],
            stationary_residual_over_r=b['original_physics']['residual_over_radius'],stationary_cpu_s=b['cpu_s'],stationary_wall_s=b['wall_s'],
            global_peak_lower=br.get('verified_global_peak_lower_after_dense_cap_scan') or br['verified_global_peak_lower'],
            true_peak=float(truth.max())))
    csvfile(res/'region_intervals.csv',allrows);csvfile(res/'objects.csv',subjects)
    (res/'verified_profiles_catalog.json').write_text(json.dumps(catalog,indent=2))
    def counts(rows):
        return dict(objects=len(set(v['case'] for v in rows)),regions=len(rows),zero=sum(v['zero_witness'] for v in rows),
            exact_hull=sum(v['width_exact'] for v in rows),continuous=sum(v['full_continuous'] for v in rows),
            covered=sum(v['covers'] for v in rows),gaussian_covered=sum(v['gaussian_covers'] for v in rows),
            originally_vacuum=sum(v['true_cells']==0 for v in rows),removable_with_real_cells=sum(v['true_cells']>0 and v['zero_witness'] for v in rows))
    s={k:counts([v for v in allrows if v['case'].startswith(k)]) for k in ('n32','n96')};s['all']=counts(allrows)
    (res/'region_summary.json').write_text(json.dumps(s,indent=2));print(s)
if __name__=='__main__':summary()
