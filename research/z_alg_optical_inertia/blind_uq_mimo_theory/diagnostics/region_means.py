"""Compute attainable ROI mean differences from previously verified full-wave maps."""
import csv
import json
import math
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from time import perf_counter, process_time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def exact_mean(values):
    # Stored real binary64 materials are dyadic. Sum and divide before rounding.
    return sum((Fraction.from_float(float(value)) for value in values), Fraction()) / len(values)


def one_sided(value, upper):
    nearest = float(value)
    offset = Fraction.from_float(nearest) - value
    if (upper and offset < 0) or (not upper and offset > 0):
        return math.nextafter(nearest, math.inf if upper else -math.inf)
    return nearest


def read_material(source):
    base, sep, patch = source.partition(' + single cell ')
    rel, _, key = base.partition('#')
    with np.load(ROOT / rel, allow_pickle=False) as data:
        material = data[key or 'material'].ravel().copy()
    if sep:
        pixel, contrast = patch.split('=')
        material[int(pixel)] = float(contrast)
    return material


def main():
    tic, cpu = perf_counter(), process_time()
    archive = ROOT / 'results'
    by_case = defaultdict(list)
    for rec in json.loads((archive / 'verified_profiles_catalog.json').read_text()):
        if not 0 <= rec['upper_over_r'] <= 1:
            raise ValueError('The catalog contains an unverified map.')
        by_case[rec['case']].append(rec)
    out = []
    for case, rows in sorted(by_case.items()):
        with np.load(ROOT / 'imaging/results' / case / 'stationary/baseline.npz') as a:
            rois = a['rois'].reshape(len(a['rois']), -1)
        values = [[] for _ in rois]
        for row in rows:
            x = read_material(row['file'])
            for i, mask in enumerate(rois):
                if float(np.max(x[mask])) != row['peaks'][i]:
                    raise ValueError('Peak mismatch reading an accepted material.')
                values[i].append((exact_mean(x[mask]), row))
        for i, means in enumerate(values):
            low, lowrec = min(means,key=lambda pair:pair[0])
            high, highrec = max(means,key=lambda pair:pair[0])
            out.append(dict(case=case,roi=i+1,cells=int(rois[i].sum()),
                guaranteed_mean_outer_low=0,guaranteed_mean_outer_high=4.5,
                low_map_mean_upper=one_sided(low,True), high_map_mean_lower=one_sided(high,False),
                proven_mean_width_lower=one_sided(high-low,False),low_map=lowrec['file'],high_map=highrec['file'],
                low_replay=lowrec['replay'],high_replay=highrec['replay']))
    with (archive/'region_means.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
    def summarize(rows):
        widths=[r['proven_mean_width_lower'] for r in rows]
        return dict(regions=len(rows),nonzero=sum(x>0 for x in widths),at_least_one=sum(x>=1 for x in widths),
                    at_least_two=sum(x>=2 for x in widths),minimum=min(widths),median=float(np.median(widths)),maximum=max(widths))
    summary={k:summarize([r for r in out if r['case'].startswith(k)]) for k in ('n32','n96','')}
    summary['all']=summary.pop('')
    summary.update(wall_s=perf_counter()-tic,cpu_s=process_time()-cpu,source='Verified maps only, no truth or held-out data')
    (archive/'region_means_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False))
    print(json.dumps(summary,ensure_ascii=False))


if __name__ == '__main__':
    main()
