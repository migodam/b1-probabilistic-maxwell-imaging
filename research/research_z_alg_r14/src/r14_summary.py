"""Combine the saved designs, acquisition records, and truth scores for reporting."""
import csv
import json
from pathlib import Path
import numpy as np
from flint import arb

from data14 import ROOT, load_cases
from measure_acquire import archived_generator
from measure_validation import load_validation_cases, private_trials
from material_controls import GRID,HI,LO

R=ROOT/'results'

def read(path):
    return json.loads(path.read_text())

def normal_limits(case,row):
    col=3 if row['coherent'] else 2
    upper=arb(max(x[col] for x in row['terminal_intervals']))
    risk=((upper/(2*arb(case.sigma))).erfc()/2).lower()
    r=arb(row['threshold_2r'])
    return dict(worst_norm_upper=float(upper),conditional_binary_error_lower=np.nextafter(float(risk),-np.inf),
        required_energy_for_disjoint_balls_lower=np.nextafter(float(((r/upper)**2).lower()),-np.inf))

def builds():
    old_records=read(R/'pivot_old.json')['cases']
    witness=read(R/'pivot_diameter_witnesses.json')['cases']
    raw=archived_generator()
    cases=load_cases()
    validation=load_validation_cases();private={t.pair_id:t for t in private_trials()}
    out=[];choices={};seq={}
    for case in (*cases,*validation):
        cid=case.case_id
        historic=cid in old_records
        if historic:
            store=R/'posterior_joint'
            true=raw[cid].truth
            old=np.load(R/'material_covers_seeded'/(cid+'.npz'))['boxes']
            h=read(R/'material_covers_seeded'/(cid+'.json'))
            n=read(R/'frequency_covers'/(cid+'.json'))
            v=read(R/'posterior'/(cid+'.json'))
            u=read(store/(cid+'.json'))
            d=old_records[cid]['methods']['r10_scalar']
            arm=n['selected']['0.5'].copy()
            arm['wall_s']=d['wall_s']+h['wall_s']+n['wall_s']
            choices[cid]=arm
            stage='first_half' if v['first_half']['reached_half'] else 'both_halves'
            seq[cid]=dict(action=('half_split_first' if stage=='first_half' else 'split_0.02'),
                certified=v[stage]['reached_half'],D_upper=v[stage]['D_upper'],
                wall_s=d['wall_s']+h['wall_s']+v['first_half']['wall_s']+
                    (v['both_halves']['wall_s'] if stage=='both_halves' else 0.))
            extra=dict(prior_old_D_upper=n['zero_new_D_upper'],prospective_D_upper=n['actions']['split_0.02']['D_upper'],
                D_lower=witness[cid]['actions']['split_0.02']['witness']['rms_lower'],
                connectedness_witness=bool(read(R/'disconnected_old.json').get(cid)),old_time=d['wall_s']+h['wall_s'])
        else:
            store=R/'validation_regions'
            z=read(store/(cid+'.json'));u=z['posterior'];v=u
            true=private[cid.rsplit('-',1)[0]].material
            old=np.load(store/(cid+'.npz'))['old']
            extra=dict(old_time=z['r10']['wall_s']+z['old']['wall_s'])
        def contains(box):
            material=LO+(HI-LO)*box/GRID
            return bool(np.any(np.all((material[:,0]<=true)&(true<=material[:,1]),axis=1)))
        with np.load(store/(cid+'.npz')) as boxes:
            out.append(dict(case_id=cid,cohort=case.cohort,old_contains=contains(old),
                    first_contains=contains(boxes['first_half']),both_contains=contains(boxes['both_halves']),
                    first_upper=v['first_half']['D_upper'],both_upper=v['both_halves']['D_upper'],
                    exact_first_upper=u['first_half']['D_upper'],exact_both_upper=u['both_halves']['D_upper'],**extra))
    paths=[(R/'region_actions.json',choices),(R/'sequential_actions.json',seq)]
    for path,payload in paths:path.write_text(json.dumps(payload,indent=2)+'\n')
    keys=list(dict.fromkeys(k for x in out for k in x))
    with (R/'case_summary.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=keys,lineterminator="\n");w.writeheader();w.writerows(out)
    inputs={c.case_id:c for c in cases}
    bands={}
    for file in sorted((R/'continuous_band').glob('*.json')):
        r=read(file)
        if r['certified_continuum']:
            bands[file.stem]=dict(rms=r['rms_lower'],coherent=r['coherent'],cpu_s=r['cpu_s'],
                intervals=len(r['terminal_intervals']),old_scores=r['analytic_old_scores_upper'],
                **normal_limits(inputs[r['case_id']],r))
    (R/'continuous_summary.json').write_text(json.dumps(bands,indent=2)+'\n')
    for group in dict.fromkeys(x['cohort'] for x in out):
        a=[x for x in out if x['cohort']==group]
        print(group,len(a),[sum(x[k] for x in a) for k in ('old_contains','first_contains','both_contains')],
              'stop',sum(x['first_upper']<=.5 for x in a),'full',sum(x['both_upper']<=.5 for x in a))
    for k,v in bands.items():
        print(k,round(v['conditional_binary_error_lower'],5),round(v['required_energy_for_disjoint_balls_lower'],3))
if __name__=='__main__':builds()
