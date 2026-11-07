"""Public receipt summary; does not open truth or numerical arrays."""
from pathlib import Path
from collections import Counter
import json
import csv

S=Path(__file__).resolve().parents[1]

def main():
    records=json.loads((S/'data/manifests/online.json').read_text())['records']
    meta={r['obs_id']:r for r in records}
    groups={}
    stops=[]
    for f in sorted((S/'anchors').glob('obs*.json')):
        d=json.loads(f.read_text())
        obs=f.stem.split('_weakreg')[0]
        r=meta[obs]
        host='i9 Windows' if str(d['public_path']).startswith('D:') else 'M4 macOS'
        key=('weakreg45' if '_weakreg'in f.stem else 'dataGN',r['split'],host)
        x=groups.setdefault(key,dict(n=0,task_seconds=0.,LU=0,complex_RHS=0,
                                      statuses=Counter(),stop_reasons=Counter(),final_stages=Counter()))
        x['n']+=1
        x['task_seconds']+=d.get('anchor_plus_feature_wall_seconds',
            d['solve_wall_seconds']+d.get('feature_wall_seconds',0.))
        x['LU']+=d['counts'].get('LU',0)
        x['complex_RHS']+=d['counts'].get('actual_complex_RHS',0)
        x['statuses'][d['status']]+=1
        x['stop_reasons'][str(d['stop_reason'])]+=1
        x['final_stages'][str(d['final_stage'])]+=1
        if d['status']=='NUMERICAL_STOP':
            stops.append(dict(obs_id=obs,split=r['split'],host=host,
                              reason=d['stop_reason'],LU=d['counts']['LU'],
                              complex_RHS=d['counts']['actual_complex_RHS'],
                              accepted_material_resolution=d['final_stage']))
    rows=[dict(method=k[0],split=k[1],host=k[2],**v)for k,v in groups.items()]
    result=dict(scope='Completed endpoint receipts; initial failed attempt is extra',
                rows=rows,numerical_stops=stops,
                note='Concurrent task seconds are not overall elapsed. Numerical termination is not a GT recovery certificate.',
                additional_failure='Initial SVD failure of obs_8b81ff10f9d25c2c896b has 30.3716685s and unknown LU/RHS; do not report an exact all-attempt total.')
    (S/'audit/ANCHOR_COMPLETED_COUNTS.json').write_text(json.dumps(result,indent=2))
    with (S/'audit/ANCHOR_COST_BY_SPLIT_PLATFORM.csv').open('w',newline='')as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader()
        for r in rows:
            w.writerow({k:json.dumps(v)if isinstance(v,dict)else v for k,v in r.items()})
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
