"""Remaining discrepancy allowed after all four fixed reference/blank choices."""
import csv
import json
from pathlib import Path

from flint import arb, ctx
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'core'))
from certificate import safe_float


def summarize(input_name, radius_key, output_table, output_summary):
    ctx.prec=128
    anchors=json.loads((HERE/'four_calibrations.json').read_text())['rows']
    inp=json.loads((HERE.parent/'manitoba'/input_name).read_text())['rows']
    data={(x['target'],x['frequency_GHz'],x['cap'],x['arm']):x for x in inp}
    rows=[]
    for row in anchors:
        for cap in (.25,.5,.75,1.):
            for model,arms in [('passive_sector',('real','positive_optical')),
                               ('exactly_lossless',('real','positive_optical','free_optical'))]:
                best=max((data[row['target'],row['frequency_GHz'],cap,a] for a in arms),
                         key=lambda x:x['verified_lower'])
                d=arb(best['verified_lower'])
                radius=arb(row[radius_key])
                rem=d-radius
                norm_lo,norm_hi=row['primary_norm']
                # For positive residual allowance the upper denominator gives a
                # lower relative allowance. A negative value is a diagnostic.
                relative=rem/(arb(norm_hi) if rem >= 0 else arb(norm_lo))
                rows.append(dict(target=row['target'],GHz=row['frequency_GHz'],cap=cap,model=model,
                                 best_eligible_arm=best['arm'],lower_distance=best['verified_lower'],
                                 four_calibration_radius_upper=row[radius_key],
                                 margin_lower=safe_float(rem,upper=False),relative_margin_lower=safe_float(relative,upper=False),
                                 positive=bool(rem>0)))
    with (HERE/output_table).open('w',newline='') as file:
        writer=csv.DictWriter(file,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    groups={}
    for model in ('passive_sector','exactly_lossless'):
        groups[model]={}
        for cap in (.25,.5,.75,1.):
            selected=[r for r in rows if r['model']==model and r['cap']==cap and r['positive']]
            groups[model][str(cap)]=dict(count=len(selected),objects=len({r['target'] for r in selected}),
                       positives=[dict(target=r['target'],GHz=r['GHz'],margin=r['relative_margin_lower']) for r in selected])
    (HERE/output_summary).write_text(json.dumps(groups,indent=2)+'\n')
    for k,v in groups.items():print(k,[(c,x['count'],x['objects']) for c,x in v.items()])

if __name__=='__main__':
    summarize('verified_scan.json','hull_radius_upper','transfer_margins.csv','summary.json')
    summarize('reciprocal_verified.json','symmetric_hull_radius_upper','symmetric_margins.csv','symmetric_summary.json')
