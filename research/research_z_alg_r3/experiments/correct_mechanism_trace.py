"""Clarify archived representatives using paid V6 traces, with zero new RHS."""
import datetime
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def main():
    source=ROOT/'results'/'holdout_globaldual_mechanism.json'
    target=ROOT/'results'/'holdout_globaldual_mechanism_v2.json'
    if target.exists():raise FileExistsError(target)
    result=json.loads(source.read_text())
    scores=json.loads((ROOT/'results'/'holdout_scores.json').read_text())
    scoremap={r['file']:r for r in scores['rows']}
    for row in result['rows']:
        prediction=json.loads((ROOT/row['file']).read_text())
        scored=scoremap[row['file']]
        tr={x['seconds']:x for x in scored['online_anytime_top1']}
        for endpoint in row['raw_endpoint_rows']:
            correction=endpoint.get('correction')
            if correction is None:continue
            correction['archived_representative_relative_residual']=correction.pop('relative_residual')
            correction['archived_representative_material_RMS']=correction.pop('material_RMS')
            correction['archive_boundary']='A near-duplicate fit endpoint can merge with an earlier archive entry. Actual global improvements below come from the exact paid-trial trace.'
            origin='globaldual_fullj_correction_'+str(endpoint['start'])
            improvements=[t for t in prediction['best_so_far'] if t['origin']==origin]
            if improvements:
                paid=min(improvements,key=lambda t:t['relative_residual'])
                correction['best_paid_global_improvement']=dict(relative_residual=paid['relative_residual'],seconds=paid['seconds'],material_RMS=tr[paid['seconds']]['top1_RMS'])
        if row['pure_endpoint_original_residual_top1'] is not None:
            start=row['pure_endpoint_original_residual_top1']['start']
            row['pure_endpoint_original_residual_top1']=next(r for r in row['raw_endpoint_rows'] if r['start']==start)
    result['schema']='R3_GLOBALDUAL_MECHANISM_V2_ARCHIVE_BOUNDARY'
    result['derived_from']=str(source.relative_to(ROOT))
    result['derived_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
    result['new_physics_rhs']=0
    result['scoring_only']=False
    result['cost_boundary']='Derived metadata clarification only, reuses all prior independent validation and offline scoring; original v1 costs remain counted once'
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(output=str(target),new_physics_rhs=0,private_truth_read=False)))


if __name__=='__main__':main()
