"""One independent total-field solve for each declared successful endpoint.

Same disk-pulse discretization, different state formulation/CPU implementation.
This is not a cross-grid or continuous-Maxwell verification, or an error-floor
test. The plan is finite; no trial-point or derivative checks are repeated.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.linalg import solve

def validate(plan):
    rows=[];beg=time.perf_counter();rhs=0
    for item in plan['cases']:
        with np.load(item['operators_path']) as z:G,E,S=z['G'],z['E0'],z['S']
        with np.load(item['observation_path']) as z:Y=z['Y']
        pred=json.loads(Path(item['prediction_path']).read_text(encoding='utf-8-sig'))
        m=np.asarray(pred.get('result',pred)['m'],float)
        # Original total-field equation (I-GX)Q=E, not the CUDA current LU.
        Q=solve(np.eye(len(m))-G*m[None,:],E,check_finite=False)
        fitted=S@(m[:,None]*Q);relative=float(np.linalg.norm(fitted-Y)/np.linalg.norm(Y))
        rhs+=E.shape[1]
        rows.append({'case_id':item['case_id'],'relative_residual':relative,'material_box_ok':bool(np.all((m>=.02)&(m<=20.))),
                     'data_limit':item.get('data_limit',.001),'data_pass':relative<=item.get('data_limit',.001)})
    return {'scope':'one CPU original total-field solve per declared endpoint, same discrete physicalmodel',
            'rows':rows,'costs':{'LU':len(rows),'source_RHS':rhs},'wall_s':time.perf_counter()-beg}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('plan');ap.add_argument('output');a=ap.parse_args()
    result=validate(json.loads(Path(a.plan).read_text()))
    Path(a.output).write_text(json.dumps(result,indent=2));print(json.dumps({'validated':len(result['rows']),'pass':sum(x['data_pass'] for x in result['rows'])}))
if __name__=='__main__':main()
