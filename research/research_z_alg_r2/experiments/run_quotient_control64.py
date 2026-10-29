"""Frozen64-start reduced Full-J coverage control on two exposed root cases."""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_key]='1'
import sys,json,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from scipy.stats import qmc
P=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(P/'src'));sys.path.insert(0,str(P/'experiments'))
from physics import Model,ReducedKernel
from validate import CountingBook
from run_quotient_benchmark import data_rows,fullj,compare,qualified
from run_branch_study import write
IDS=['adaptive_E1_LOW6_02_slice_12','high4_E2_HIGH4_01']

def control_starts(row):
    chosen=[np.asarray(x).copy() for x in row['starts']]
    lower,upper=row['box']
    for point in qmc.Sobol(4,scramble=False).random_base2(8):
        x=lower+(upper-lower)*point
        if all(np.linalg.norm(x-old)>1e-12 for old in chosen):chosen.append(x)
        if len(chosen)==64:break
    return np.asarray(chosen)

def main():
    path=P/'experiments/QUOTIENT_CONTROL64_CONTRACT.json'
    rows,adaptive=data_rows();rows=[r for r in rows if r['id'] in IDS]
    if '--freeze' in sys.argv:
        if path.exists():raise RuntimeError('Control64 contract already frozen')
        write(path,dict(schema='Z-ALG-R2-REDUCED-FULLJ64-CONTROL-v1',frozen_utc=datetime.now(timezone.utc).isoformat(),
            exposure='EXPOSED_DATA_DIAGNOSTIC_NOT_FRESH_BLIND',tasks=IDS,method='ReducedFullJ64',
            starts={row['id']:control_starts(row) for row in rows},
            initialization='Exactoriginal8 starts thenfirst56 additionaldistinctunscrambledSobol4 points',
            max_nfev_per_start=120,ftol=1e-12,xtol=1e-12,gtol=1e-12,bounds='sameoriginalbox',
            kernel='exactsameglobalmidpointkernel plusSchur restriction;analyticreducedJacobian',
            validation='freshfullnonlinear solve eachunique legal candidate,response1e-8,state1e-10,sigma1e-6',
            root_merge_RMS=1e-7,truth='notread orused',comparison='rootsets versusquotient;nottruthscores'))
        print('CONTROL64_FROZEN',flush=True);return
    contract=json.loads(path.read_text())
    if list((P/'results').glob('quotient_control64_*.json')):raise RuntimeError('Preserve firstcontrol64batch')
    outputs=[];start=time.perf_counter()
    for row in rows:
        row['starts']=np.asarray(contract['starts'][row['id']])
        book=CountingBook();model=Model(book,**row['geometry'])
        kernel=ReducedKernel(model,np.full(model.n,np.mean(row['box'])),active=np.arange(model.n))
        if model.n==6:kernel=kernel.restrict([0,1,2,3],row['fixed'])
        out=fullj(model,row,kernel);out['method']='ReducedFullJ64'
        old=json.loads((P/'results'/f'quotient_bench_coverage_{row["id"]}_quotient.json').read_text())
        out.update(accounting=book.summary(model.ns,model.nr),model_counts=dict(model.counts),
                   vs_quotient=compare(old,out),exposure=contract['exposure'])
        outputs.append(out);write(P/'results'/f'quotient_control64_{row["id"]}.json',out)
        print('CONTROL64',row['id'],'roots',len(qualified(out)),out['vs_quotient'],flush=True)
    anchor=adaptive['anchors']['E1_LOW6_02']['original_cost']
    summary=dict(schema='Z-ALG-R2-REDUCED64-SUMMARY-v1',cases=[dict(id=o['id'],roots=len(qualified(o)),vs_quotient=o['vs_quotient']) for o in outputs],
        all_quotient_roots_recovered=all(o['vs_quotient']['same_root_set'] for o in outputs),
        new_actual_RHS=sum(o['accounting']['total_complex_rhs'] for o in outputs),
        reused_anchor_cost=anchor,wall_including_logging=time.perf_counter()-start,
        scope='64starts is afinitecoveragecontrol;noexhaustivenesscertificate orGT claim')
    write(P/'results/quotient_control64_summary.json',summary)
    print('CONTROL64_COMPLETE',json.dumps(summary),flush=True)

if __name__=='__main__':main()
