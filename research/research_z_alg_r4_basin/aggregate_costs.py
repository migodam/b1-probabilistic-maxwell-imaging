"""Single-pass accounting from completed manifests, without replaying physics."""
import csv,json
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def main():
    rows=[]
    def add(ident,kind,source,counts,wall=None,cpu=None,**extra):
        counts=counts or {}
        rows.append(dict(id=ident,kind=kind,source=source,full_LU=counts.get('LU',0),
           exposed_complex_LU_RHS=counts.get('actual_complex_RHS',0),source_RHS=counts.get('source_RHS',0),
           adjoint_RHS=counts.get('adjoint_RHS',0),receiver_transpose_RHS=counts.get('receiver_transpose_RHS',0),
           exact_HVP_tangent_RHS=counts.get('tangent_RHS',0),exact_HVP_adjoint_RHS=counts.get('incremental_adjoint_RHS',0),
           dense_ray_kernel_RHS=counts.get('ray_kernel_RHS',0),dense_Schur_576=counts.get('dense_complex_schur',0),
           transfer_triangular_RHS_576=counts.get('triangular_small_transfer_RHS',0),
           ray_response_evaluations=counts.get('ray_transfer_evaluations',0),FullJ=counts.get('FullJ',0),
           Gram=counts.get('Gram',0),Gram_RHS_12=counts.get('Gram_RHS',0),
           HVP=counts.get('HVP',0),real_J_SVD=counts.get('diagnostic_SVD',0),
           wall_seconds=wall,cpu_seconds=cpu,**extra))
    for path in sorted((ROOT/'results').glob('*/MANIFEST.json')):
        d=json.loads(path.read_text())
        for row in d.get('rows',[]):
            add(row['output'],'solver',row['output'],row.get('costcounts') if isinstance(row.get('costcounts'),dict) else {},
                row.get('wall_seconds'),row.get('cpu_seconds'),run=d['run_name'],case=row['case_id'],method=row['method'],
                status=row['status'],scope='FRESH' if d['run_name']=='fresh_v3' else 'DEVELOPMENT',
                note='All completed and failed paid work included; serialization/process startup outside solver timer.')
    p=ROOT/'diagnostics/SUMMARY.json';d=json.loads(p.read_text())['costs']
    add('geometry-paths','diagnostic',str(p.relative_to(ROOT)),d['counts'],d['wall_seconds'],d['cpu_seconds'],
        status='COMPLETED',note='LAPACK zgecon internal RHS are opaque, not zero; no historical replay cost added.')
    p=ROOT/'diagnostics/operator_certificate.json';d=json.loads(p.read_text())
    add('operator-certificate','diagnostic',str(p.relative_to(ROOT)),d['counts'],d['wall_seconds'],d['cpu_seconds'],status=d['status'])
    p=ROOT/'diagnostics/hessian_geometry.json';d=json.loads(p.read_text())
    for case in d['cases']:
        cost=case['costs'];add('hessian-'+case['case'],'diagnostic',str(p.relative_to(ROOT)),cost['counts'],cost['wall_seconds'],cost['cpu_seconds'],
            status=case['eigensolver'].get('status','RECORDED'),note='Includes ARPACK attempts that exhausted the declared HVP cap.')
        follow=case.get('dense_freeface_followup',{})
        if 'costs' in follow:
            cost=follow['costs'];add('hessian-dense-'+case['case'],'diagnostic',str(p.relative_to(ROOT)),cost['counts'],cost['wall_seconds'],cost['cpu_seconds'],status=follow['status'])
    for split in ('FACTOR','FRESH'):
        p=ROOT/'audit'/('GENERATION_'+split+'.json');d=json.loads(p.read_text())
        add('generation-'+split.lower(),'generation',str(p.relative_to(ROOT)),d['model_counts'],d['wall_seconds'],d['cpu_seconds'],status='COMPLETED',note='Reused preexisting public Green/operators; no Green regeneration.')
    add('failed-font-generation','failure','audit/GENERATOR_FONT_FAILURE.json',{},.591291458,None,status='FAILED_BEFORE_FORWARD',note='No GT saved and no physical solve; process CPU not captured.')
    for path in sorted((ROOT/'results').glob('*/SCORES.json')):
        add('score-'+path.parent.name,'scoring',str(path.relative_to(ROOT)),{},None,None,status='COMPLETED',note='Array scoring and plotting have zero new Maxwell solves; separate timers not captured.')
    fields=sorted(set().union(*(r.keys() for r in rows)))
    with (ROOT/'COST_LEDGER.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    numeric=['full_LU','exposed_complex_LU_RHS','dense_Schur_576','transfer_triangular_RHS_576','FullJ','Gram','HVP','real_J_SVD']
    totals={k:sum(r[k] for r in rows) for k in numeric}
    totals.update(solver_wall_seconds=sum(r['wall_seconds'] or 0 for r in rows if r['kind']=='solver'),
       solver_cpu_seconds=sum(r['cpu_seconds'] or 0 for r in rows if r['kind']=='solver'),
       known_all_wall_seconds_sum=sum(r['wall_seconds'] or 0 for r in rows),known_all_cpu_seconds_sum=sum(r['cpu_seconds'] or 0 for r in rows),
       solver_outputs=sum(r['kind']=='solver' for r in rows),ledger_rows=len(rows))
    methods={}
    for r in rows:
        if r.get('run')=='fresh_v3':
            m=methods.setdefault(r['method'],dict(LU=0,complex_LU_RHS=0,triangular_transfer_RHS=0,Schur=0,wall_seconds=0.,cpu_seconds=0.))
            for key,col in [('LU','full_LU'),('complex_LU_RHS','exposed_complex_LU_RHS'),('triangular_transfer_RHS','transfer_triangular_RHS_576'),('Schur','dense_Schur_576'),('wall_seconds','wall_seconds'),('cpu_seconds','cpu_seconds')]:m[key]+=r[col] or 0
    result=dict(totals=totals,fresh_method_totals=methods,
      limitations=['Numerical solver timers include physical and optimizer work, not per-job JSON/NPZ serialization or process startup.',
       'Triangular transfer RHS are full dimension576 and are charged separately, not called free reduced solves.',
       'The max_rhs budget bounds exposed physical-LU/kernel RHS; separate transfer triangular RHS are limited by wall, not that counter.',
       'LAPACK estimator and SciPy internal factorization work are opaque but their inclusive wall/CPU is captured.',
       'Old R3 costs are provenance references, not charged again as newly run R4 work.',
       'Concurrent worker times must not be summed as sprint elapsed time.'])
    (ROOT/'COST_TOTALS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
