"""Qualified-density Gaussian baseline metrics on regular mixed-stratum references."""
from run_v3_strata import *
def main():
    t=Target('A','A1_baseline_final',output_dir=OUT/'physics/A1_baseline_final');model=loadmodel(t,(0,64));g=grid(model,48);rows=[]
    for ni in range(3):
        name=f'S1_P1_N{ni}'
        with np.load(OUT/f'v2/datasets/{name}.npz') as z:d=z['d'];sig=float(z['sigma'])
        dat=json.loads((S/f'{name}_MODES.json').read_text());m=dict(z=np.array(dat['z']),C=np.array(dat['C']));m['r']=model.full(m['z'],'final_baseline_paid_mode');m['J']=model.jac(m['z'],'final_baseline_paid_J')/sig
        ref,raw=summary_parts(model,m,d,sig,21,g);bm,bs,H=baseline_metrics(model,m,d,sig,ref,raw)
        for k,v in bm.items():v['density_metric_qualification']='NUMERICAL_REFERENCE_CONVERGED; finite integration error remains, no rigorous global TV certificate'
        row=dict(observation=name,baselines=bm,status=bs,cost=t.book.totals());rows.append(row);report(S/f'{name}_GAUSSIAN_BASELINES.json',row);print(json.dumps(plain(row)),flush=True)
    report(S/'GAUSSIAN_BASELINES_SUMMARY.json',rows)
if __name__=='__main__':main()
