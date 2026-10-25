"""Complete all18 logical observation attempts; adaptive original-coordinate diagnostics."""
from run_v3_strata import *
D=OUT/'remaining_spaces';D.mkdir(exist_ok=True)
def main():
    t=Target('A','A1_remaining_v1',output_dir=OUT/'physics/A1_remaining_v1');rows=[]
    for si,pair in [(0,(0,1)),(2,(12,27))]:
        model=loadmodel(t,pair);grid=[]
        if si==2:
            for a in np.linspace(-3.999999,3.999999,9):
                for b in np.linspace(-3.999999,3.999999,9):grid.append(model.full(np.array([a,b]),'remaining_seed_grid'))
        cases=[(1,2)] if si==0 else [(p,n) for p in range(2) for n in range(3)]
        for pi,ni in cases:
            name=f'S{si}_P{pi}_N{ni}'
            if (D/f'{name}_RESULT.json').exists():rows.append(json.loads((D/f'{name}_RESULT.json').read_text()));continue
            with np.load(OUT/f'v2/datasets/{name}.npz') as z:d=z['d'];sig=float(z['sigma'])
            before=t.book.totals();old=OUT/f'v1_failed_preserved/physical/{name}_MODES.json'
            if old.exists():
                data=json.loads(old.read_text());m=data['modes'][0];m['z']=np.array(m['z']);m['C']=np.array(m['C']);m['r']=model.full(m['z'],'remaining_paid_mode');m['J']=model.jac(m['z'],'remaining_mode_J')/sig;modes=[m]
            else:modes,attempts=discover(model,grid,d,sig);m=modes[0]
            results={};raw=None
            for n in ([15,21] if si==0 else [9,13]):
                ss,raw=original(model,modes,d,sig,n);results[str(n)]=ss;report(D/f'{name}_ORIGINAL_{n}.json',ss)
            ref=list(results.values())[-1];conv=difference(ref,list(results.values())[-2]);metrics,status,H=baseline_metrics(model,m,d,sig,ref,raw)
            rt=rto(Conditional(model,m,d,sig),8 if si==0 else 4);report(D/f'{name}_QR_RTO.json',rt)
            # Declared saturation-exposure probe: measured from reference-node decoder only.
            x0=model.V@m['z'];_,aux=t.decoder(x0);frac=aux['fraction'][t.prior.mask]
            row=dict(observation=name,space=si,parent=pi,noise_multiplier=[1,.5,.25][ni],original_convergence=conv,baseline_metrics=metrics,baseline_status=status,RTO_successes=rt['succeeded'],RTO_status=rt['status'],mode_fraction_host=float(np.mean(frac==0)),mode_fraction_pure=float(np.mean(frac==1)),mode_fraction_transition=float(np.mean((frac>0)&(frac<1))),physical_response_chart='NOT_RUN: remaining fixed stage budget prioritized all18 original-coordinate attempts and coupled mixed-rank study',qualification='NUMERICAL_REFERENCE_UNQUALIFIED_UNLESS_CONVERGED',scale_up='HOLD',cost_before=before,cost=t.book.totals())
            rows.append(row);report(D/f'{name}_RESULT.json',row);report(D/'RUN_STATUS.json',dict(status='RUNNING',results=rows,cost=t.book.totals()));print(json.dumps(plain(row)),flush=True);model.close()
    report(D/'RUN_STATUS.json',dict(status='CLOSED_NUMERICAL_NOT_CERTIFIED',results=rows,cost=t.book.totals()))
if __name__=='__main__':
    try:main()
    except BaseException as e:report(D/'RUN_FAILURE.json',dict(error=repr(e),traceback=traceback.format_exc()));raise
