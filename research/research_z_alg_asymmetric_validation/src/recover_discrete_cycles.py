"""Complete interrupted cycle reporting; never repeat the 33 paid multistarts."""
import json
import numpy as np
from accounting import P,Book,write
from physics import model_for,quality
from charts import MaterialChart
from run_discrete import unpack

def main():
    out=P/'results/E2';old=json.loads((out/'INTERRUPTED.json').read_text());contract=json.loads((P/'EXPERIMENT_CONTRACT.json').read_text())['E2']
    assert len(old['rows'])==33 and 'multiple values' in old['error']
    journal=[json.loads(x) for x in (P/'COST_LEDGER.jsonl').read_text().splitlines()]
    oldrun=[r for r in journal if r['run']=='E2_ASYMMETRIC_DISCRETE_BRANCHES']
    begins={r['action']:r for r in oldrun if r['event']=='BEGIN'}
    saves=[r for r in oldrun if r['event']=='SAVE']
    last_state=next(r for r in reversed(saves) if begins[r['action']]['kind']=='FULL_FORWARD')
    last_J=next(r for r in reversed(saves) if begins[r['action']]['kind']=='FULL_MATERIAL_J_ADJOINT')
    assert begins[last_state['action']]['label']=='CYCLE_FULLJ'
    with np.load(P/last_state['path']) as a:cached={k:a[k].copy() for k in a.files}
    with np.load(P/last_J['path']) as a:cachedJ=a['Jreal'].copy();assert np.array_equal(a['m'],cached['m'])
    with np.load(out/'START_000.npz') as a:m0=a['m'].copy();target=a['target'].copy()
    b=Book('E2_CYCLE_REPORTING_RECOVERY');model=model_for(b,'E2');base=model.pack(target)
    rng=np.random.default_rng(contract['seed']);rng.uniform(.5,4,size=(contract['multistarts'],4));directions=rng.normal(size=(2,len(base)));directions/=np.linalg.norm(directions,axis=1)[:,None]
    dep=np.arange(4);free=np.array([],dtype=int);cycles=[]
    try:
        for k,d in enumerate(directions):
            m=m0.copy();nodes=[];complete=True
            for step in range(1,33):
                t=step/32;wanted=base+.02*np.linalg.norm(base)*d*np.sin(2*np.pi*t);yy=unpack(wanted,model)
                before=b.run_used
                if k==0 and step==1:
                    state=dict(m=cached['m'],j=cached['j'],E=cached['E'],y=cached['y'],J=cachedJ,
                        state_relative_residual=float(cached['state_relative_residual']),resolvent_sigma_min=float(cached['resolvent_singular_values'][-1]))
                    m=state['m'];reused=True
                else:
                    chart=MaterialChart(model,yy,m0,dep,free);r=chart.correct(m,np.array([]),'RECOVERED_CYCLE_'+str(k)+'_'+str(step),maxiter=16)
                    state=model.forward(r['m'],jac=True,label='CYCLE_FULLJ');m=state['m'];reused=False
                q=quality(model,state,yy,m0);sv=np.linalg.svd(state['J'],compute_uv=False)
                node=dict(step=step,t=t,material=m.tolist(),response_sigma_min=float(sv[-1]),new_RHS=b.run_used-before,reused_paid_first_node=reused,**q)
                nodes.append(node)
                np.savez_compressed(out/('CYCLE_'+str(k)+'_'+str(step).zfill(2)+'.npz'),m=m,y=state['y'],j=state['j'],J=state['J'],target=yy)
                if not q['qualified']:complete=False;break
            cycles.append(dict(cycle=k,completed=complete,nodes=nodes,return_RMS=float(np.linalg.norm(m-m0)/2),interpretation='Retraced data path; no nontrivial monodromy claim'))
        result=dict(status='COMPLETE_AFTER_REPORTING_FIX',roots=old['roots'],root_count_found=len(old['roots']),starts=old['rows'],cycles=cycles,
            root_count_qualification='Found roots only, no bounded-domain uniqueness or full-branch certificate',
            algebraic_homotopy='NOT_RUN: finite-cover/discriminant preconditions not established',permutation_audit=model.permutation_audit(),
            paid_first_run_RHS=old['paid'],paid_recovery_RHS=b.run_used,total_paid_RHS=old['paid']+b.run_used,
            correction='Rename response_sigma_min to avoid duplicate state-resolvent sigma_min JSON keyword; no algorithm/threshold change; paid first node reused')
        write(out/'SUMMARY.json',result);b.end();print(json.dumps({k:v for k,v in result.items() if k not in ('starts','cycles')},indent=2))
    except BaseException as exc:b.end('FAILED',error=repr(exc));raise
if __name__=='__main__':main()
