import json,time
import numpy as np
from accounting import P,Book,write
from physics import model_for,quality
from charts import MaterialChart

def unpack(v,model):
    a=np.asarray(v).reshape(model.ns,model.nr,2)
    return (a[:,:,0]+1j*a[:,:,1]).T

def main():
    c=json.loads((P/'EXPERIMENT_CONTRACT.json').read_text())['E2']
    if not all((P/'results/E1'/('P'+str(i))/m/'SUMMARY.json').exists() for i in range(3) for m in ('A','B')):raise RuntimeError('Complete E1 first')
    b=Book('E2_ASYMMETRIC_DISCRETE_BRANCHES');model=model_for(b,'E2');out=P/'results/E2';out.mkdir(parents=True,exist_ok=False)
    m0=np.array(c['public_reference_material']);ref=model.forward(m0,jac=True,label='public_reference4');y=ref['y']
    dep=np.arange(model.n);free=np.array([],dtype=int);chart=MaterialChart(model,y,m0,dep,free)
    rng=np.random.default_rng(c['seed']);starts=[m0.copy(),*rng.uniform(.5,4,size=(c['multistarts'],4))]
    rows=[];roots=[];minimum_singular=[]
    try:
        for i,initial in enumerate(starts):
            before=b.run_used;t=time.perf_counter();r=chart.correct(initial,np.array([]),'MULTISTART_'+str(i),maxiter=c['maximum_iterations_per_start'])
            state=model.forward(r['m'],jac=True,label='ENDPOINT_FULLJ_'+str(i));J=state['J'];sv=np.linalg.svd(J,compute_uv=False)
            verified=model.forward(r['m'],independent=True,label='ENDPOINT_INDEPENDENT_'+str(i));q=quality(model,verified,y,m0)
            root_id=None
            if q['qualified']:
                for k,old in enumerate(roots):
                    if np.linalg.norm(old-r['m'])/2<=c['material_root_merge_RMS']:root_id=k;break
                if root_id is None:root_id=len(roots);roots.append(r['m'].copy())
            row=dict(start=i,initial=initial.tolist(),material=r['m'].tolist(),iterations=r.get('iterations'),corrector_success=r['success'],reason=r.get('reason','CONVERGED'),
                numerical_rank=int(np.sum(sv>1e-10*sv[0])),singular_values=sv.tolist(),root_id=root_id,**q,actual_RHS=b.run_used-before,wall_seconds=time.perf_counter()-t)
            row['classification']='NUMERICALLY_SAME_FULL_LABELED_RESPONSE' if q['qualified'] else ('FINITE_NOISE_COMPETITIVE_ONLY' if q['pair_whitened']<=1 else 'NOT_A_RESPONSE_MATCH')
            rows.append(row);np.savez_compressed(out/('START_'+str(i).zfill(3)+'.npz'),initial=initial,m=r['m'],y=verified['y'],j=verified['j'],E=verified['E'],J=J,target=y)
            write(out/('START_'+str(i).zfill(3)+'.json'),dict(result=row,history=r['history']))
            write(out/'CHECKPOINT.json',dict(completed=len(rows),roots=[x.tolist() for x in roots],paid_RHS=b.run_used))
        cycles=[];base=model.pack(y);directions=rng.normal(size=(2,len(base)));directions/=np.linalg.norm(directions,axis=1)[:,None]
        for k,direction in enumerate(directions):
            m=m0.copy();nodes=[];complete=True
            for step in range(1,33):
                t=step/32;wanted=base+.02*np.linalg.norm(base)*direction*np.sin(2*np.pi*t)
                yy=unpack(wanted,model);cy=MaterialChart(model,yy,m0,dep,free)
                before=b.run_used;r=cy.correct(m,np.array([]),'CYCLE_'+str(k)+'_'+str(step),maxiter=16)
                verified=model.forward(r['m'],jac=True,label='CYCLE_FULLJ');q=quality(model,verified,yy,m0);sv=np.linalg.svd(verified['J'],compute_uv=False)
                node=dict(step=step,t=t,qualified=q['qualified'],material=r['m'].tolist(),sigma_min=float(sv[-1]),RHS=b.run_used-before,**{a:v for a,v in q.items() if a!='qualified'})
                nodes.append(node)
                np.savez_compressed(out/('CYCLE_'+str(k)+'_'+str(step).zfill(2)+'.npz'),m=r['m'],y=verified['y'],j=verified['j'],J=verified['J'],target=yy)
                if not q['qualified']:complete=False;break
                m=r['m']
            cycles.append(dict(cycle=k,completed=complete,nodes=nodes,return_RMS=float(np.linalg.norm(m-m0)/2),
                interpretation='Retraced data line; identity expected under unique regular bounded lifting, not nontrivial monodromy'))
        result=dict(status='COMPLETE',roots=[x.tolist() for x in roots],root_count_found=len(roots),root_count_qualification='Lower bound only; no bounded-domain exclusion or completeness certificate',
            starts=rows,cycles=cycles,permutation_audit=model.permutation_audit(),
            algebraic_homotopy='NOT_RUN: no finite-cover root count or discriminant control established',source_labels='fixed, never permuted',paid_combined_RHS=b.run_used)
        write(out/'SUMMARY.json',result);b.end();print(json.dumps({k:v for k,v in result.items() if k not in ('starts','cycles')},indent=2))
    except BaseException as exc:
        write(out/'INTERRUPTED.json',dict(error=repr(exc),rows=rows,roots=[x.tolist() for x in roots],paid=b.run_used));b.end('FAILED_OR_BUDGET_STOP',error=repr(exc));raise
if __name__=='__main__':main()
