import argparse,json,time
from pathlib import Path
import numpy as np
from accounting import Book,P,write,BudgetStop
from physics import model_for,quality
from charts import MaterialChart,CurrentChart,partition

def run(parent,method):
    contract=json.loads((P/'EXPERIMENT_CONTRACT.json').read_text());c=contract['E1']
    if json.loads((P/'audit/INITIAL_PHYSICAL_PROBE.json').read_text())['status']!='PASS':raise RuntimeError('Probe gate not passed')
    book=Book('E1_P'+str(parent)+'_'+method);out=P/'results/E1'/('P'+str(parent))/method;out.mkdir(parents=True,exist_ok=False)
    model=model_for(book);mref=np.array(c['public_reference_materials'][parent]);reference=model.forward(mref,jac=True,label='COMMON_REFERENCE')
    dep,free=partition(reference['J']);target=reference['y'];common_paid=book.run_used
    write(out/'CHART_IDENTITY.json',dict(parent=parent,method=method,dep=dep,free=free,shared_reference_RHS=common_paid,target_role='public physical existence reference',new_prior='P_Z_real_product_linear'))
    setup_tick=time.perf_counter()
    chart=MaterialChart(model,target,mref,dep,free) if method=='A' else CurrentChart(model,target,reference,dep,free)
    zref=None if method=='A' else chart.encode(reference['j'],mref)
    setup_wall=time.perf_counter()-setup_tick
    np.savez_compressed(out/'REFERENCE.npz',m=mref,j=reference['j'],E=reference['E'],y=target,J=reference['J'],G=model.G,S=model.S,E0=model.E0,dep=dep,free=free)
    rows=[];accepted=[dict(free=mref[free].copy(),m=mref.copy(),z=zref,j=reference['j'])]

    def one(values,key,kind,seed=None):
        if seed is None:seed=min(accepted,key=lambda a:np.linalg.norm(a['free']-values))
        before=book.run_used;tick=time.perf_counter()
        result=chart.correct(seed['m'] if method=='A' else seed['z'],values,key)
        row=dict(parent=parent,method=method,key=key,kind=kind,free_values=values.tolist(),corrector_success=result['success'],iterations=result.get('iterations',0),reason=result.get('reason','CONVERGED'))
        folder=out/key;folder.mkdir(parents=True,exist_ok=False)
        write(folder/'ITERATIONS.json',result['history'])
        legal=result.get('m') is not None and np.all((result['m']>=.5)&(result['m']<=4))
        entry=None
        if legal:
            # New independent full forward, never fitted current as its own truth.
            verified=model.forward(result['m'],label='INDEPENDENT_ENDPOINT_'+key,independent=True)
            q=quality(model,verified,target,mref);row.update(q)
            arrays=dict(m=result['m'],j_full=verified['j'],E_full=verified['E'],y_full=verified['y'],target=target,free_values=values)
            if method=='B':
                r=np.linalg.norm(result['j']-result['m'][:,None]*result['E'])/max(np.linalg.norm(result['j']),1e-300)
                denom=np.sum(abs(result['E'])**2,axis=1)
                validnodes=denom>1e-20
                recovered=np.full(model.n,np.nan+1j*np.nan,dtype=np.complex128)
                recovered[validnodes]=np.sum(result['E'][validnodes].conj()*result['j'][validnodes],axis=1)/denom[validnodes]
                row.update(shared_material_relative=float(r),receiver_null_relative=float(np.linalg.norm(model.S@result['j']-target)/np.linalg.norm(target)),minimum_field_squared=float(np.min(denom)),node_cells=np.flatnonzero(~validnodes).tolist(),recovered_real_material_error=float(np.max(abs(recovered[validnodes]-result['m'][validnodes]))) if np.any(validnodes) else None)
                # Explicit material variables persist at nodes. Do not divide by zero.
                row['qualified']=bool(row['qualified'] and r<=1e-10 and row['receiver_null_relative']<=1e-10)
                arrays.update(z=result['z'],j_chart=result['j'],E_chart=result['E'],recovered_material=recovered)
            np.savez_compressed(folder/'ENDPOINT.npz',**arrays)
            if row['qualified']:
                entry=dict(free=values.copy(),m=result['m'].copy(),z=result.get('z'),j=verified['j'])
                accepted.append(entry)
        else:row.update(qualified=False,legal=False,material_RMS=None,pair_relative=None)
        row.update(actual_RHS=book.run_used-before,controlled_wall_seconds=time.perf_counter()-tick)
        write(folder/'RESULT.json',row);rows.append(row)
        write(out/'CHECKPOINT.json',dict(completed=len(rows),qualified=sum(r['qualified'] for r in rows),paid_RHS=book.run_used,latest=key))
        return entry,row

    try:
        offsets=c['local_grid']['offsets'];grid=sorted([(a,b) for a in offsets for b in offsets],key=lambda p:(p[0]**2+p[1]**2,p))
        for k,offset in enumerate(grid):one(mref[free]+offset,'GRID_'+str(k).zfill(3),'fixed_grid')
        for d in range(c['radial_paths']['directions']):
            direction=np.array([np.cos(2*np.pi*d/8),np.sin(2*np.pi*d/8)]);seed=accepted[0]
            for j in range(1,c['radial_paths']['max_steps']+1):
                values=mref[free]+j*c['radial_paths']['material_free_step']*direction
                entry,row=one(values,'RADIAL_'+str(d)+'_'+str(j).zfill(2),'radial_continuation',seed)
                if entry is None:break
                seed=entry
        summary=dict(parent=parent,method=method,attempts=len(rows),qualified=sum(r['qualified'] for r in rows),grid_qualified=sum(r['qualified'] for r in rows if r['kind']=='fixed_grid'),
            current_setup_wall=setup_wall if method=='B' else 0.,method_setup_wall=setup_wall,common_reference_RHS=common_paid,
            total_actual_RHS=book.run_used,rows=rows,chart_switches=0,
            coverage='Fixed sampled grid plus legal radial segments; no global component coverage certificate',
            chart_switch_status='No unambiguous singular-chart event certified; boundary/stagnation ends preserved, no changed-free-domain counted')
        write(out/'SUMMARY.json',summary);book.end('COMPLETE',qualified=summary['qualified'],attempts=len(rows));print(json.dumps({k:v for k,v in summary.items() if k!='rows'}))
    except BaseException as e:
        write(out/'INTERRUPTED.json',dict(error=repr(e),rows=rows,paid_RHS=book.run_used));book.end('FAILED_OR_BUDGET_STOP',error=repr(e));raise

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--parent',type=int,choices=range(3),required=True);a.add_argument('--method',choices=['A','B'],required=True);v=a.parse_args();run(v.parent,v.method)
