import json,time
import numpy as np
from numpy.polynomial.legendre import leggauss
from accounting import P,Book,write
from physics import model_for,quality
from charts import MaterialChart,CurrentChart,partition

def prior(m):return float(np.prod(2*np.asarray(m)/(4**2-.5**2)))

def main():
    contract=json.loads((P/'EXPERIMENT_CONTRACT.json').read_text());c=contract['E4_local']
    b=Book('E4_LOCAL_CONDITIONAL_MATERIAL_INTEGRATION');model=model_for(b);mref=np.array(contract['E1']['public_reference_materials'][0])
    ref=model.forward(mref,jac=True,label='LOCAL_REFERENCE');dep,free=partition(ref['J']);target=ref['y']
    ma=MaterialChart(model,target,mref,dep,free);cu=CurrentChart(model,target,ref,dep,free)
    z0=cu.encode(ref['j'],mref);out=P/'results/E4_local';out.mkdir(parents=True,exist_ok=False)
    write(out/'CONDITIONAL_TARGET.json',dict(target='LOCAL-FIBER-BOX under the original P_Z coarea kernel, renormalized only over this declared graph tile',
        reference=mref,dependent=dep,free=free,halfwidth=c['halfwidth'],event='first selected free MATERIAL coordinate exceeds its public reference value; split integration at zero',
        global_missing_mass='NO_BOUND',posterior='NOT_FULL_NOISY_POSTERIOR',future_or_truth_inputs=False))
    results={};allrows=[];lastarrays=None
    try:
        for n in c['orders']:
            node,ww=leggauss(n);width=c['halfwidth'];rows=[]
            # Discontinuous material event is integrated on its two exact half-tiles.
            points=[]
            for half in (-1,1):
                lo,hi=(-width,0.) if half<0 else (0.,width)
                for i,x in enumerate(node):
                    for j,y in enumerate(node):points.append((np.array([(lo+hi)/2+(hi-lo)*x/2,width*y]),(hi-lo)*ww[i]/2*width*ww[j],half))
            mstart=mref.copy();zstart=z0.copy();materials=[];weightsA=[];weightsB=[];baseweights=[];events=[];meansA=[];meansB=[]
            for index,(offset,w,half) in enumerate(points):
                wanted=mref[free]+offset;row=dict(order=n,index=index,offset=offset.tolist(),integration_weight=float(w))
                before=b.run_used;t=time.perf_counter();ra=ma.correct(mstart,wanted,'LOCAL_A_'+str(n)+'_'+str(index))
                if not ra['success']:raise RuntimeError('Direct local reference root not converged')
                JA=model.jacobian(ra['state']);a=ra['m'];wa=prior(a)/abs(np.linalg.det(JA[:,dep]))
                va=model.forward(a,independent=True,label='LOCAL_A_INDEPENDENT');qa=quality(model,va,target,mref)
                row.update(A_RHS=b.run_used-before,A_wall=time.perf_counter()-t,A_pair=qa['pair_relative'])
                before=b.run_used;t=time.perf_counter();rb=cu.correct(zstart,wanted,'LOCAL_B_'+str(n)+'_'+str(index))
                if not rb['success']:raise RuntimeError('Current local reference root not converged')
                z=rb['m'];tangent=cu.material_tangent(rb['z'],wanted)
                vb=model.forward(z,jac=True,label='LOCAL_B_FULLJ_VERIFICATION');JB=vb['J']
                area=np.sqrt(np.linalg.det(tangent.T@tangent));normal=np.prod(np.linalg.svd(JB,compute_uv=False));wb=prior(z)*area/normal
                vi=model.forward(z,independent=True,label='LOCAL_B_INDEPENDENT');qb=quality(model,vi,target,mref)
                row.update(B_RHS=b.run_used-before,B_wall=time.perf_counter()-t,B_pair=qb['pair_relative'],material_pair_RMS=float(np.linalg.norm(a-z)/np.sqrt(6)),direct_density=wa,current_density=wb,area=float(area),normal_Jacobian=float(normal))
                if not qa['qualified'] or not qb['qualified']:raise RuntimeError('Local probability node failed independent physics gate')
                materials.append(a);weightsA.append(w*wa);weightsB.append(w*wb);baseweights.append(w);events.append(float(half>0));mstart=a;zstart=rb['z']
                np.savez_compressed(out/('NODE_'+str(n)+'_'+str(index).zfill(3)+'.npz'),m_A=a,m_B=z,j_B_chart=rb['j'],j_A_full=va['j'],j_B_full=vi['j'],y_A=va['y'],y_B=vi['y'],J_A=JA,J_B=JB,tangent_B=tangent,target=target)
                rows.append(row)
            mats=np.array(materials);ea=np.array(events);wA=np.array(weightsA);wB=np.array(weightsB);w0=np.array(baseweights)
            def summary(w):
                return dict(coarea_tile_density_mass=float(w.sum()),mean=(w@mats/w.sum()).tolist(),second_moment=(w@(mats**2)/w.sum()).tolist(),event_probability=float(w@ea/w.sum()))
            results[str(n)]=dict(A=summary(wA),B=summary(wB),wrong_uniform_free=summary(w0),wrong_equal_nodes=summary(np.ones(len(mats))),
                node_count=len(mats),maximum_density_relative_discrepancy=float(np.max(abs(wA-wB)/wA)),maximum_material_pair_RMS=max(r['material_pair_RMS'] for r in rows),rows=rows)
            lastarrays=dict(m=mats,weights=wA/wA.sum(),free=free,dep=dep)
            allrows.extend(rows);write(out/'CHECKPOINT.json',dict(completed_orders=list(results),paid=b.run_used))
        high=str(c['orders'][-1]);previous=str(c['orders'][-2]);s=results[high]['A'];r=results[previous]['A'];bhi=results[high]['B']
        err=dict(mass_relative=abs(s['coarea_tile_density_mass']/r['coarea_tile_density_mass']-1),mean_max=float(np.max(abs(np.array(s['mean'])-r['mean']))),event_abs=abs(s['event_probability']-r['event_probability']))
        cross=dict(mass_relative=abs(s['coarea_tile_density_mass']/bhi['coarea_tile_density_mass']-1),mean_max=float(np.max(abs(np.array(s['mean'])-bhi['mean']))),event_abs=abs(s['event_probability']-bhi['event_probability']))
        gate=all(v['mass_relative']<=1e-4 and v['mean_max']<=1e-4 and v['event_abs']<=1e-3 for v in (err,cross))
        # Controlled independentization of the SAME finite weighted quadrature law.
        rng=np.random.default_rng(2026102904);inds=rng.choice(len(lastarrays['m']),size=16,p=lastarrays['weights']);joint=lastarrays['m'][inds]
        product=joint.copy()
        for i in range(6):product[:,i]=joint[rng.permutation(len(joint)),i]
        bad=[]
        for i,m in enumerate(product[:8]):
            state=model.forward(m,independent=True,label='INDEPENDENTIZED_LOCAL_CONTROL');bad.append(quality(model,state,target,mref))
        np.savez_compressed(out/'INDEPENDENTIZED_CONTROL.npz',joint=joint,product=product)
        result=dict(status='LOCAL_NUMERICAL_REFERENCE_PASS' if gate else 'LOCAL_REFERENCE_UNQUALIFIED',orders=results,
            convergence=err,cross_chart=cross,global_missing_mass='UNBOUNDED',global_prior_support='unchanged P_Z, but reported kernel explicitly conditioned to this tile',
            original_noisy_posterior='NOT_COMPUTED',independentization=dict(scope='Controlled finite quadrature law; NOT claimed exact full conditional prior marginals',proposals=bad),paid=b.run_used)
        write(out/'SUMMARY.json',result);b.end(result['status']);print(json.dumps({k:v for k,v in result.items() if k not in ('orders','independentization')},indent=2))
    except BaseException as exc:
        write(out/'INTERRUPTED.json',dict(error=repr(exc),completed_orders=results,rows=allrows,paid=b.run_used));b.end('FAILED',error=repr(exc));raise
if __name__=='__main__':main()
