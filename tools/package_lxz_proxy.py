"""Expose the actual corrected proxy chains as material-event packages.

No full physics. Representative predictions make a separately counted proxy
response call, never rebranded as new Maxwell validation or posterior samples.
"""
from pathlib import Path
import json,sys,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from lxz.proxy import OUT,PAID,IDENTITY,dump
from b1_joint.material import event_ids,prior_event_masses,decode,image,geometry
from b1_joint.response import Response

def main():
    if (OUT/'PROXY_PACKAGE_RECEIPT.json').exists():raise RuntimeError('already built; no silently repeated proxy evaluation')
    start=time.perf_counter();response=Response(PAID/'COMMON_RESPONSE_BANK.npz');count=0
    data=json.loads((PAID/'DEVELOPMENT_CASES.json').read_text())['cases']
    seeds=json.loads((OUT/'configs/PROTOCOL.json').read_text())['E1']['seeds'];prior=prior_event_masses()
    for case in ['solid','large_cavity','annular_cavity','combined']:
        chunks=[];lik=[];perseed=[]
        for seed in seeds:
            z=np.load(OUT/f'results/round1/{case}_global_{seed}.npz');burn=len(z['latent'])//4
            x=z['latent'][burn:];chunks.append(x);lik.append(z['log_likelihood'][burn:]);perseed.append(np.bincount(event_ids(x),minlength=12)/len(x))
        x=np.concatenate(chunks);ll=np.concatenate(lik);ids=event_ids(x);p=np.mean(perseed,axis=0);eorder=np.argsort(-p)
        selected=[];represented=0.
        for ev in eorder:
            if p[ev]==0:continue
            selected.append(int(ev));represented+=p[ev]
            if represented>=.95:break
        hy=[];reps=[];item=next(i for i in data if i['id']==case)
        for ev in selected:
            where=np.flatnonzero(ids==ev);best=where[np.argmax(ll[where])];v=x[where];g=np.c_[geometry(v),v[:,2:]]
            representative=x[best];reps.append(representative);desc=decode(representative)
            f=response.predict(representative,(1.6,))[0];count+=1
            estimates=np.array(perseed)[:,ev]
            hy.append(dict(event_id=ev,event_definition=f'original central-radius bin{ev//3} and annular-halfwidth bin{ev%3}; disjoint original P2 partition',
                posterior_mass=float(p[ev]),mass_error=dict(type='four chain-seed descriptive range, not confidence or numerical certificate',range=estimates[[np.argmin(estimates),np.argmax(estimates)]],
                    independent_reference='ADAPTIVE_PHYSICAL_REFERENCE.json; finite cubature error separate'),prior_mass=float(prior[ev]),
                posterior_to_prior_odds_factor=float((p[ev]/max(1-p[ev],1e-300))/(prior[ev]/max(1-prior[ev],1e-300))),
                representative=dict(latent=representative,radii=desc['radii'],epsilon_real=desc['eps_layers'].real,epsilon_imag=desc['eps_layers'].imag,
                    rule='maximum saved likelihood sample within event, no truth'),
                material_uncertainty=dict(coordinates=['central_radius','annular_halfwidth','epsilon_real','epsilon_imag'],mean=g.mean(0),SD=g.std(0),quantiles=np.quantile(g,[.05,.5,.95],axis=0)),
                forward_prediction=dict(target=IDENTITY,kind='paid-bank proxy response, not fresh full-Mie',real_channels=f),
                full_Maxwell_validation=None,separated_mode_status='NOT_ESTABLISHED',posterior_image_point_mass=None))
        target=OUT/'MaterialHypothesisPackage'/f'PROXY_{case}.json'
        dump(target,dict(observation_id=case,target_identity=IDENTITY,observation=np.array(item['y'])[:32],sigma=np.array(item['sigma'])[:32],
            inference='global corrected independence MH, four seeds, each4096 likelihood including64setup,25percent burn',hypotheses=hy,
            represented_mass=float(represented),unshown_modeled_mass=float(1-represented),unexplored_mass_status='finite-chain unquantified; remainder is empirical modeled-event mass, not bound on missed modes',
            outside_model_support_mass=None,numerical_status='PROXY_TARGET_MECHANISM_ONLY',
            posterior_samples=f'../results/PROXY_{case}_SAMPLES.npz',no_posterior_renormalization_of_displayed_K=True))
        np.savez_compressed(OUT/f'results/PROXY_{case}_SAMPLES.npz',latent=x,weights=np.full(len(x),1/len(x)),event_ids=ids,log_likelihood=ll)
        fig,axes=plt.subplots(2,len(reps),figsize=(3.2*len(reps),5),squeeze=False,layout='constrained')
        for j,(rep,h) in enumerate(zip(reps,hy)):
            axis,eps=image(rep)
            for row,field,limits in [(0,eps.real,(1,2.8)),(1,eps.imag,(0,.16))]:
                axes[row,j].imshow(field,origin='lower',extent=[-1.05,1.05,-1.05,1.05],vmin=limits[0],vmax=limits[1],cmap='viridis' if row==0 else 'magma')
                axes[row,j].set_xlabel('x/R');axes[row,j].set_ylabel('y/R')
            axes[0,j].set_title(f"event {h['event_id']}: {h['posterior_mass']:.3f}\nprior {h['prior_mass']:.3f}")
        fig.suptitle(f'{case}: corrected PROXY event representatives, not separated modes\nRe(epsilon) top / Im(epsilon) bottom; displayed mass={represented:.3f}; remainder={1-represented:.3f}')
        fig.savefig(OUT/f'figures/proxy_{case}.png',dpi=150);fig.savefig(OUT/f'figures/proxy_{case}.pdf');plt.close(fig)
    with (OUT/'COST_LEDGER.jsonl').open('a') as f:f.write(json.dumps(dict(event='DISPLAY_PROXY_RESPONSE',auxiliary_proxy_response_attempts=count,proxy_likelihood_attempts=0,new_full_states=0,new_RHS=0,wall_seconds=time.perf_counter()-start))+'\n')
    dump(OUT/'PROXY_PACKAGE_RECEIPT.json',dict(auxiliary_proxy_responses=count,new_full=0,posterior_resampling=0,source='savedcorrectedchains',wall_seconds=time.perf_counter()-start))
    p=OUT/'GALLERY.html';s=p.read_text()
    for case in ['solid','large_cavity','annular_cavity','combined']:
        s=s.replace('</html>',f'<section><h2>{case}：proxy材料事件代表</h2><p>有事件质量不等于分离mode，未展示质量保留，fullMie未新增验证。<a href="MaterialHypothesisPackage/PROXY_{case}.json">机器可读概率与样本</a></p><img src="figures/proxy_{case}.png" alt="proxy材料事件"></section></html>')
    p.write_text(s)

if __name__=='__main__':main()
