"""Metered selected endpoint packaging and independent-sensor experiment."""
from pathlib import Path
import sys,json,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.b1_joint.accounting import Ledger,dump
from src.b1_joint.forward import MaterialForward
from src.b1_joint.response import Response
from src.b1_joint.reference import Reference
from src.b1_joint.material import event_ids,void_volume,geometry
from src.b1_joint.likelihood import context_loglik
from src.b1_joint.package import build_packages,update_context

def decisions(package,x,w,c,ref_event,truth_event):
    updated=update_context(x,w,lambda z:context_loglik(z,c))
    shown=np.array([a['event_id'] for a in package['hypotheses']]);rep=np.array([a['representative_latent'] for a in package['hypotheses']]);p=np.array([a['posterior_mass'] for a in package['hypotheses']])
    values=np.log(p)+context_loglik(rep,c);top=int(shown[np.argmax(values)])
    variants=dict(single_MAP=int(package['single_MAP']['event_id']),TopK_representatives=top,whole_weighted_posterior=int(np.argmax(updated['event_mass'])),TopK_conditional_samples_plus_remainder=int(np.argmax(updated['event_mass'])))
    rows={k:dict(selected_event=v,posterior_decision_regret=float(ref_event.max()-ref_event[v]),realized_zero_one_loss=int(v!=truth_event)) for k,v in variants.items()}
    # The rich package retains all conditional particles, including the tail;
    # its equality to the whole sample update is a representation identity.
    return rows,dict(event_mass=updated['event_mass'],ESS=updated['ESS'],max_event_error=float(abs(updated['event_mass']-ref_event).max()),TopK_point_update_status='approximate action only; unknown context likelihood of unshown mass means no normalized full posterior claimed',rich_package_identity=True)

def main():
    out=ROOT/'research/joint_p2/results';ledger=Ledger(out);response=Response(out/'COMMON_RESPONSE_BANK.npz');ref=Reference(response);forward=MaterialForward(ROOT,ledger);ledger.stage='fresh'
    with ledger.job('p2-context-partition-refinement-resume'):
        paths=list((out/'packages').glob('*.json'));packages={p.stem:json.loads(p.read_text()) for p in paths}
        if len(packages)!=7:packages=build_packages(out,response,ref,forward)
        rng=np.random.default_rng(261013217);cases={};raw=json.loads((out/'DEVELOPMENT_CASES.json').read_text())['cases'];corrected=json.loads((out/'FULL_MIE_CORRECTED_REFERENCE.json').read_text())['cases']
        frozen=out/'CONTEXT_SENSOR_MANIFEST.json'
        if not frozen.exists():
            # Preserve the exact stream already used before the first context
            # refinement failure; freeze it and separate all bootstrap RNGs.
            sensor_rng=np.random.default_rng(26101328);values=[]
            for item in raw:
                values.append(float(void_volume(np.array(item['latent']))+.01*sensor_rng.normal()))
                if item['id'] in corrected:
                    n=corrected[item['id']]['n_per_event']
                    for j in range(12):sensor_rng.choice(np.zeros(n),(2000,n),replace=True)
            dump(frozen,dict(values=values,case_order=[a['id'] for a in raw],seed=26101328,stream='original pre-failure sequence retained; independent bootstrap stream from resume',sigma=.01,no_selection=True))
        contexts=json.loads(frozen.read_text())['values']
        for item in raw:
            name=item['id'];truth=np.array(item['latent']);c=contexts[len(cases)];y=np.array(item['y'])[:32];sd=np.array(item['sigma'])[:32];r=ref.solve(y,sd,order=17,context=c);adaptive=ref.adaptive_physical(r,max_cells=8192);p=np.array(adaptive['event_mass']);numerical='proxy continuous independent adaptive reference';MC=None
            if not adaptive['converged']:raise RuntimeError('context reference did not converge')
            if name in corrected:
                a=np.load(out/f'{name}_FULL_CORRECTION.npz');x=a['latent'].reshape(12,-1,4);d=a['residual_weights'];sensor=np.exp(context_loglik(x.reshape(-1,4),c)).reshape(d.shape);base=p*np.exp(adaptive['log_evidence']-float(a['likelihood_shift']));z=base+(d*sensor).mean(1)
                if np.any(z<=0):raise FloatingPointError('context full-residual event evidence unresolved')
                p=z/z.sum();boots=np.column_stack([base[j]+rng.choice(d[j]*sensor[j],(2000,d.shape[1]),replace=True).mean(1) for j in range(12)])
                if np.any(boots<=0):raise FloatingPointError('context residual bootstrap nonpositive')
                bp=boots/boots.sum(1)[:,None];MC=dict(SE=bp.std(0,ddof=1),interval95=np.quantile(bp,[.025,.975],axis=0).T,scope='reweighted full-minus-proxy residual MC only; cubature separate');numerical='full-Mie residual-control-variate reference, paid points reused with independent context likelihood'
            pkg=packages[name];sample=np.load(out/'packages'/f'{name}_samples.npz');decision,approx=decisions(pkg,sample['latent'],sample['weights'],c,p,int(event_ids(truth)))
            row=dict(sensor_value=c,sensor_sigma=.01,true_void_volume=float(void_volume(truth)),event_mass_before=pkg['all_event_mass'],event_mass_after=p,context_independence='declared simulated independent noisy geometry sensor given material; no y reuse, truth not supplied to inference',numerical_status=numerical,context_reference=adaptive,physical_MC=MC,decisions=decision,weighted_sample_update_diagnostics=approx,new_physics=0)
            cases[name]=row;pkg['downstream_context_update']['actual_demonstration']=row;dump(out/'packages'/f'{name}.json',pkg)
            print(json.dumps(dict(case=name,context=c,before=pkg['all_event_mass'],after=p.tolist(),decisions=decision)),flush=True)
        # Fresh parents: reference inference has no learned model-selection role.
        final=np.load(out/'C1_FINAL_DATA.npz');truth=final['truth'];parents=final['parent_ids'];panel=[]
        ref_rng=np.random.default_rng(261013216);contexts=void_volume(truth)+.01*ref_rng.normal(size=len(truth))
        for i,yy in enumerate(final['y']):
            c=float(contexts[i]);r=ref.solve(yy[:32],sd,order=13,context=c);adaptive=ref.adaptive_physical(r,max_cells=8192)
            if not adaptive['converged']:raise RuntimeError('fresh context reference unresolved, stop claims')
            p=np.array(adaptive['event_mass']);x=final['reference_samples'][i];w=np.full(len(x),1/len(x));before=final['reference_event_mass'][i];eid=event_ids(x)
            # Exact base event masses, node approximation only within event.
            for j in range(12):
                mask=eid==j
                if mask.any():w[mask]=before[j]/mask.sum()
            w/=w.sum();update=update_context(x,w,lambda z:context_loglik(z,c));shown=np.argsort(-before);shown=shown[:min(12,1+np.searchsorted(np.cumsum(before[shown]),.95))];rep=[]
            for j in shown:
                a=x[eid==j];rep.append(a[np.argmin(np.sum(((np.column_stack((geometry(a),a[:,2:]))-np.column_stack((geometry(a),a[:,2:])).mean(0))/np.array([.3,.12,.6,.12]))**2,axis=1))])
            point=int(shown[np.argmax(np.log(before[shown])+context_loglik(np.array(rep),c))]);map_event=int(event_ids(x[np.argmax(final['reference_log_density'][i])]))
            choices=dict(single_MAP_particle=map_event,TopK_representatives=point,whole_weighted_posterior=int(np.argmax(update['event_mass'])),TopK_conditional_samples_plus_remainder=int(np.argmax(update['event_mass'])),reference_optimal=int(np.argmax(p)),no_new_information=int(np.argmax(before)))
            panel.append(dict(observation=i,parent=int(parents[i]),sensor_value=c,truth_event=int(event_ids(truth[i])),event_mass_before=before,event_mass_after=p,reference_converged=adaptive['converged'],estimated_reference_error=adaptive['estimated_relative_error'],sample_update_ESS=update['ESS'],sample_event_max_error=float(abs(update['event_mass']-p).max()),choices={k:dict(selected_event=v,posterior_decision_regret=float(p.max()-p[v]),realized_zero_one_loss=int(v!=int(event_ids(truth[i])))) for k,v in choices.items()}))
            if i%16==0:print(f'fresh independent context {i}/{len(truth)}',flush=True)
        # Paired bootstrap unit is parent, never two correlated observations.
        methods=list(panel[0]['choices']);bootstrap={};brng=np.random.default_rng(26101329)
        for m in methods:
            for metric in ['posterior_decision_regret','realized_zero_one_loss']:
                values=np.array([np.mean([a['choices'][m][metric] for a in panel if a['parent']==j]) for j in np.unique(parents)]);sample=brng.choice(values,(2000,len(values)),replace=True).mean(1);bootstrap.setdefault(m,{})[metric]=dict(mean=float(values.mean()),CI95=np.quantile(sample,[.025,.975]).tolist())
        delta=np.array([np.mean([a['choices']['whole_weighted_posterior']['posterior_decision_regret']-a['choices']['single_MAP_particle']['posterior_decision_regret'] for a in panel if a['parent']==j]) for j in np.unique(parents)]);pair=brng.choice(delta,(2000,len(delta)),replace=True).mean(1)
        dump(out/'DOWNSTREAM.json',dict(cases=cases,fresh_panel=panel,summary=bootstrap,paired_weighted_minus_MAP_regret=dict(mean=float(delta.mean()),CI95=np.quantile(pair,[.025,.975]).tolist()),loss='zero-one loss for frozen twelve material events; regret against declared context-updated reference',fresh_MAP_status='highest latent-density saved quadrature particle, not continuous/global MAP',sensor_seed_demo=26101328,sensor_seed_fresh=261013216,bootstrap_seed=26101329,new_physics_for_context=0,conditional_independence=True,no_context_identity='update_context callback identically0 preserves whole weights to roundoff; no additional information cannot resolve ambiguity',rich_package_and_whole_particle_equality='same support/weights retained, not an independent algorithm advantage'))
    dump(out/'RESOURCE_SUMMARY.json',ledger.summary())
if __name__=='__main__':main()
