"""Parent-run source-sequential tempered SMC, with a strict local pilot cap.

All particles start from the proper original prior. A paid DEV030 anchor is
read only for an energy-landscape baseline, never as particle initialization.
"""
from pathlib import Path
import argparse,json,os,sys
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ.setdefault(name,'1')
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'tools'))
from p3r2.core import Book,Target,write,save
from p3r2.non_gaussian import source_sequential_smc
from run_p3r2_posterior import finite,freeze


def landscape(prior,q,weights,energies,anchor_energies=None):
    weights=np.asarray(weights,float);weights/=weights.sum();q=np.asarray(q,float)
    fields=prior.config['field_mean']+q[:,:64]@prior.A.T;inside=fields[:,prior.mask]
    host=(inside<=-prior.temperature).mean(axis=1);pure=(inside>=prior.temperature).mean(axis=1)
    row=dict(weighted_total_Phi=float(weights@energies.sum(axis=1)),minimum_total_Phi=float(energies.sum(axis=1).min()),
        weighted_Phi_per_source=weights@energies,minimum_Phi_per_source=energies.min(axis=0),
        weighted_native_property_mean=weights@q[:,64:],native_property_minimum=q[:,64:].min(axis=0),native_property_maximum=q[:,64:].max(axis=0),
        ROI_host_plateau_fraction=float(weights@host),ROI_pure_plateau_fraction=float(weights@pure),
        ROI_transition_fraction=float(weights@(1-host-pure)))
    if anchor_energies is not None: row.update(paid_anchor_Phi_per_source=anchor_energies,paid_anchor_total_Phi=float(anchor_energies.sum()),minimum_full_energy_excess_over_paid_anchor=float(energies.sum(axis=1).min()-anchor_energies.sum()))
    return row


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--job',required=True);parser.add_argument('--parent',default='DEV_030')
    parser.add_argument('--anchor',type=Path,default=ROOT.parent/'b1-p3l-local-posterior/research/p3l_local_posterior_imaging/A/anchors/DEV_030.npz')
    parser.add_argument('--particles',type=int,default=64);parser.add_argument('--mutation-steps',type=int,default=2)
    parser.add_argument('--pcn-beta',type=float,default=.04);parser.add_argument('--property-sd',type=float,default=.02)
    parser.add_argument('--rhs-cap',type=int,default=262144);parser.add_argument('--max-stages',type=int,default=256)
    parser.add_argument('--max-wall',type=int,default=7200);parser.add_argument('--seed',type=int,default=2026102293)
    args=parser.parse_args()
    if args.rhs_cap>262144 or args.rhs_cap<args.particles*8: raise ValueError('bounded new source-SMC pilot required')
    book=Book('B',args.job,max_wall=args.max_wall);status='COMPLETE';error=None;calls=0;events=None
    try:
        target=Target(book);y=target.observation(args.parent);dest=book.jobdir
        anchor_energies=None
        if args.anchor.exists():
            if args.parent!='DEV_030': raise ValueError('DEV030 paid baseline must not be applied to another observation')
            with np.load(args.anchor,allow_pickle=False) as z:
                q=z['q'].copy();prediction=z['prediction'].copy();chi=z['chi'].copy()
            if not np.allclose(target.prior.decode(q)[0],chi,atol=1e-12,rtol=0): raise ValueError('paid baseline decoder identity')
            anchor_energies=np.sum(abs(prediction-y)**2,axis=0)/target.sigma**2
        protocol=dict(job=args.job,stage='B',parent=args.parent,target_id=target.cfg['target_id'],GT_used=False,
            online_observation_source_isolated=True,target='Original P1 Gaussian64 x Uniform2 and all8 FullMaxwell sources',
            source_order=list(range(8)),source_energy='sum_receiver(abs(prediction[:,s]-y[:,s])**2)/sigma**2',
            algorithm='Source-sequential adaptive CESS tempering; all previous source powers1, current sourcebeta, future sources0',
            particle_initialization='Independent normalized original prior draws: z~N(0,I64), u~Uniform([0,1]^2)',
            no_anchor_particle_initialization=True,anchor_baseline_only=str(args.anchor.resolve()),
            original_FULL_only_if_all_source_powers_one=True,log_likelihood_constants_omitted=True,
            evidence_scope='Unnormalized source likelihood constants omitted; no absolute evidence claim',particles=args.particles,mutation_steps=args.mutation_steps,
            cess_target=.8,resample_threshold=.6,pcn_beta=args.pcn_beta,property_sd=args.property_sd,
            mutation='Frozen pCN field plus symmetric native property reflection; exact MH at current source target',
            cost='Each proposal returns all8 exact source energies through a new FullMaxwell forward costing8 RHS',
            pilot_RHS_cap=args.rhs_cap,maximum_full_evaluations=args.rhs_cap//8,max_stages=args.max_stages,
            seed=args.seed,reference_qualification=False,posterior_qualification=False,
            budget_stop='Stop before next full evaluation; save exact last source coefficients and all particle/weight/ancestry arrays',
            no_replacement_of_old_SMC_results=True,paid_anchor_Phi_per_source=anchor_energies)
        freeze(dest/'PROTOCOL.json',protocol);events=(dest/'EVENTS.jsonl').open('x');landscapes=[]
        def source_energy(q):
            nonlocal calls
            if (calls+1)*8>args.rhs_cap: raise RuntimeError('wrapper pilot cap')
            calls+=1;p=target.evaluate(q,purpose='p3r2_source_sequential_exactFull');e=np.sum(abs(p.state.data-y)**2,axis=0)/target.sigma**2
            # identity verifies the source decomposition against the original
            # full likelihood before releasing the paid state.
            if not np.isclose(e.sum(),target.phi(p,y),rtol=2e-13,atol=1e-10): raise ValueError('source energy does not sum to original Phi')
            del p
            if target.model is not None: target.model._transfers.clear()
            return e
        def callback(row):
            events.write(json.dumps(finite(row),allow_nan=False,separators=(',',':'))+'\n');events.flush()
            if row['event'] in ('INITIALIZED','STAGE_COMPLETE','FINISHED'):
                label=f"stage_{row['stage']:04d}_{row['event'].lower()}"
                save(dest/'particles'/f'{label}.npz',q=row['particles'],source_energies=row['source_energies'],
                    logweights=row['logweights'],ancestor_ids=row['ancestor_ids'],target_coefficients=row['target_coefficients'])
                metrics=landscape(target.prior,row['particles'],np.exp(row['logweights']),row['source_energies'],anchor_energies)
                summary=dict(event=row['event'],stage=row['stage'],source=row['source'],beta=row['beta'],
                    target_coefficients=row['target_coefficients'],new_RHS=calls*8,distinct_prior_ancestors=int(np.unique(row['ancestor_ids']).size),
                    **metrics)
                if 'history' in row: summary['mutation']=row['history']['mutation'];summary['weight_ESS_before_resampling']=row['history']['weight_ESS_before_resampling']
                landscapes.append(summary);write(dest/'LANDSCAPE_PROGRESS.json',finite(dict(rows=landscapes,reference_qualification=False)))
                print(json.dumps(finite(summary)),flush=True)
        result=source_sequential_smc(target.prior.draw,source_energy,rng=np.random.default_rng(args.seed),
            N=args.particles,source_count=8,mutation_steps=args.mutation_steps,pcn_beta=args.pcn_beta,
            property_sd=args.property_sd,cess_target=.8,resample_threshold=.6,max_stages=args.max_stages,
            max_evaluations=args.rhs_cap//8,callback=callback)
        save(dest/'PARTICLES.npz',**{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
        chi=np.array([target.prior.decode(q)[0] for q in result['particles']]);w=result['weights'];mean=w@chi;variance=w@abs(chi-mean)**2
        save(dest/'EXPLICIT_TARGET_MATERIAL_MOMENTS.npz',chi=chi,weights=w,mean=mean,variance=variance,target_coefficients=result['target_coefficients'])
        summary={k:v for k,v in result.items() if not isinstance(v,np.ndarray)}
        summary.update(protocol=protocol,new_states=calls,new_RHS=calls*8,
            target_coefficients=result['target_coefficients'],final_landscape=landscapes[-1],
            weighted_native_property_mean=w@result['particles'][:,64:],
            full_target_reached=result['full_target_reached'],posterior_qualification=False,reference_qualification=False,
            target_interpretation='Original FULL target' if result['full_target_reached'] else 'Explicit partial-source power target; not original FULL')
        write(dest/'RESULT.json',finite(summary));print(json.dumps(finite({k:summary[k] for k in ('status','new_RHS','target_coefficients','weight_ESS','distinct_prior_ancestors','target_interpretation')})),flush=True)
    except BaseException as exc:
        status='ERROR';error=repr(exc);write(book.jobdir/'ERROR.json',dict(error=error,new_RHS=calls*8,partial_outputs_preserved=True));raise
    finally:
        if events is not None: events.close()
        book.close(status,error)


if __name__=='__main__': main()
