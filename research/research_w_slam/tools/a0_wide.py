"""Optional four-parent unknown 3x3 support extension; no work on import."""
from pathlib import Path
import argparse
import gc
import json
import time

import numpy as np

from wslam.accounting import BudgetBook, BudgetExceeded, json_default
from wslam.forward import Evaluator
from wslam.geometry import enumerate_supports, validate_geometry
from wslam.physics import acquisition, pack
from wslam.world import World
from a0_run import common_pool, run_fit, write_json, save_npz
from a0_score import offline_metrics, paid_prediction, complex_json


def generate_public_prior_data(root,out,book,supplement):
    """Generator-only function: returns no labels or stratum to the inverse."""
    online,offline=out/'online',out/'offline'
    online.mkdir();offline.mkdir()
    choices=enumerate_supports(3,.65)
    groups={True:[],False:[]}
    for support,mass in choices:
        groups[validate_geometry(support)['holes']>=1].append((support,mass))
    if not all(groups.values()):raise RuntimeError('Both frozen public hole strata must be nonempty')
    rng=np.random.default_rng(supplement['truth_seed'])
    noise=np.random.default_rng(supplement['noise_seed'])
    tx,rx=acquisition(outer=True)
    observations=[];labels=[];means=[]
    def truncated():
        while True:
            value=float(rng.normal(0,.005))
            if abs(value)<=.01:return value
    for parent,stratum in enumerate(supplement['parent_strata']):
        options=groups[stratum=='at_least_one_hole']
        masses=np.array([p for _,p in options]);masses/=masses.sum()
        support=options[int(rng.choice(len(options),p=masses))][0].copy()
        world=World(support,rng.uniform(1.5,6),rng.uniform(0,.25),
                    [rng.uniform(-1,1),rng.uniform(1.6,2.6)],radius=.12,contrast=3.,target_loss=.02)
        gain,phase=truncated(),truncated()
        theta=np.r_[world.epsinf,world.ell,world.pose,gain,phase]
        mean,metadata=paid_prediction(world,theta,tx,rx,book)
        y=mean+.01*(noise.normal(size=mean.shape)+1j*noise.normal(size=mean.shape))
        observations.append(y);means.append(mean)
        labels.append(dict(parent_id=parent,truth_stratum=stratum,**world.metadata(),
                           instrument=dict(amplitude_offset=gain,phase=phase),solver_metadata=metadata))
        write_json(offline/'TRUTHS.json',labels)
        save_npz(offline/'CLEAN_MEANS.npz',dict(data=np.stack(means)))
        # All real observations, not idealized curves, remain inspectable.
        save_npz(online/'BASE.npz',dict(base=np.stack(observations)[:,:2],parent_id=np.arange(len(observations)),noise_repeat=np.zeros(len(observations),int),tx=tx[:2],rx=rx))
        save_npz(offline/'HOLDOUT.npz',dict(data=np.stack(observations)[:,2:],tx=tx[2:],rx=rx))
    return dict(parents=4,repeats=1,online_path=str(online/'BASE.npz'))


def fit_online_only(root,out,book):
    """This loader reads online measurements and shared online bank only."""
    with np.load(root/'inputs/online/COARSE_BANK.npz',allow_pickle=False) as raw:
        if int(raw['completed'])!=int(raw['total']):raise RuntimeError('Shared bank incomplete')
        bank=raw['data'][:,:2].copy()
    with np.load(root/'inputs/online/COARSE_GRID.npz',allow_pickle=False) as raw:grid={k:raw[k].copy() for k in raw.files}
    with np.load(out/'online/BASE.npz',allow_pickle=False) as raw:
        data=raw['base'].copy();ids=raw['parent_id'].copy()
    rows=[dict(parent_id=int(pid),status='NOT_RUN',selection=None) for pid in ids]
    write_json(out/'ONLINE_RESULTS.json',dict(status='ONLINE_RUNNING',rows=rows,offline_read=False))
    stopped=False
    for y,row in zip(data,rows):
        if stopped:
            row['reason']='SHARED_A0_BUDGET_STOP';continue
        directory=out/f"fit_parent_{row['parent_id']:02d}";directory.mkdir()
        try:
            pool,_,metadata=common_pool(y,bank,grid)
            write_json(directory/'INITIALIZER_POOL.json',dict(pool=pool,metadata=metadata))
            result=run_fit(y,pool[0],6,book,directory,'FULL_J_TR')
            row.update(status=result['fit_status'],selection=result if result['fit_status']!='FAILED' else None,result=result)
        except BudgetExceeded as exc:
            row.update(status='STOPPED_BUDGET',error=str(exc));stopped=True
        except Exception as exc:
            row.update(status='FAILED',error_type=type(exc).__name__,error=str(exc))
        write_json(out/'ONLINE_RESULTS.json',dict(status='ONLINE_RUNNING',rows=rows,offline_read=False))
    write_json(out/'ONLINE_RESULTS.json',dict(status='ALL_ONLINE_TERMINAL_DURABLE',rows=rows,offline_read=False,budget_stopped=stopped))
    return stopped


def image_actual(path,metrics,label,selection,base_y,base_pred,outer_y,outer_pred):
    """Offline-only material maps and actual measured/predicted components."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    theta=np.asarray(selection['theta'])
    truth=validate_geometry(label['support'],label['epsinf'],label['ell'])
    candidate=validate_geometry(selection['support'],theta[0],theta[1])
    fig,axes=plt.subplots(2,3,figsize=(13,7),layout='constrained')
    for ax,field,title in [(axes[0,0],truth['epsinf'],'Offline true epsilon_inf'),
                           (axes[0,1],candidate['epsinf'],'Selected epsilon_inf')]:
        ax.imshow(field,origin='lower',extent=(-1.5,1.5,.5,1.25),vmin=1,vmax=6,aspect='auto');ax.set_title(title)
    axes[0,2].imshow(candidate['support'].astype(int)-truth['support'].astype(int),origin='lower',vmin=-1,vmax=1,cmap='coolwarm',aspect='auto')
    axes[0,2].set_title('Support difference: selected minus true')
    for ax,observed,predicted,title in [(axes[1,0],base_y,base_pred,'Base actual measurement'),(axes[1,1],outer_y,outer_pred,'Independent outer measurement')]:
        ax.plot(observed.ravel().real,label='observed real');ax.plot(predicted.ravel().real,'--',label='predicted real')
        ax.plot(observed.ravel().imag,label='observed imag');ax.plot(predicted.ravel().imag,'--',label='predicted imag');ax.set_title(title);ax.legend(fontsize=6)
    axes[1,2].scatter(*label['target_pose'],marker='x',label='offline true target')
    axes[1,2].scatter(*theta[2:4],marker='o',label='selected target');axes[1,2].set_xlim(-1,1);axes[1,2].set_ylim(1.6,2.6);axes[1,2].legend(fontsize=7)
    fig.suptitle('Exploratory 3x3 parent; raw observations and selected fit; no posterior claim')
    fig.savefig(path,dpi=150);plt.close(fig)


def score_after_all_online(root,out,book):
    saved=json.loads((out/'ONLINE_RESULTS.json').read_text())
    if saved['status']!='ALL_ONLINE_TERMINAL_DURABLE':raise RuntimeError('All online fits must be saved before labels are opened')
    labels=json.loads((out/'offline/TRUTHS.json').read_text())
    label_by_id={int(r['parent_id']):r for r in labels}
    with np.load(out/'online/BASE.npz',allow_pickle=False) as raw:base=raw['base'].copy()
    with np.load(out/'offline/HOLDOUT.npz',allow_pickle=False) as raw:
        outer=raw['data'].copy();tx=raw['tx'].copy();rx=raw['rx'].copy()
    rule=json.loads((root/'results/a0/preparation/DISCREPANCY_RULE.json').read_text())
    upper=float(rule['upper_noise_energy'])
    rows=[dict(parent_id=r['parent_id'],fit_status=r['status'],selected_fit_success=r.get('selection',{}).get('success') if r.get('selection') else None,
               status='NOT_RUN' if r.get('selection') else 'NO_SELECTED_CANDIDATE',finer_status='NOT_RUN',online_error=r.get('error')) for r in saved['rows']]
    stopped=saved.get('budget_stopped',False)
    for row,fitrow in zip(rows,saved['rows']):
        selection=fitrow.get('selection');pid=row['parent_id']
        if selection is None:continue
        if stopped:row['reason']='SHARED_A0_BUDGET_STOP';continue
        before=book.totals('A0')
        try:
            theta=np.asarray(selection['theta'])
            world=World(selection['support'],theta[0],theta[1],theta[2:4],radius=.12,contrast=3.,target_loss=.02)
            pred,meta=paid_prediction(world,theta,tx,rx,book)
            train_pred=complex_json(selection['pred']);train_white=pack(train_pred-base[pid]);white=pack(pred-outer[pid])
            train_energy=float(train_white@train_white);energy=float(white@white)
            row.update(offline_metrics(selection,label_by_id[pid]),status='SCORED',truth_stratum=label_by_id[pid]['truth_stratum'],
                training_noise_energy=train_energy,holdout_noise_energy=energy,
                training_noise_envelope_compatible=train_energy<=upper,holdout_noise_envelope_compatible=energy<=upper,
                selected_fit_nonconvergence_preserved=not selection['success'],primary_metadata=meta)
            save_npz(out/f'parent_{pid:02d}_primary.npz',dict(theta=theta,support=np.asarray(selection['support']),base_observed=base[pid],base_predicted=train_pred,
                outer_observed=outer[pid],outer_predicted=pred,train_white=train_white,holdout_white=white))
            try:
                image_actual(out/f'parent_{pid:02d}_actual.png',row,label_by_id[pid],selection,base[pid],train_pred,outer[pid],pred)
                row['image_status']='SAVED'
            except Exception as image_error:
                row.update(image_status='FAILED',image_error=str(image_error))
            fine,fine_meta=paid_prediction(world,theta,tx,rx,book,finer=True)
            fine_white=pack(fine-outer[pid]);difference=pack(fine-pred)
            save_npz(out/f'parent_{pid:02d}_independent.npz',dict(primary=pred,independent=fine,difference_white=difference,holdout_white=fine_white,
                metadata_json=np.array(json.dumps(fine_meta,default=json_default))))
            row.update(finer_status='SCORED',finer_prediction_difference_white=float(np.linalg.norm(difference)),finer_holdout_noise_energy=float(fine_white@fine_white))
        except BudgetExceeded as exc:
            stopped=True;row.update(error=str(exc),budget_stopped=True)
            if row['status']=='NOT_RUN':row['status']='STOPPED_BUDGET'
        except Exception as exc:
            row.update(error_type=type(exc).__name__,error=str(exc))
            if row['status']=='NOT_RUN':row['status']='FAILED'
            elif row['finer_status']=='NOT_RUN':row['finer_status']='FAILED'
        finally:
            row['score_stage_cost']={k:book.totals('A0')[k]-before[k] for k in before}
            write_json(out/'OFFLINE_SCORES.json',dict(rows=rows,noise_rule=rule,offline_read_after_all_online=True,budget_stopped=stopped,
                larger_ROI='NOT_RUN',general_A0_PASS=False,posterior_claim=False))
    return stopped


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);args=parser.parse_args();root=args.root.resolve()
    supplement=json.loads((root/'A0_WIDE_SUPPLEMENT.json').read_text())
    out=root/'results/a0/wide_extension';out.mkdir(parents=True,exist_ok=False)
    book=BudgetBook(root/'results/budget_ledger.json')
    remaining=book.snapshot('A0')['A0']['remaining']
    report=dict(status='NOT_RUN',supplement=supplement,headroom=remaining,larger_ROI='NOT_RUN',general_A0_PASS=False,posterior_claim=False,
                planned_online_fits=[dict(parent_id=i,status='NOT_RUN') for i in range(4)])
    prerequisites=[root/'results/a0/scoring/RESULTS.json',root/'results/a0/intrinsic_stress/RESULTS.json']
    pending=[]
    for path in prerequisites:
        if not path.exists():pending.append(str(path));continue
        status=json.loads(path.read_text()).get('status','')
        if status in ('','RUNNING','INITIALIZED','ONLINE_RUNNING','ALL_ONLINE_TERMINAL_DURABLE','OFFLINE_SCORING'):pending.append(str(path))
    if pending:
        report.update(reason='MAIN_AND_UNKNOWN_INTRINSIC_MUST_BE_TERMINAL_FIRST',pending=pending);write_json(out/'RESULTS.json',report);return
    if remaining['rhs']<500 or remaining['factorizations']<40 or remaining['process_cpu_seconds']<=0:
        report['reason']='INSUFFICIENT_SHARED_A0_HEADROOM_NO_BORROWING';write_json(out/'RESULTS.json',report);return
    if json.loads((root/'results/e0/gate.json').read_text()).get('decision')!='PASS':raise RuntimeError('E0 PASS required')
    began=time.perf_counter();report.update(status='GENERATING',budget_start=book.totals('A0'));write_json(out/'RESULTS.json',report)
    try:
        report['generation']=generate_public_prior_data(root,out,book,supplement)
        report['status']='ONLINE_FITTING';write_json(out/'RESULTS.json',report)
        stopped=fit_online_only(root,out,book)
        terminal=json.loads((out/'ONLINE_RESULTS.json').read_text())
        report['planned_online_fits']=[dict(parent_id=r['parent_id'],status=r['status']) for r in terminal['rows']]
        report['status']='ALL_ONLINE_TERMINAL_DURABLE';write_json(out/'RESULTS.json',report)
        stopped=score_after_all_online(root,out,book) or stopped
        report['status']='STOPPED_A0_BUDGET' if stopped else 'EXPLORATORY_WIDE_EXECUTION_COMPLETE_PARENT_REVIEW_REQUIRED'
    except BudgetExceeded as exc:report.update(status='STOPPED_A0_BUDGET',error=str(exc))
    except Exception as exc:report.update(status='FAILED',error_type=type(exc).__name__,error=str(exc))
    finally:
        gc.collect();report.update(budget_end=book.totals('A0'),wall_seconds=time.perf_counter()-began)
        write_json(out/'RESULTS.json',report)
        write_json(out/'gate.json',dict(decision='INCONCLUSIVE',execution_status=report['status'],general_A0_PASS=False,
            scope='Four exploratory stratified 3x3 parents, one noise draw each; larger ROI NOT_RUN'))


if __name__=='__main__':main()
