"""Frozen classic Austria extension: prepare, queue, solve, seal, offline score.

No imports of truth helpers occur in solve or queue. Existing R8 public Y is
read directly. Scoring checks every requested task before constructing truth.
"""
import time
PROCESS_STARTED, PROCESS_CPU_STARTED = time.perf_counter(), time.process_time()
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
sys.path.insert(0, str(REPO))
OUT = ROOT/'RESULTS/classic_extension'
CONFIG = ROOT/'configs/classic_extension_main.json'
# Preserve the active environment entrypoint; resolving a venv symlink loses
# its import environment on this host. This also makes packaged queue portable.
PYTHON = sys.executable
ARMS = ['cie_single','cie_dual','global_dual','a_material','a_plus','physical_lbfgs',
        'paper_nie','paper_lsie']


def dump(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(obj, indent=2, allow_nan=False))
    temporary.replace(path)


def freeze(args):
    if CONFIG.exists():
        raise FileExistsError('Frozen manifest exists; do not overwrite: '+str(CONFIG))
    cases = []
    r8_repo = Path(args.r8_repo) if args.r8_repo else REPO
    r8_study = r8_repo/'research/research_max_struct_r8'
    for key, eps, beta, mf in [
        ('austria_eps3p5_400MHz',3.5,[6,2,.5],[5,7,10]),
        ('austria_eps4_400MHz',4.,[8,2,.5],[5,8,10]),
    ]:
        source = r8_study/'RESULTS/paper_benchmarks'/key/'PUBLIC_DATA.npz'
        if not source.exists(): raise FileNotFoundError(source)
        target = ROOT/'data_public/classic_extension'/key/'PUBLIC_DATA.npz'
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(source,target)
        shutil.copy2(source.parent/'DATA_RECEIPT.json',target.parent/'R8_DATA_RECEIPT.json')
        cases.append(dict(case_id=key, epsilon_r=eps, frequency_Hz=400e6,
                          beta=beta, mf=mf, public_observation_file=str(target.relative_to(REPO)),
                          original_R8_provenance_file=str(source),
                          source_size=source.stat().st_size,
                          source_mtime_ns=source.stat().st_mtime_ns,
                          role='PAPER_POSITIVE' if eps==3.5 else 'PAPER_KNOWN_LIMIT'))
    shutil.copy2(r8_study/'delegated/paper_baselines/PAPER_CONTRACTS.json',
                 ROOT/'data_public/classic_extension/ORIGINAL_PAPER_CONTRACTS.json')
    tasks = [dict(task_id=f'{case["case_id"]}__{arm}', case_id=case['case_id'], arm=arm)
             for case in cases for arm in ARMS]
    manifest = dict(status='FROZEN_BEFORE_NEW_PREDICTIONS_AND_OFFLINE_SCORING',
        UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), cases=cases, tasks=tasks,
        contract=dict(dataset_status='PUBLIC_PAPER_DEVELOPMENT; not blind confirmation',
            input='unmodified R8 public Y; no regeneration', forward_grid=[150,150],
            inverse_grid=[64,64], independent_real_material_cells=4096,
            chi_box=[0,20], initialization='air chi=0, APCS coefficients=0',
            source_count=20, receiver_count=40, domain_width_m=2., receiver_radius_m=3.75,
            frequency_Hz=400e6, noise='inherited fixed-Frobenius-norm10%; realization unchanged',
            dpcs_rank=15, rank_noise_supplement='NOT_RUN; main engineering prioritized',
            auxiliary_maxiter_per_stage=500, delta2D=[.01,.001,.001],
            transport='exact frequency-key coefficient zero-padding between nested masks',
            single_current_beta='material beta at each stage',
            dual_and_global_and_A_current_beta='fixed first-stage beta',
            ordinary_global_A_Aplus_current='same ordinary scalar GlobalDual PR/current operator',
            A_material='onehop material objective; unchanged analytic gradient/readout',
            Aplus_material='inherited real Mobius state-stationary candidates; prior MM acceptance',
            matched_six_prior=dict(smooth=1e-4,tv=1e-4),
            independent_paper_reference='no material prior; real-box adaptation; not author code',
            source_weights='data 1/||Y_s||²; state 1/||E_s||²; no GT',
            candidate_selection='minimum original physical Phi among initial + all paid stage endpoints',
            corrector=dict(algorithm='continuous L-BFGS-B physical objective', maxiter=300,
                maxfun=1200, maxls=30, ftol=1e-12, objective_scale=1e6,
                scaled_gtol=.01, checkpoint=80, direct_same_corrector=True),
            physical_backend=dict(operator='R8 equal-area disk scalar TM; padded FFT',
                GMRES_rtol=1e-9, GMRES_restart=100, GMRES_maxiter=10,
                checked_relative_residual=2e-9, fallback='same operator exact dense LU',
                engineering_change='cache immutable dense Green after first fallback; all arms identical'),
            stopping='record COEFFICIENT_TOL vs MAXITER per stage; physical CAP explicit',
            offline_metrics='paper epsilon-relative RMS + contrast relativeL2/amplitudeNRMSE/IoU',
            success_rule='NONE_ADDED; no overall binary success without pre-existing threshold',
            phase1_costs='FROZEN_113_INVERSIONS; phase2 separate ledger; no reuse double-charge'))
    dump(CONFIG, manifest)
    dump(OUT/'PREPARATION_RECEIPT.json', dict(status='COMPLETED', tasks=len(tasks),
        cpu_seconds=time.process_time()-PROCESS_CPU_STARTED,
        wall_seconds=time.perf_counter()-PROCESS_STARTED, physical_solves=0, GT_loaded=False))
    print(str(CONFIG), flush=True)


def solve(args):
    manifest = json.loads(CONFIG.read_text())
    task = next(t for t in manifest['tasks'] if t['task_id']==args.task_id)
    dest = OUT/task['task_id']
    if (dest/'RECEIPT.json').exists() or (dest/'FAILURE.json').exists():
        raise FileExistsError('Terminal task already exists; do not overwrite '+str(dest))
    from research.research_max_struct_r9.src.classic_extension import solve as run
    receipt = run(manifest, task, dest, PROCESS_STARTED, PROCESS_CPU_STARTED)
    print(json.dumps(dict(task_id=task['task_id'], status=receipt['status'],
        CPU=receipt['process_cpu_seconds'], Phi=receipt['original_Phi'],
        corrector=receipt['corrector'])), flush=True)


def queue(args):
    if not 1 <= args.lanes <= 2: raise ValueError('At most two one-thread physical workers')
    manifest = json.loads(CONFIG.read_text())
    folder = OUT/'queue'
    folder.mkdir(parents=True, exist_ok=True)
    if (folder/'QUEUE_STARTED.json').exists():
        raise FileExistsError('Queue already launched; keep costs and inspect before resuming')
    # Start one case on each lane; both controls continue through every arm.
    ordered = [t for arm in ARMS for t in manifest['tasks'] if t['arm']==arm]
    started = time.time()
    dump(folder/'QUEUE_STARTED.json', dict(status='STARTED', tasks=ordered, lanes=args.lanes,
        threads_per_child=1, UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        epoch_seconds=started, python=PYTHON, GT_access='NONE'))
    pending, active, rows = list(ordered), {}, []
    while pending or active:
        while pending and len(active)<args.lanes:
            task = pending.pop(0)
            dest = OUT/task['task_id']
            if (dest/'RECEIPT.json').exists() or (dest/'FAILURE.json').exists():
                raise FileExistsError('Unexpected terminal output '+str(dest))
            command = [PYTHON, str(Path(__file__).resolve()), 'solve', '--task-id', task['task_id']]
            log = open(folder/(task['task_id']+'.log'), 'w')
            wall = time.time()
            process = subprocess.Popen(command, cwd=REPO, stdout=log, stderr=subprocess.STDOUT)
            active[process.pid] = dict(process=process, task=task, log=log, epoch_start=wall,
                                       command=command)
            print('START', task['task_id'], process.pid, flush=True)
        for pid, item in list(active.items()):
            done, status, usage = os.wait4(pid, os.WNOHANG)
            if not done: continue
            item['process'].returncode = os.waitstatus_to_exitcode(status)
            item['log'].close()
            row = dict(task_id=item['task']['task_id'], command=item['command'],
                returncode=item['process'].returncode,
                epoch_start=item['epoch_start'], epoch_end=time.time(),
                child_process_CPU_seconds=usage.ru_utime+usage.ru_stime,
                peak_rss_bytes=usage.ru_maxrss,
                output_terminal='RECEIPT' if (OUT/item['task']['task_id']/'RECEIPT.json').exists()
                    else 'FAILURE' if (OUT/item['task']['task_id']/'FAILURE.json').exists() else 'LAUNCH_FAILURE')
            row['process_wall_seconds'] = row['epoch_end']-row['epoch_start']
            rows.append(row)
            dump(folder/'TASK_PROCESS_RECEIPTS.json', rows)
            print('FINISH', row['task_id'], row['returncode'], row['child_process_CPU_seconds'], flush=True)
            del active[pid]
        time.sleep(.25)
    dump(folder/'QUEUE_COMPLETE.json', dict(status='COMPLETED', tasks=len(rows),
        completed=sum(r['returncode']==0 for r in rows),
        failed=sum(r['returncode']!=0 for r in rows),
        UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        makespan_seconds=time.time()-started, epoch_start=started, epoch_end=time.time(),
        queue_process_CPU_seconds=time.process_time()-PROCESS_CPU_STARTED,
        actual_child_CPU_seconds=sum(r['child_process_CPU_seconds'] for r in rows),
        note='Actual elapsed is queue makespan, not sum of task walls; no reused R8 cost'))


def seal(args):
    manifest = json.loads(CONFIG.read_text())
    if OUT.name=='classic_extension_xinan':
        control=ROOT/'delegated/classic_extension/MAC_MIGRATION_CONTROL.json'
        record=json.loads(control.read_text())
        for child in record['active_children']:
            task_id=child['command'].split('--task-id ',1)[1].strip()
            local=ROOT/'RESULTS/classic_extension'/task_id
            if not (local/'RECEIPT.json').exists() and not (local/'FAILURE.json').exists():
                raise RuntimeError('Original Mac child still nonterminal: '+task_id)
    rows = []
    for task in manifest['tasks']:
        d = OUT/task['task_id']
        if (d/'RECEIPT.json').exists() and (d/'PREDICTION_SEAL.json').exists() and (d/'PREDICTION.npz').exists():
            receipt = json.loads((d/'RECEIPT.json').read_text())
            seal_record = json.loads((d/'PREDICTION_SEAL.json').read_text())
            if receipt['task_id']!=task['task_id'] or seal_record['task_id']!=task['task_id']:
                raise ValueError('Task identity mismatch')
            rows.append(dict(task_id=task['task_id'], status='SEALED'))
        elif (d/'FAILURE.json').exists():
            rows.append(dict(task_id=task['task_id'], status='FAILED'))
        else: raise RuntimeError('Not terminal/sealed: '+task['task_id'])
    result = dict(status='ALL_REQUESTED_PREDICTIONS_TERMINAL_BEFORE_GT', tasks=rows,
                  GT_loaded=False, UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
    dump(OUT/'ALL_PREDICTIONS_SEALED.json', result)
    return manifest, rows


def paid_prefix_selections(manifest):
    """Select only existing saved states by paid CPU and public Phi, before GT."""
    protocol=json.loads((ROOT/'configs/classic_extension_paid_budgets.json').read_text())
    process_costs={r['task_id']:r for r in json.loads(
        (OUT/'queue/TASK_PROCESS_RECEIPTS.json').read_text())}
    rows=[]
    for task in manifest['tasks']:
        if task['arm'] not in protocol['matched_arms']:continue
        d=OUT/task['task_id']
        if not (d/'RECEIPT.json').exists():
            rows.extend(dict(**task,CPU_budget_seconds=b,status='FAILED',GT_loaded=False)
                        for b in protocol['inclusive_process_CPU_budgets_seconds'])
            continue
        receipt=json.loads((d/'RECEIPT.json').read_text())
        total=process_costs[task['task_id']]['child_process_CPU_seconds']
        shift=max(0.,total-receipt['process_cpu_seconds'])
        curve=json.loads((d/'PHI_CPU_CURVE.json').read_text())
        initial=next(p for p in curve if p['phase']=='INITIAL')
        candidates=[dict(checkpoint='initial',file=None,key=None,
            source='actual prescribed all-zero initialization, recorded initial Phi',
            Phi=initial['Phi'],CPU_upper_bound=initial['process_cpu_seconds']+shift)]
        for st in json.loads((d/'STAGES.json').read_text()) if (d/'STAGES.json').exists() else []:
            candidates.append(dict(checkpoint=f"stage_{st['stage']:02d}",
                file=f"stage_{st['stage']:02d}.npz",key='m',Phi=st['Phi'],
                CPU_upper_bound=st['process_cpu_seconds']+shift))
        aux=next(p for p in curve if p['phase']=='CORRECTOR_START')
        candidates.append(dict(checkpoint='selected_auxiliary',file='AUXILIARY_PREDICTION.npz',
            key='m',Phi=aux['Phi'],CPU_upper_bound=aux['process_cpu_seconds']+shift))
        if (d/'CHECKPOINT_080.json').exists():
            checkpoint=json.loads((d/'CHECKPOINT_080.json').read_text())
            candidates.append(dict(checkpoint='corrector80',file='CHECKPOINT_080.npz',key='m',
                Phi=checkpoint['Phi'],CPU_upper_bound=checkpoint['process_cpu_seconds']+shift))
        candidates.append(dict(checkpoint='final',file='PREDICTION.npz',key='m',
            Phi=receipt['original_Phi'],CPU_upper_bound=total))
        for budget in protocol['inclusive_process_CPU_budgets_seconds']:
            eligible=[c for c in candidates if c['CPU_upper_bound']<=budget]
            if not eligible:
                rows.append(dict(**task,CPU_budget_seconds=budget,status='NO_SAVED_STATE_WITHIN_BUDGET',
                                 GT_loaded=False));continue
            chosen=min(eligible,key=lambda c:c['Phi'])
            latest=max(eligible,key=lambda c:c['CPU_upper_bound'])
            rows.append(dict(**task,status='PAID_SAVED_STATE_ONLY',CPU_budget_seconds=budget,
                actual_paid_prefix_CPU_seconds=min(total,float(budget)),total_process_CPU_seconds=total,
                selected_public_Phi_incumbent=chosen,latest_available_saved_checkpoint=latest,
                coverage_gap_seconds=min(total,float(budget))-latest['CPU_upper_bound'],
                conservative_CPU_shift_seconds=shift,CPU_shift_semantics='UPPER_SHIFT_including_output_tail',
                caveat='actual unsaved material incumbent between checkpoints is unknown; no interpolation',
                GT_loaded=False))
    dump(OUT/'PAID_PREFIX_PUBLIC_SELECTIONS_BEFORE_GT.json',dict(
        status='SELECTED_USING_SAVED_PUBLIC_PHI_AND_CPU_ONLY_BEFORE_GT',rows=rows,
        protocol='configs/classic_extension_paid_budgets.json',GT_loaded=False))
    return rows


def score(args):
    manifest, sealed = seal(args)  # No truth construction/import before this.
    if (OUT/'OFFLINE_SCORES.json').exists(): raise FileExistsError('Scoring output exists; retain it')
    prefix_rows=paid_prefix_selections(manifest)
    import numpy as np
    from research.research_max_struct_r8.src.paper_physics import acquisition, austria
    from research.research_max_struct_r8.src.paper_scores import material_scores
    rows = []
    truths = {}
    process_records=json.loads((OUT/'queue/TASK_PROCESS_RECEIPTS.json').read_text())
    process_costs={r['task_id']:r for r in process_records}
    for case in manifest['cases']:
        truth = austria(acquisition(64,400e6,dense=False)['centers'],case['epsilon_r'])
        truths[case['case_id']] = truth
        truth_file=ROOT/'sealed_truth/classic_extension'/(case['case_id']+'.npz')
        truth_file.parent.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(truth_file,m=truth)
    dump(ROOT/'sealed_truth/classic_extension/score_metadata.json',dict(
        status='PUBLIC_PAPER_GEOMETRY_CONSTRUCTED_ONLY_AFTER_ALL16_TERMINAL',
        inverse_grid=[64,64], source='NIE2016 Austria geometry; R8 austria function',
        first_access_UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
    for task, terminal in zip(manifest['tasks'], sealed):
        d = OUT/task['task_id']
        case = next(c for c in manifest['cases'] if c['case_id']==task['case_id'])
        if terminal['status']=='FAILED':
            rows.append(dict(**task, status='FAILED', receipt=json.loads((d/'FAILURE.json').read_text())))
            continue
        r = json.loads((d/'RECEIPT.json').read_text())
        p = np.load(d/'PREDICTION.npz')
        truth = truths[task['case_id']]
        def metrics(m):
            return dict(**material_scores(m,truth,case['epsilon_r']-1),
                paper_epsilon_relative_RMS=float(np.sqrt(np.mean(abs((m-truth)/(1+truth))**2))),
                paper_epsilon_relative_RMS_scatterer=float(np.sqrt(np.mean(
                    abs((m[truth>0]-truth[truth>0])/(1+truth[truth>0]))**2))))
        snapshots = dict(final=metrics(p['m']), auxiliary=metrics(p['m_auxiliary']))
        if (d/'CHECKPOINT_080.npz').exists():
            snapshots['corrector80'] = metrics(np.load(d/'CHECKPOINT_080.npz')['m'])
        for path in sorted(d.glob('stage_*.npz')):
            snapshots[path.stem] = metrics(np.load(path)['m'])
        rows.append(dict(**task, status='COMPLETED', **snapshots['final'],
            snapshots=snapshots, Phi=r['original_Phi'], relative_data_residual=r['relative_data_residual'],
            process_CPU_seconds=process_costs[task['task_id']]['child_process_CPU_seconds'],
            wall_seconds=process_costs[task['task_id']]['process_wall_seconds'],
            internal_process_timer_CPU_seconds=r['process_cpu_seconds'],
            cost_origin=process_costs[task['task_id']].get('external_CPU_semantics',
                'Mac wait4 inclusive child CPU; all imports/setup/output included'),
            physical_cost=r['physical_cost'], auxiliary_cost=r['auxiliary'], corrector=r['corrector'],
            stages=r['stages'],
            auxiliary_stage_CAP_count=sum(s.get('status')=='MAXITER' for s in r['stages']),
            all_auxiliary_stages_reached_paper_stop=(all(s.get('status')=='COEFFICIENT_TOL'
                for s in r['stages']) if r['stages'] else None),
            prior=r['material_prior'], public_development=True))
    dump(OUT/'OFFLINE_SCORES.json', rows)
    prefix_scores=[]
    for row in prefix_rows:
        if row['status']!='PAID_SAVED_STATE_ONLY':
            prefix_scores.append(row);continue
        case=next(c for c in manifest['cases'] if c['case_id']==row['case_id'])
        truth=truths[row['case_id']]
        scored={}
        for name in ['selected_public_Phi_incumbent','latest_available_saved_checkpoint']:
            point=row[name]
            material=(np.zeros(4096) if point['file'] is None else
                np.load(OUT/row['task_id']/point['file'])[point['key']])
            scored[name]=dict(**material_scores(material,truth,case['epsilon_r']-1),
                paper_epsilon_relative_RMS=float(np.sqrt(np.mean(abs((material-truth)/(1+truth))**2))))
        prefix_scores.append(dict(**row,offline_metrics=scored,
                                  GT_loaded_for_offline_scoring_after_public_selection=True))
    dump(OUT/'PAID_PREFIX_OFFLINE_SCORES.json',prefix_scores)
    if OUT.name=='classic_extension_xinan':
        original=ROOT/'RESULTS/classic_extension'
        external={x['task_id']:x for x in json.loads((original/'queue/TASK_PROCESS_RECEIPTS.json').read_text())}
        diagnostics=[]
        for task in manifest['tasks']:
            d=original/task['task_id']
            if not (d/'RECEIPT.json').exists():
                diagnostics.append(dict(**task,status='NOT_RUN_MAC_USER_HARDWARE_MIGRATION'));continue
            r=json.loads((d/'RECEIPT.json').read_text());p=np.load(d/'PREDICTION.npz')
            case=next(c for c in manifest['cases'] if c['case_id']==task['case_id'])
            truth=truths[task['case_id']]
            snapshots={}
            for name,material in [('final',p['m']),('auxiliary',p['m_auxiliary'])]:
                snapshots[name]=dict(**material_scores(material,truth,case['epsilon_r']-1),
                    paper_epsilon_relative_RMS=float(np.sqrt(np.mean(abs((material-truth)/(1+truth))**2))))
            if (d/'CHECKPOINT_080.npz').exists():
                material=np.load(d/'CHECKPOINT_080.npz')['m']
                snapshots['corrector80']=dict(**material_scores(material,truth,case['epsilon_r']-1),
                    paper_epsilon_relative_RMS=float(np.sqrt(np.mean(abs((material-truth)/(1+truth))**2))))
            e=external.get(task['task_id'])
            diagnostics.append(dict(**task,status='COMPLETED_PARTIAL_PLATFORM_DIAGNOSTIC',
                snapshots=snapshots,**snapshots['final'],Phi=r['original_Phi'],
                CPU_known_seconds=e['child_process_CPU_seconds'] if e else r['process_cpu_seconds'],
                CPU_semantics='wait4 inclusive' if e else 'internal process lowerbound; tail unknown',
                corrector=r['corrector'],scientific_comparison='SEPARATE_PARTIAL_MAC; do not pool with16Windows'))
        dump(original/'MAC_PARTIAL_OFFLINE_SCORES.json',dict(
            status='SCORED_ONLY_AFTER_ALL16_WINDOWS_AND_2_LEGACY_MAC_TERMINAL',
            rows=diagnostics,selection_or_budget_tuning=False))
    dump(OUT/'SCORING_RECEIPT.json', dict(status='COMPLETED', rows=len(rows),
        cpu_seconds=time.process_time()-PROCESS_CPU_STARTED,
        wall_seconds=time.perf_counter()-PROCESS_STARTED,
        GT_first_access='after all16 terminal seal in this invocation', physical_solves=0))
    print(json.dumps([{k:r.get(k) for k in ['case_id','arm','amplitude_normalized_nrmse',
        'otsu_iou','Phi','process_CPU_seconds']} for r in rows], indent=2))


def figures(args):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from research.research_max_struct_r8.src.paper_physics import acquisition, austria
    if not (OUT/'OFFLINE_SCORES.json').exists(): raise RuntimeError('Score after seal first')
    manifest = json.loads(CONFIG.read_text())
    folder = ROOT/'FIGURES/classic_extension'
    folder.mkdir(parents=True,exist_ok=True)
    for case in manifest['cases']:
        truth = austria(acquisition(64,400e6,dense=False)['centers'],case['epsilon_r'])
        panels = [('Truth', truth)]
        for arm in ARMS:
            file = OUT/f'{case["case_id"]}__{arm}'/'PREDICTION.npz'
            panels.append((arm,np.load(file)['m'] if file.exists() else None))
        # The declared box gives one full scale for every stage and every arm,
        # including auxiliary spikes larger than a final candidate's maximum.
        vmax = 20.
        for scale, top in [('full_range',vmax),('morphology',case['epsilon_r']-1)]:
            fig,axes = plt.subplots(3,3,figsize=(12,12), constrained_layout=True)
            for ax,(name,m) in zip(axes.flat,panels):
                if m is None:
                    ax.text(.5,.5,'FAILED',ha='center',va='center');ax.set_title(name)
                    continue
                im=ax.imshow(1+m.reshape(64,64),origin='lower',extent=(-1,1,-1,1),
                             vmin=1,vmax=1+top,cmap='viridis')
                ax.set_title(f'{name}\nmax εr={1+m.max():.3g}')
                ax.set_xlabel('x (m)');ax.set_ylabel('y (m)')
            fig.colorbar(im,ax=list(axes.flat),label='εr')
            fig.suptitle(f'{case["case_id"]}: {scale}; common scale 1–{1+top:.3g}\n'
                'Public development; morphology scale saturates amplitudes above GT')
            fig.savefig(folder/f'{case["case_id"]}_{scale}.png',dpi=150)
            plt.close(fig)
        fig,axes=plt.subplots(8,4,figsize=(12,24),constrained_layout=True)
        for i,arm in enumerate(ARMS):
            directory=OUT/f'{case["case_id"]}__{arm}'
            for j,name in enumerate(['stage_00','stage_01','stage_02','PREDICTION']):
                ax=axes[i,j];file=directory/(name+'.npz')
                if file.exists():
                    m=np.load(file)['m']
                    im=ax.imshow(1+m.reshape(64,64),origin='lower',extent=(-1,1,-1,1),
                                 vmin=1,vmax=1+vmax,cmap='viridis')
                    ax.set_title(f'{arm}: {name}\nmax εr={1+m.max():.3g}',fontsize=9)
                else:ax.text(.5,.5,'N/A',ha='center',va='center')
                ax.set_xticks([]);ax.set_yticks([])
        fig.colorbar(im,ax=list(axes.flat),label='εr; common full range')
        fig.suptitle(case['case_id']+' all auxiliary endpoints and final; full amplitude')
        fig.savefig(folder/f'{case["case_id"]}_stages_full_range.png',dpi=130)
        plt.close(fig)
        fig,ax=plt.subplots(figsize=(10,6),constrained_layout=True)
        for arm in ARMS:
            path=OUT/f'{case["case_id"]}__{arm}'/'PHI_CPU_CURVE.json'
            if not path.exists():continue
            points=json.loads(path.read_text())
            ax.plot([p['process_cpu_seconds'] for p in points],
                    [p['Phi'] for p in points],
                    '--' if arm.startswith('paper_') else '-',label=arm)
        ax.set_yscale('log');ax.set_xlabel('Paid process CPU seconds (imports/setup included)')
        ax.set_ylabel('Original physical Phi (reference arms omit material prior)')
        ax.set_title(case['case_id']+'; paid stages and continuous physical corrector')
        ax.legend(fontsize=8);ax.grid(alpha=.2)
        fig.savefig(folder/f'{case["case_id"]}_phi_cpu.png',dpi=150)
        plt.close(fig)
    score_rows=json.loads((OUT/'OFFLINE_SCORES.json').read_text())
    fig,axes=plt.subplots(1,3,figsize=(14,5),constrained_layout=True)
    columns=[('process_CPU_seconds','Process CPU (s)'),
             ('amplitude_normalized_nrmse','Amplitude NRMSE'),('otsu_iou','Otsu IoU')]
    for ax,(key,label) in zip(axes,columns):
        values=np.full((len(ARMS),2),np.nan)
        for row in score_rows:
            if row['status']=='COMPLETED':
                values[ARMS.index(row['arm']),[c['case_id'] for c in manifest['cases']].index(row['case_id'])]=row[key]
        im=ax.imshow(values,aspect='auto',cmap='viridis')
        for i in range(len(ARMS)):
            for j in range(2):
                ax.text(j,i,f'{values[i,j]:.3g}',ha='center',va='center',color='white',fontsize=9)
        ax.set_xticks([0,1],['εr=3.5','εr=4']);ax.set_yticks(range(len(ARMS)),ARMS)
        ax.set_title(label);fig.colorbar(im,ax=ax)
    fig.suptitle('Classic public development; common iteration ceilings, measured unequal CPU\n'
                 'Independent paper references omit material prior; no binary success threshold')
    fig.savefig(folder/'paired_cost_material_matrix.png',dpi=150)
    plt.close(fig)
    dump(OUT/'PLOT_RECEIPT.json',dict(status='COMPLETED',physical_solves=0,
        cpu_seconds=time.process_time()-PROCESS_CPU_STARTED,
        wall_seconds=time.perf_counter()-PROCESS_STARTED))


def ledger(args):
    manifest = json.loads(CONFIG.read_text())
    process_file = OUT/'queue/TASK_PROCESS_RECEIPTS.json'
    rows = json.loads(process_file.read_text()) if process_file.exists() else []
    costs = []
    counts = {}
    auxiliary_counts = {}
    for row in rows:
        d = OUT/row['task_id']
        p = d/'RECEIPT.json' if (d/'RECEIPT.json').exists() else d/'FAILURE.json'
        receipt = json.loads(p.read_text()) if p.exists() else {}
        costs.append(dict(**row, solver_receipt=str(p), physical_cost=receipt.get('physical_cost',{}),
                          auxiliary_cost=receipt.get('auxiliary',{}),
                          requested_source_RHS=20*receipt.get('physical_cost',{}).get('FULL_FORWARD',0),
                          completed_adjoint_RHS_requests=20*receipt.get('physical_cost',{}).get('DATA_ADJOINT',0),
                          RHS_semantics='actual_complex_RHS includes failed iterative work and fallback repeated RHS; requested counts are derived separately'))
        for key,value in receipt.get('physical_cost',{}).items():
            if any(t in key for t in ['maximum','bytes']): continue
            if isinstance(value,(int,float)): counts[key] = counts.get(key,0)+value
        for key,value in receipt.get('auxiliary',{}).get('operation_counts',{}).items():
            if any(t in key for t in ['peak','bytes','seconds']): continue
            if isinstance(value,(int,float)):
                auxiliary_counts[key]=auxiliary_counts.get(key,0)+value
    overhead = []
    for name in ['PREPARATION_RECEIPT.json','SCORING_RECEIPT.json','PLOT_RECEIPT.json']:
        path=OUT/name
        if path.exists():overhead.append(dict(path=str(path),**json.loads(path.read_text())))
    archived=ROOT/'delegated/classic_extension/pre_run_absolute_path_manifest/PREPARATION_RECEIPT.json'
    if archived.exists():
        overhead.append(dict(path=str(archived),**json.loads(archived.read_text())))
    original_prep=ROOT/'RESULTS/classic_extension/PREPARATION_RECEIPT.json'
    if OUT.name=='classic_extension_xinan' and original_prep.exists():
        overhead.append(dict(path=str(original_prep),**json.loads(original_prep.read_text())))
    q=OUT/'queue/QUEUE_COMPLETE.json'
    queue_receipt=json.loads(q.read_text()) if q.exists() else None
    actual=sum(r['child_process_CPU_seconds'] for r in rows)+sum(r['cpu_seconds'] for r in overhead)
    if queue_receipt:actual+=queue_receipt['queue_process_CPU_seconds']
    primary_CPU=actual
    audit=ROOT/'delegated/classic_extension'
    mac_file=audit/'MAC_PARTIAL_EXECUTION.json'
    mac_partial=json.loads(mac_file.read_text()) if OUT.name=='classic_extension_xinan' and mac_file.exists() else None
    if mac_partial:actual+=mac_partial['known_child_CPU_seconds_sum']
    all_counts=dict(counts)
    if mac_partial:
        for key,value in mac_partial['physical_counts'].items():
            all_counts[key]=all_counts.get(key,0)+value
    transport=[]
    if OUT.name=='classic_extension_xinan':
        for path in sorted(audit.glob('REMOTE_*_RECEIPT.json')):
            data=json.loads(path.read_text());data=data if isinstance(data,list) else [data]
            for record in data:
                transport.append(dict(source=str(path),**record))
                actual+=record.get('local_controller_CPU_seconds',0)+record.get('local_transport_child_CPU_seconds',0)+record.get('remote_helper_CPU_seconds',0)
    compiler_cpu=time.process_time()-PROCESS_CPU_STARTED
    actual+=compiler_cpu
    dump(OUT/'PHASE2_COST_LEDGER.json',dict(status='KNOWN_ACTUAL_CPU_LOWER_BOUND_WITH_EXPLICIT_GAPS',
        scope='classic2D Windows16 primary + Mac9 partial if present + transport/analysis; root-owned3D separate',
        expected_tasks=len(manifest['tasks']), measured_tasks=len(rows), tasks=costs,
        physical_counts=counts, auxiliary_operation_counts=auxiliary_counts,
        Green_transform_counts=dict(
            FFT_2D_source_transforms=counts.get('Green_FFT_2D_transforms',0)+auxiliary_counts.get('Green_FFT_2D_transforms',0),
            IFFT_2D_source_transforms=counts.get('Green_IFFT_2D_transforms',0)+auxiliary_counts.get('Green_IFFT_2D_transforms',0),
            note='physical+auxiliary nonoverlapping Green calls; coefficient transforms separate in auxiliary counters'),
        overhead=overhead,queue=queue_receipt,
        actual_new_known_CPU_seconds=actual, monetary_cost=None,
        primary_cohort_plus_shared_analysis_known_CPU_seconds=primary_CPU,
        Mac_partial_separate=mac_partial,transport_and_remote_helpers=transport,
        all_actual_physical_counts=all_counts,
        actual_inversion_attempts=len(rows)+(mac_partial['completed']+mac_partial['failed'] if mac_partial else 0),
        actual_inversions_completed=sum(r['returncode']==0 for r in rows)+(mac_partial['completed'] if mac_partial else 0),
        phase1_113_inversions='FROZEN; not included or overwritten',
        historical_public_Y_generation='REUSED_R8; zero new generation, historical fee not recharged',
        note='Auxiliary counters and physical counters separate; component CPU subset timers not added twice',
        compiler_CPU_seconds=compiler_cpu,
        unmeasured_costs=['interactive code authoring/inspection and py_compile CPU not instrumented; no physical solves',
            'Mac final2 worker startup/output tails and stopped dispatcher shutdown exactCPU unknown',
            'first remote mkdir/unpack/status helper CPU and shared SSH server overhead not instrumented',
            'engineering syntax and optional psutil probe failures have0physics but unknown CPU'],
        engineering_failures=(json.loads((audit/'MIGRATION_ENGINEERING_FAILURES.json').read_text())
                              if (audit/'MIGRATION_ENGINEERING_FAILURES.json').exists() else None)))
    print('PHASE2 actual known CPU',actual,'measured tasks',len(rows),flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('mode',choices=['freeze','solve','queue','seal','score','figures','ledger'])
    parser.add_argument('--task-id')
    parser.add_argument('--lanes',type=int,default=2)
    parser.add_argument('--r8-repo',help='R8 repository root for initial public-data packaging only')
    parser.add_argument('--output-root',type=Path,help='Separate platform cohort output, relative to repository')
    args=parser.parse_args()
    if args.output_root is not None:
        OUT=args.output_root if args.output_root.is_absolute() else REPO/args.output_root
    globals()[args.mode](args)
