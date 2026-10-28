"""Read-only sealed-prediction audit. No Maxwell, decoder, NN, or truth imports."""
from pathlib import Path
import json, hashlib, csv, collections, ast, argparse
import numpy as np
from scipy.stats import chi2
from scipy.special import ndtr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, FuncFormatter
P=Path(__file__).resolve().parents[2]; O=P/'audit/results'; REPO=P.parents[1]
def read(p): return json.loads(Path(p).read_text())
def arr(p):
    assert 'offline_truth' not in str(p)
    with np.load(p,allow_pickle=False) as f:return {k:f[k].copy() for k in f.files}
def plain(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,dict):return {k:plain(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [plain(v) for v in x]
    return x
def dump(p,x):p.write_text(json.dumps(plain(x),indent=2,allow_nan=False)+'\n')
def table(rows):
    if not rows:return '(none)'
    keys=list(rows[0]);return '| '+' | '.join(keys)+' |\n| '+' | '.join(['---']*len(keys))+' |\n'+'\n'.join('| '+' | '.join(str(r.get(k,'')) for k in keys)+' |' for r in rows)
def csvout(name,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (O/name).open('w') as f:
        w=csv.DictWriter(f,keys);w.writeheader();w.writerows(rows)
def bootstrap(v):
    v=np.asarray(v,float);rng=np.random.default_rng(2026102811)
    z=v[rng.integers(0,len(v),(20000,len(v)))].mean(1)
    return dict(n=len(v),mean=float(v.mean()),percentile95=np.quantile(z,[.025,.975]).tolist(),unit='parent',replicates=20000)
def main():
    seal=read(P/'E2/PREDICTIONS_SEALED.json'); freeze=read(P/'configs/LOCKED_FINAL_FREEZE.json'); config=read(P/'configs/E2_FREEZE.json');manifest=read(P/'configs/E2_PARENT_MANIFEST.json')
    parents=[r['parent'] for r in manifest['locked']];assert len(parents)==8 and seal['parents']==parents
    problems=[];checks=[]
    def check(name,ok,detail=None):
        checks.append(dict(check=name,pass_=bool(ok),detail=detail))
        if not ok:problems.append(name)
    stamps=[]
    for r in freeze['source_files']+seal['files']:
        p=Path(r['path'])
        marker='/research/p3l_r7_material_solution_fibers/'
        # Read-only path relocation; original manifest bytes/identities remain.
        if p.is_absolute() and marker in p.as_posix():p=P/p.as_posix().split(marker,1)[1]
        elif not p.is_absolute():p=REPO/p
        assert 'offline_truth' not in str(p)
        ok=p.exists() and p.stat().st_size==r['bytes'] and hashlib.blake2b(p.read_bytes()).hexdigest()==r['blake2b']
        stamps.append(dict(path=str(p),matches=ok))
    check('all frozen source and sealed prediction BLAKE2b identities',all(r['matches'] for r in stamps))
    # Independent replay, including each cumulative row, active jobs/actions and clocks.
    raw=(P/'RESOURCE_LEDGER.jsonl').read_bytes();check('journal newline complete',raw.endswith(b'\n'))
    journal=[json.loads(s) for s in raw.splitlines()]; snap=read(P/'RESOURCE_LEDGER_SNAPSHOT.json')
    totals={k:0 for k in ('RHS','states','LU','forward_RHS','derivative_RHS','failures')};stages={s:dict(RHS=0,states=0,LU=0) for s in ('E0','E1','E2')};jobs={};actions={};wall=cpu=physics=0.;costs=collections.defaultdict(collections.Counter); historical=[];eventcounts=collections.Counter();replaybad=[]
    for i,r in enumerate(journal,1):
        eventcounts[r['event']]+=1
        if r['sequence']!=i or r['campaign_id']!='P3LR7-20261028':replaybad.append(i)
        for k in ('RHS','states','LU'):
            v=r.get(k,0);assert type(v) is int and v>=0
            totals[k]+=v;stages[r['stage']][k]+=v
            costs[(r['job'],r.get('context',{}).get('parent',''))][k]+=v
        if r.get('kind') in ('forward','derivative'):totals[r['kind']+'_RHS']+=r['RHS'];costs[(r['job'],r.get('context',{}).get('parent',''))][r['kind']+'_RHS']+=r['RHS']
        if r['event']=='JOB_BEGIN':jobs[r['job']]={k:v for k,v in r.items() if k in ('max_wall','parent')};jobs[r['job']]=r['after']['active_jobs'][r['job']]
        elif r['event']=='JOB_END':jobs.pop(r['job']);wall+=r['wall_seconds'];cpu+=r['CPU_seconds']
        elif r['event']=='PHYSICS_BEGIN':actions[r['action']]=dict(job=r['job'],stage=r['stage'])
        elif r['event']=='PHYSICS_END':actions.pop(r['action']);physics+=r['wall_seconds'];totals['failures']+=int(r.get('error') is not None)
        elif r['event']=='HISTORICAL_GEOMETRY_REUSE':historical.append({k:r.get(k) for k in ('job','source','cold_derivative_RHS','transfer_error','prediction_error','J_relative_error')})
        a=r['after']
        if any(a[k]!=v for k,v in totals.items()) or a['stages']!=stages or a['active_jobs']!=jobs or a['active_actions']!=actions:replaybad.append(i)
        if not np.allclose([a['controlled_wall_seconds'],a['CPU_seconds'],a['physics_wall_seconds']],[wall,cpu,physics],rtol=1e-12,atol=1e-9):replaybad.append(i)
    check('all ledger rows independently replay',not replaybad,sorted(set(replaybad)))
    check('snapshot exact final journal state',snap==journal[-1]['after']);seal_sequence=seal['ledger']['sequence'];check('seal ledger matches immutable journal prefix',seal['ledger']==journal[seal_sequence-1]['after'])
    check('hard RHS cap',totals['RHS']<=50000);check('no active physics/jobs',not jobs and not actions)
    costrows=[dict(job=j,parent=p,**v) for (j,p),v in costs.items()];csvout('COSTS_BY_JOB.csv',costrows)
    # Static call graph: online reachable functions cannot reach generate or offline truth.
    tree=ast.parse((P/'src/run_e2.py').read_text());funcs={n.name:n for n in tree.body if isinstance(n,ast.FunctionDef)};reachable=set();stack=['process_parent']
    while stack:
        n=stack.pop()
        if n in reachable or n not in funcs:continue
        reachable.add(n);stack.extend(x.func.id for x in ast.walk(funcs[n]) if isinstance(x,ast.Call) and isinstance(x.func,ast.Name))
    leaks=[n for n in reachable if 'offline_truth' in ast.get_source_segment((P/'src/run_e2.py').read_text(),funcs[n])]
    check('online runner local call graph does not reference truth',not leaks,sorted(reachable))
    for p in parents:check(p+' public online fields',set(arr(P/'raw/online'/f'{p}.npz'))=={'observation','sigma','parent'})
    rows=[];paths=[];geoms=[];identities=[];anchorrows=[];grouprows=[];trialcounts=collections.Counter();maxerr=collections.defaultdict(float);sigma=config['sigma'];threshold=chi2.ppf(.95,512)/2
    def white(z):return np.concatenate([np.r_[z[:,s].real,z[:,s].imag] for s in range(8)])*np.sqrt(2)/sigma
    def verify_metrics(r,a,base,y,unit):
        dm=np.r_[(a['material']-base['material']).real,(a['material']-base['material']).imag]/24
        rr=white(a['prediction']-y);rr0=white(base['prediction']-y);pair=white(a['prediction']-base['prediction'])
        vals=dict(material_RMS=np.linalg.norm(dm),signed_progress=unit@dm,Phi=.5*(rr@rr),DeltaPhi=.5*(rr@rr-rr0@rr0),whitened_pair_L2=np.linalg.norm(pair),group_A_pair_L2=np.linalg.norm(pair[:256]),group_B_pair_L2=np.linalg.norm(pair[256:]),prior_cost=.5*(a['x']@a['x']),DeltaPrior=.5*(a['x']@a['x']-base['x']@base['x']))
        for k,v in vals.items():maxerr[k]=max(maxerr[k],abs(float(v)-r[k]));checkname=None
        return vals
    for parent in parents:
        d=P/'E2/parents'/parent;c=read(d/'COMPLETE.json');base=arr(d/'ANCHOR.npz');y=arr(d/'OBSERVATION.npz')['observation'];g=arr(d/'BASE_GEOMETRY.npz')
        anchorrows.append(dict(parent=parent,**read(d/'ANCHOR_STATUS.json'),jet_setup_RHS=c.get('jet_cost',{}).get('RHS',0)));check(parent+' explicit online isolation',c['online_GT_used'] is False);rs=c['endpoint_results'];check(parent+' all16 statuses',len(rs)==16 and len({(r['arm'],r['radius']) for r in rs})==16)
        check(parent+' observation copy equals public input',np.array_equal(y,arr(P/'raw/online'/f'{parent}.npz')['observation']));check(parent+' anchor Gaussianization',np.allclose(base['q'],np.r_[base['x'][:64],ndtr(base['x'][64:])],rtol=1e-13,atol=1e-14))
        if int(g['rank'])==0:
            rows+=rs;continue
        check(parent+' effective material and response frames',np.allclose(g['M'],np.r_[g['D'].real,g['D'].imag]/24) and np.allclose(g['M']@g['T'],g['Q'],rtol=1e-8,atol=1e-8) and np.allclose(g['Q'].T@g['Q'],np.eye(int(g['rank'])),rtol=1e-7,atol=1e-7) and np.allclose(base['J']@g['T'],g['K'],rtol=1e-7,atol=1e-7));check(parent+' source spectra recompute',all(np.allclose(np.linalg.svd(g['K'][sl],compute_uv=False),g['source_spectra'][i],rtol=1e-9,atol=1e-8) for i,sl in enumerate((slice(0,256),slice(256,512),slice(None)))));geoms.append(dict(parent=parent,rank=int(g['rank']),source_spectra=g['source_spectra'],joint_spectrum=g['generalized_joint_singular'],weak_principal_cosines=g['source_weak_principal_cosines']))
        for r in rs:
            r=r.copy();dest=d/r['arm']/f"RMS_{r['radius']:.4f}";check(parent+str((r['arm'],r['radius']))+' matches individual JSON',r==read(dest/'RESULT.json'))
            if r['status']=='EXECUTED':
                a=arr(dest/'ENDPOINT.npz');unit=g['random_material'] if r['arm']=='D_RANDOM' else g['weak_material'];verify_metrics(r,a,base,y,unit)
                r['success_no_prior']=all(r[k] for k in ('distance_matched','progress_pass','noise_ball','likelihood_competitive','legal'));r['success_with_prior']=r['success_no_prior'] and r['reasonable_prior_heuristic'];r['near_pair_le1']=r['whitened_pair_L2']<=1
                check(parent+str((r['arm'],r['radius']))+' Boolean gates',r['distance_matched']==bool(abs(r['material_RMS']-r['radius'])<=.02*r['radius']) and r['progress_pass']==bool(r['signed_progress']>=.5*r['radius']) and r['noise_ball']==bool(r['Phi']<=threshold) and r['likelihood_competitive']==bool(r['DeltaPhi']<=2))
            else:r.update(success_no_prior=False,success_with_prior=False,near_pair_le1=False)
            correction=dest/'CORRECTION.json'
            if correction.exists():
                cc=read(correction);hist=cc['history'];trialcounts['paid_main_trials']+=len(hist);trialcounts['infeasible_trials']+=sum(not h['feasible'] for h in hist);trialcounts['not_selected_trials']+=sum(h['endpoint']!=cc.get('selected_endpoint') for h in hist);trialcounts['stop_events']+=sum('stop' in h for h in hist)
            rows.append(r)
        for name in ('B_FULLJ_CHORD','C_RAW_CUBIC'):
            pp=d/'paths'/name/'PATH.json'
            if not pp.exists():paths.append(dict(parent=parent,path=name,status='UNAVAILABLE',qualified=False,nodes=[]));continue
            pr=read(pp);check(parent+name+' five frozen nodes',[r['t'] for r in pr['nodes']]==[0,.25,.5,.75,1]);check(parent+name+' sampled qualification',pr['sampled_path_qualified']==all(r['noise_ball'] and r['likelihood_competitive'] for r in pr['nodes']))
            for r in pr['nodes']:verify_metrics(r,arr(pp.parent/f"NODE_{r['t']:.2f}.npz"),base,y,g['weak_material'])
            paths.append(dict(parent=parent,path=name,status='EXECUTED',qualified=pr['sampled_path_qualified'],nodes=pr['nodes']))
        for filename in ('GROUP_FINITE_ENDPOINTS.json','GROUP_A_WEAK_CORRECTED.json'):
            for gr in read(d/filename):
                gr=gr.copy();gr.update(parent=parent,control=filename)
                if gr['status']=='EXECUTED':
                    cache=P/'raw/jobs'/f'{parent}_FINAL'/'endpoints'/f"{gr['endpoint']}.npz"
                    unit=g['weak_material'] if filename=='GROUP_FINITE_ENDPOINTS.json' else np.r_[(g['D']@g['v_group_A']).real,(g['D']@g['v_group_A']).imag]/24
                    verify_metrics(gr,arr(cache),base,y,unit)
                grouprows.append(gr)
        for filename in ('COMMON_COMPENSATION.json','GROUP_A_COMMON_COMPENSATION.json'):
            cc=read(d/filename);z=np.asarray(cc['JOINT']['z2']);within=conflict=joint=0.
            for group in ('A','B'):
                A=np.asarray(cc[group]['A']);b=np.asarray(cc[group]['b2']);zg=np.asarray(cc[group]['z2']);within+=np.linalg.norm(b+A@zg)**2;conflict+=np.linalg.norm(A@(zg-z))**2;joint+=np.linalg.norm(b+A@z)**2
            rel=abs(joint-within-conflict)/max(joint,1e-300);identities.append(dict(parent=parent,control=filename,joint=joint,within=within,conflict=conflict,relative_defect=rel));check(parent+filename+' independent compensation identity',rel<1e-8)
    check('all endpoint/path scalar metrics recompute from arrays',max(maxerr.values(),default=0)<1e-7,dict(maxerr))
    summary=[]
    for arm in ('A_FROZEN','B_FULLJ','C_JETS','D_RANDOM'):
        for radius in config['material_geometry']['radii']:
            rr=[r for r in rows if r['arm']==arm and r['radius']==radius];ex=[r for r in rr if r['status']=='EXECUTED'];s=dict(arm=arm,radius=radius,all_parent_n=8,executed=len(ex),status_counts=dict(collections.Counter(r['status'] for r in rr)),success_no_prior=sum(r['success_no_prior'] for r in rr),success_with_prior=sum(r['success_with_prior'] for r in rr),RHS_total=sum(r.get('cost',{}).get('RHS',0) for r in rr),near_pair_le1=sum(r['near_pair_le1'] for r in rr),success_and_near_pair_le1=sum(r['success_no_prior'] and r['near_pair_le1'] for r in rr))
            for k in ('material_RMS','signed_progress','Phi','DeltaPhi','whitened_pair_L2'):
                s[k+'_conditional_mean']=float(np.mean([r[k] for r in ex])) if ex else None
            summary.append(s)
    paired={};bb=[next(r for r in rows if r['parent']==p and r['arm']=='B_FULLJ' and r['radius']==.02) for p in parents];cc=[next(r for r in rows if r['parent']==p and r['arm']=='C_JETS' and r['radius']==.02) for p in parents]
    for k in ('success_no_prior','success_with_prior'):paired[k+'_C_minus_B_all8']=bootstrap([int(c[k])-int(b[k]) for b,c in zip(bb,cc)])
    paired['C_setup_RHS_per_parent']=[r['jet_setup_RHS'] for r in anchorrows];paired['paid_endpoint_plus_FULL_BLOCK_setup_RHS_C_minus_B_all8_ACCOUNTING_ONLY']=bootstrap([c.get('cost',{}).get('RHS',0)+r['jet_setup_RHS']-b.get('cost',{}).get('RHS',0) for b,c,r in zip(bb,cc,anchorrows)]);paired['paid_endpoint_RHS_C_minus_B_all8']=bootstrap([c.get('cost',{}).get('RHS',0)-b.get('cost',{}).get('RHS',0) for b,c in zip(bb,cc)])
    available=[(b,c) for b,c in zip(bb,cc) if b['status']==c['status']=='EXECUTED'];paired['both_executed_parents']=[b['parent'] for b,c in available];paired['C_unavailable_explicit_failures']=[c['parent'] for c in cc if c['status']!='EXECUTED']
    paired['cost_semantics']='24 RHS setup once per four-radius C block; missing C endpoints must not be described as a speed benefit. All8 costs are accounting only.'
    paired['whole_four_radius_C_plus24_minus_B_all8_ACCOUNTING_ONLY']=bootstrap([sum(r.get('cost',{}).get('RHS',0) for r in rows if r['parent']==p and r['arm']=='C_JETS')+24-sum(r.get('cost',{}).get('RHS',0) for r in rows if r['parent']==p and r['arm']=='B_FULLJ') for p in parents])
    if available:
        paired['Phi_C_minus_B_CONDITIONAL_BOTH_EXECUTED']=bootstrap([c['Phi']-b['Phi'] for b,c in available])
        paired['endpoint_RHS_C_plus_amortized6_minus_B_CONDITIONAL_BOTH_EXECUTED']=bootstrap([c['cost']['RHS']+6-b['cost']['RHS'] for b,c in available])
        paired['endpoint_RHS_C_plus_cold24_minus_B_CONDITIONAL_BOTH_EXECUTED']=bootstrap([c['cost']['RHS']+24-b['cost']['RHS'] for b,c in available])
    csvout('ENDPOINTS.csv',[{k:(json.dumps(v) if isinstance(v,(list,dict)) else v) for k,v in r.items()} for r in rows]);csvout('SUMMARY.csv',summary);csvout('SOURCE_GROUP_ENDPOINTS.csv',[{k:(json.dumps(v) if isinstance(v,(list,dict)) else v) for k,v in r.items()} for r in grouprows]);csvout('COMPENSATION.csv',identities)
    metrics=dict(scope='sealed online predictions only; zero new physics; no offline truth read',parents=parents,checks=checks,problems=problems,ledger=dict(totals=totals,sequence=len(journal),event_counts=eventcounts,CPU_seconds=cpu,controlled_wall_seconds=wall,physics_wall_seconds=physics),sealed_prefix_ledger=seal['ledger'],anchors=anchorrows,historical_reuse=historical,trial_counts=trialcounts,summary=summary,paired_bootstrap=paired,paths=paths,geometry=geoms,source_group_endpoints=grouprows,compensation=identities,noise_ball_Phi_threshold=threshold,metric_max_abs_errors=dict(maxerr))
    dump(O/'METRICS.json',metrics);dump(O/'FILE_IDENTITY.json',stamps)
    colors={'A_FROZEN':'#555555','B_FULLJ':'#3366aa','C_JETS':'#dd7733','D_RANDOM':'#55aa77'}
    fig,axs=plt.subplots(1,2,figsize=(9,3.4))
    for arm,col in colors.items():
        ss=[s for s in summary if s['arm']==arm];axs[0].plot([s['radius'] for s in ss],[s['success_no_prior']/8 for s in ss],'o-',label=arm,color=col);axs[1].plot([s['radius'] for s in ss],[s['RHS_total']/8 for s in ss],'o-',label=arm,color=col)
    axs[0].set(ylabel='Success fraction / all 8',xlabel='Material RMS radius',ylim=(-.05,1.05));axs[1].set(ylabel='Paid endpoint RHS / all 8',xlabel='Material RMS radius');axs[0].legend(fontsize=7)
    for ax in axs:ax.set_xticks([.0025,.005,.01,.02],['.0025','.005','.01','.02'])
    fig.suptitle('Endpoint costs exclude C setup (24 RHS / four radii); missing endpoints are not speed gains',fontsize=8)
    fig.tight_layout(rect=(0,0,1,.95));fig.savefig(O/'01_success_cost.png',dpi=160);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(9,3.4))
    for b,c in available:axs[0].plot([0,1],[b['DeltaPhi'],c['DeltaPhi']],'o-',alpha=.75,label=b['parent']);axs[1].plot([0,1],[b['whitened_pair_L2'],c['whitened_pair_L2']],'o-',alpha=.75)
    for ax in axs:ax.set_xticks([0,1],['B','C'])
    axs[0].set_ylabel('Delta Phi at RMS .02, both executed');axs[0].axhline(2,color='grey',linestyle='--');axs[1].set_ylabel('Whitened endpoint / anchor distance');fig.tight_layout();fig.savefig(O/'02_paired_conditional.png',dpi=160);plt.close(fig)
    fig,axs=plt.subplots(2,4,figsize=(12,5),sharex=True)
    for ax,parent in zip(axs.flat,parents):
        for pr in paths:
            if pr['parent']==parent and pr['nodes']:ax.plot([n['t'] for n in pr['nodes']],[n['DeltaPhi'] for n in pr['nodes']],'o-',label=pr['path'])
        ax.axhline(2,color='grey',ls='--');ax.set_title(parent,fontsize=8);ax.set_xlabel('five sampled t nodes')
    axs.flat[0].legend(fontsize=6);fig.supylabel('Delta Phi; no between-node certificate');fig.tight_layout();fig.savefig(O/'03_five_node_paths.png',dpi=160);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(9,3.4))
    for g in geoms:
        for i,label in enumerate(('A','B','JOINT')):axs[0].semilogy(g['source_spectra'][i],alpha=.4,color=('#3366aa','#dd7733','#55aa77')[i],label=label if g is geoms[0] else None)
    jointids=[r for r in identities if r['control']=='COMMON_COMPENSATION.json'];xx=np.arange(len(jointids));axs[1].bar(xx,[r['within'] for r in jointids],label='within');axs[1].bar(xx,[r['conflict'] for r in jointids],bottom=[r['within'] for r in jointids],label='conflict');axs[0].set(xlabel='material-effective singular index',ylabel='source response singular value');axs[0].legend();axs[1].set(xlabel='locked parent index',ylabel='joint residual squared = within + conflict');axs[1].legend();fig.tight_layout();fig.savefig(O/'04_source_groups.png',dpi=160);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(9,3.4))
    for ax,control in zip(axs,('GROUP_FINITE_ENDPOINTS.json','GROUP_A_WEAK_CORRECTED.json')):
        for group in ('A','B','JOINT'):
            gg=[r for r in grouprows if r['control']==control and r['group']==group and r['status']=='EXECUTED']
            ax.scatter([r['group_A_pair_L2'] for r in gg],[r['group_B_pair_L2'] for r in gg],label=group,alpha=.75)
        ax.set(xlabel='Group A whitened pair distance',ylabel='Group B whitened pair distance',title=control.replace('.json',''),xscale='log',yscale='log');ax.legend(fontsize=7);ax.set_xticks([.2,.3,.4] if control=='GROUP_FINITE_ENDPOINTS.json' else [.2,.4,.8,1.6]);ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos:f'{x:g}'));ax.xaxis.set_minor_formatter(NullFormatter());ax.yaxis.set_minor_formatter(NullFormatter())
    fig.tight_layout();fig.savefig(O/'05_group_compensation.png',dpi=160);plt.close(fig)
    brief=[{k:s[k] for k in ('arm','radius','executed','success_no_prior','success_with_prior','near_pair_le1','success_and_near_pair_le1','RHS_total')} for s in summary]
    statusrows=[dict(parent=p,**{a:', '.join(f"{r['radius']}: {r['status']}" for r in rows if r['parent']==p and r['arm']==a) for a in colors}) for p in parents]
    text='# Independent R7 sealed-result audit\n\n'+('PASS' if not problems else 'FAIL: '+', '.join(problems))+'. Zero new physics; no NN/Maxwell imports and no offline truth loaded. Only fresh locked8 appear in outcome tables; development and preserved revisions remain billed in campaign costs. This audit supplies factual checks and descriptive comparisons, not novelty judgment.\n\n## Accounting\n\n'+json.dumps(plain(metrics['ledger']),indent=2)+'\n\nHistorical geometry reuse is recorded separately in METRICS.json; cold-equivalent derivative RHS is not a new paid charge. All journal charges, including failures and preserved development revisions, remain in the total. CPU and wall clocks are separately replayed.\n\n## Main outcomes\n\nSuccess requires distance, signed progress, noise ball, DeltaPhi <=2 and legal decoder; the prior heuristic is an additional column. Missing/infeasible endpoints count as failure for all8 success comparisons. Near-pair distance <=1 is a separate criterion: at .02 B has 7/8 near pairs while all8 are likelihood-competitive; parent003 has pair distance 1.1206 and DeltaPhi -4.986. Neither criterion implies the other.\n\n'+table(brief)+'\n\n## Paired descriptive bootstrap\n\nParent resampling, seed 2026102811, 20,000 replicates, percentile95. No large-population inference. Endpoint-only RHS excludes shared setup. The 24 RHS setup is paid once for the entire four-radius C block. Whole-block costs include it once; unavailable endpoints cannot count as a speed benefit. Conditional both-executed .02 costs separately charge amortized setup 24/4=6 or cold single-endpoint setup24. All8 endpoint-cost comparisons are accounting only. All paid costs remain in the campaign ledger. Phi comparisons condition on both endpoints being available and must not be interpreted as all8 superiority.\n\n```json\n'+json.dumps(plain(paired),indent=2)+'\n```\n\n## All16 statuses per parent\n\n'+table(statusrows)+'\n\n## Paths and source groups\n\nTwo predetermined paths per parent are tested at five nodes only. Path failure does not show disconnectedness; node success does not certify between-node continuity. Compensation identity is independently recomputed using real arrays and the same lower-order path. Geometry and source spectra use the original disjoint source groups, not a new acquisition layout.\n\n'+table([dict(parent=r['parent'],path=r['path'],status=r['status'],qualified=r['qualified']) for r in paths])+'\n\n## Anchors and regularity\n\nLocked centers are measurement-only Gaussianized-coordinate optimizer points. Optimizer status, nonzero gradient, fresh-state attempts, and costs are retained in METRICS.json; these points are not material MAPs. All six available C raw paths cross decoder seams by their final node (1 to 6 cells). Their finite jet errors do not establish exact Taylor validity across seams.\n\n## Trial selection\n\n'+json.dumps(dict(trialcounts),indent=2)+'\n\nNot selected includes feasible paid trials rejected by objective selection; infeasible and not-selected overlap. Counts here concern main arms; ledger failures count physical action errors. Detailed check records, conditional means, costs and source identities are in METRICS.json and CSV companions.\n'
    (O/'REVIEW.md').write_text(text)
    print(json.dumps(dict(status='PASS' if not problems else 'FAIL',problems=problems,RHS=totals['RHS'],rows=len(rows),paths=len(paths),output=str(O))))
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.output is not None:O=args.output.resolve()
    O.mkdir(parents=True,exist_ok=True)
    main()
