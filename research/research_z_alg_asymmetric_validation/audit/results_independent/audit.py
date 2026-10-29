"""Refreshable saved-array/ledger audit. No physics imports, Model/Book, or state solves.

Only reads stored arrays, computes products/norms/SVD/determinants, and writes here.
"""
from pathlib import Path
import argparse,collections,csv,json,math,time
import numpy as np
audit_t0=time.perf_counter()
HERE=Path(__file__).resolve().parent
P=HERE.parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,default=HERE/'final_refresh',help='Audit output directory; only this directory is written')
OUTPUT=parser.parse_args().output.resolve()
OUTPUT.mkdir(parents=True,exist_ok=True)
def read(p):return json.loads(Path(p).read_text())
def arrays(p):
    with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def norm(x):return float(np.linalg.norm(x))
def relative(x,base):return norm(x)/max(norm(base),1e-300)
def rms(x):return norm(x)/math.sqrt(len(x))
def write(name,x):
    def default(a):
        if isinstance(a,np.generic):return a.item()
        if isinstance(a,np.ndarray):return a.tolist()
        raise TypeError(type(a).__name__)
    (OUTPUT/name).write_text(json.dumps(x,indent=2,ensure_ascii=False,default=default,allow_nan=False)+'\n')
def csvwrite(name,rows):
    if not rows:return
    cols=list(dict.fromkeys(k for row in rows for k in row))
    with (OUTPUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(rows)
def metrics(m,j,y,target,G,S,E0,mref,E=None):
    field=E0+G@j;L=np.eye(len(m))-m[:,None]*G;rhs=m[:,None]*E0
    pair=norm(y-target);sv=np.linalg.svd(L,compute_uv=False)
    out=dict(material_RMS=rms(m-mref),pair_relative=pair/max(norm(target),1e-300),pair_whitened=pair/.001,Phi=pair**2/(2*.001**2),
      state_relative=relative(L@j-rhs,rhs),sigma_min=float(sv[-1]),legal=bool(np.isrealobj(m) and np.all((m>=.5-1e-12)&(m<=4+1e-12))),
      receiver_field_consistency=relative(y-S@j,y),constitutive_relative=relative(j-m[:,None]*field,j))
    if E is not None:out['field_consistency']=relative(E-field,field)
    out['qualified_recomputed']=bool(out['pair_relative']<=1e-10 and out['state_relative']<=1e-10 and out['legal'] and out['sigma_min']>=1e-6)
    return out
ledger_text=(P/'COST_LEDGER.jsonl').read_text();ledger=[json.loads(s) for s in ledger_text.splitlines() if s.strip()]
issues=[]
def check(ok,where,detail):
    if not ok:issues.append(dict(where=where,detail=detail))
# Ledger counted once per BEGIN. Never add quoted result metadata to these totals.
runs={};begins={};ends={}
for r in ledger:
    run=r['run'];a=runs.setdefault(run,dict(run=run,actual_complex_RHS=0,closed_form_source_RHS_equivalent=0,LU=0,legacy_complex_factorizations=0,legacy_real_factorizations=0,legacy_complex_LU_subset=0,algebraic_real_factorizations=0,algebraic_complex_factorizations=0,algebraic_real_SVD=0,algebraic_real_RHS=0,algebraic_complex_RHS=0,BEGIN_count=0,zero_RHS_BEGIN_count=0,END_kernel_wall_seconds=0.,zero_RHS_END_kernel_wall_seconds=0.,current_state_kernel_seconds=0.,RUN_START_count=0,RUN_END_count=0))
    key=(run,r.get('action'))
    if r['event']=='RUN_START':a['RUN_START_count']+=1
    elif r['event']=='BEGIN':
        check(key not in begins,run,'duplicate BEGIN action');begins[key]=r;a['BEGIN_count']+=1
        for k in ['actual_complex_RHS','closed_form_source_RHS_equivalent','LU','algebraic_real_factorizations','algebraic_complex_factorizations','algebraic_real_SVD','algebraic_real_RHS','algebraic_complex_RHS']:a[k]+=r.get(k,0)
        if r['kind'].startswith('legacy_complex'):a['legacy_complex_LU_subset']+=r.get('LU',0)
        a['algebraic_real_RHS']+=r.get('real_algebra_RHS',0)
        if 'factorizations' in r:a['legacy_complex_factorizations' if r['kind'].startswith('legacy_complex') else 'legacy_real_factorizations']+=r['factorizations']
        if r.get('actual_complex_RHS',0)==0 and r.get('closed_form_source_RHS_equivalent',0)==0:a['zero_RHS_BEGIN_count']+=1
    elif r['event']=='END':
        check(key not in ends,run,'duplicate END action');ends[key]=r;a['END_kernel_wall_seconds']+=r['wall_seconds']
        b=begins.get(key,{})
        if b.get('actual_complex_RHS',0)==0 and b.get('closed_form_source_RHS_equivalent',0)==0:a['zero_RHS_END_kernel_wall_seconds']+=r['wall_seconds']
    elif r['event']=='CURRENT_ALGEBRA_STATE':a['current_state_kernel_seconds']+=r.get('kernel_seconds',0.)
    elif r['event']=='RUN_END':
        a['RUN_END_count']+=1;a['status']=r['status'];a['controlled_job_wall_seconds']=r['controlled_wall_seconds'];a['controlled_CPU_seconds']=r['controlled_CPU_seconds'];a['reported_run_combined']=r['run_combined_budget_used']
for name,r in runs.items():
    r['combined_RHS']=r['actual_complex_RHS']+r['closed_form_source_RHS_equivalent']
    check(r['RUN_START_count']==1,name,'run name not single-use')
    if r['RUN_END_count']:
        check(r['combined_RHS']==r['reported_run_combined'],name,'RUN_END budget differs from BEGIN sum')
        r['job_minus_END_kernel_seconds']=r['controlled_job_wall_seconds']-r['END_kernel_wall_seconds']
        r['job_minus_END_and_current_state_seconds']=r['job_minus_END_kernel_seconds']-r['current_state_kernel_seconds']
    r['unclosed_BEGIN_count']=sum(k not in ends for k in begins if k[0]==name)
metadata_notes=[]
orphans=[dict(run=r['run'],action=r['action'],path=r['path']) for r in ledger if r['event']=='SAVE' and (r['run'],r['action']) not in begins]
if orphans:metadata_notes.append(dict(code='SAVE_ACTION_REFERS_TO_END_SEQ',details=orphans,impact='Raw files remain readable and paid BEGIN RHS totals are unaffected; preserve ledger and add a separate mapping receipt if needed.'))
for name,r in runs.items():check(r['RUN_END_count']==1 and r['unclosed_BEGIN_count']==0,name,'unfinished run in ledger snapshot')
link_file=P/'audit/HIGH_PRECISION_SNAPSHOT_LINK_ADDENDUM.json'
link_addendum={'status':'NOT_PRESENT'}
if link_file.exists():
    link_addendum=read(link_file);events=[r for r in ledger if r['event']=='COMBINED_HIGH_PRECISION_SNAPSHOT_LINK_ADDENDUM']
    check(len(events)==1 and events[0]['links']==link_addendum['links'],'highprecision provenance','addendum file/event differs')
    for link in link_addendum['links']:
        paid=[begins.get(('INDEPENDENT_70DPS_FULL_MAXWELL',a),{}) for a in link['paid_source_actions']]
        check(len(paid)==2 and sum(r.get('actual_complex_RHS',0) for r in paid)==2 and (P/link['path']).exists(),'highprecision provenance','snapshot paid-action mapping invalid')
    check({(x['action'],x['path']) for x in orphans}=={(x['original_snapshot'],x['path']) for x in link_addendum['links']},'highprecision provenance','orphan mappings incomplete')
    if orphans:metadata_notes[0]['status']='EXPLICIT_APPEND_ONLY_LINK_ADDENDUM_PRESENT'
# Every E1 attempted target, including failed terminal radial targets.
endpoint_rows=[];pair_rows=[];context_rows=[];refs={};summaries={};raw_by_parent={}
for parent in range(3):
    for method in ['A','B']:
        root=P/'results/E1'/f'P{parent}'/method
        if not (root/'SUMMARY.json').exists():continue
        summary=read(root/'SUMMARY.json');summaries[parent,method]=summary;ref=arrays(root/'REFERENCE.npz');refs[parent,method]=ref
        ident=read(root/'CHART_IDENTITY.json')
        rrun=runs[f'E1_P{parent}_{method}']
        check(summary['total_actual_RHS']==rrun['actual_complex_RHS'],str(root),'summary RHS differs from ledger')
        check(sum(x['actual_RHS'] for x in summary['rows'])+summary['common_reference_RHS']==rrun['actual_complex_RHS'],str(root),'endpoint+common RHS do not sum')
        rr=[]
        for row in summary['rows']:
            file=root/row['key']/'ENDPOINT.npz';out=dict(parent=parent,method=method,key=row['key'],kind=row['kind'],claimed_qualified=row['qualified'],endpoint_exists=file.exists())
            check(read(root/row['key']/'RESULT.json')==row,str(file),'summary/RESULT differs')
            if file.exists():
                z=arrays(file);q=metrics(z['m'],z['j_full'],z['y_full'],z['target'],ref['G'],ref['S'],ref['E0'],ref['m'],z['E_full']);out.update(q)
                check(np.array_equal(z['target'],ref['y']),str(file),'target changed')
                check(np.allclose(z['free_values'],row['free_values'],rtol=0,atol=1e-14),str(file),'free target differs from row')
                check(np.allclose(z['m'][ref['free']],z['free_values'],rtol=0,atol=1e-14),str(file),'material free coordinate not prescribed')
                for k in ['material_RMS','pair_relative','pair_whitened','Phi','state_relative','sigma_min']:
                    check(abs(q[k]-row[k])<=1e-12*max(1,abs(row[k])),str(file),'reported '+k+' differs from arrays')
                qualified=q['qualified_recomputed']
                if method=='B':
                    jc,Ec=z['j_chart'],z['E_chart'];den=np.sum(abs(Ec)**2,axis=1);valid=den>1e-20
                    shared=relative(jc-z['m'][:,None]*Ec,jc);null=relative(ref['S']@jc-z['target'],z['target'])
                    recovered=np.sum(Ec[valid].conj()*jc[valid],axis=1)/den[valid]
                    out.update(shared_material_relative=shared,receiver_null_relative=null,chart_field_consistency=relative(Ec-ref['E0']-ref['G']@jc,Ec),minimum_field_squared=float(min(den)),node_count=int(sum(~valid)),recovered_real_material_error=float(max(abs(recovered-z['m'][valid]))) if np.any(valid) else None,chart_vs_full_current_relative=relative(jc-z['j_full'],z['j_full']))
                    check(np.flatnonzero(~valid).tolist()==row['node_cells'],str(file),'node status differs')
                    for k in ['shared_material_relative','receiver_null_relative','minimum_field_squared']:
                        check(abs(out[k]-row[k])<=1e-12*max(1,abs(row[k])),str(file),'current report '+k+' differs')
                    qualified=qualified and shared<=1e-10 and null<=1e-10
                check(bool(qualified)==row['qualified'],str(file),'qualified flag differs from saved arrays')
                for k in ['receiver_field_consistency','constitutive_relative','field_consistency']:
                    check(q[k]<1e-10,str(file),'saved full state fails '+k)
            else:check(not row['qualified'],str(file),'qualified but no ENDPOINT')
            endpoint_rows.append(out);rr.append(out)
        context_rows.append(dict(parent=parent,method=method,attempts=len(rr),qualified=sum(x['claimed_qualified'] for x in rr),grid_qualified=sum(x['claimed_qualified'] for x in rr if x['kind']=='fixed_grid'),radial_qualified=sum(x['claimed_qualified'] for x in rr if x['kind']=='radial_continuation'),arrays=sum(x['endpoint_exists'] for x in rr),actual_RHS=rrun['actual_complex_RHS'],LU=rrun['LU'],controlled_job_wall_seconds=rrun['controlled_job_wall_seconds'],logged_END_kernel_seconds=rrun['END_kernel_wall_seconds'],zero_RHS_END_kernel_seconds=rrun['zero_RHS_END_kernel_wall_seconds'],current_state_kernel_seconds=rrun['current_state_kernel_seconds'],endpoint_wall_subtotal=sum(x['controlled_wall_seconds'] for x in summary['rows'])))
    if (parent,'A') not in summaries or (parent,'B') not in summaries:continue
    A,B=summaries[parent,'A'],summaries[parent,'B'];ra,rb=refs[parent,'A'],refs[parent,'B']
    check(all(np.array_equal(ra[k],rb[k]) for k in ['G','S','E0','m','y','dep','free']),f'parent{parent}','A/B reference or chart identity changed')
    ar={x['key']:x for x in A['rows']};br={x['key']:x for x in B['rows']}
    check(ar.keys()==br.keys(),f'parent{parent}','A/B attempted target key sets differ')
    for key in sorted(ar.keys()|br.keys()):
        a,b=ar.get(key),br.get(key)
        if a is None or b is None:continue
        root=P/'results/E1'/f'P{parent}';both=a['qualified'] and b['qualified']
        out=dict(parent=parent,key=key,kind=a['kind'],A_qualified=a['qualified'],B_qualified=b['qualified'],same_free_target=bool(np.allclose(a['free_values'],b['free_values'],rtol=0,atol=1e-14)),both_qualified=both,A_RHS=a['actual_RHS'],B_RHS=b['actual_RHS'])
        check(out['same_free_target'],f'parent{parent}/{key}','A/B free targets differ')
        check(a['qualified']==b['qualified'],f'parent{parent}/{key}','A/B coverage differs')
        if both:
            za,zb=arrays(root/'A'/key/'ENDPOINT.npz'),arrays(root/'B'/key/'ENDPOINT.npz');out['material_pair_RMS']=rms(za['m']-zb['m']);out['response_pair_relative']=relative(za['y_full']-zb['y_full'],za['target'])
        pair_rows.append(out)
# E2 arithmetic audit uses the declared E1 first-four coefficient submodel; no solve.
e2={"status":"NOT_PRESENT"};e2rows=[]
f=P/'results/E2/SUMMARY.json'
if f.exists():
    s=read(f);ref=refs[0,'A'];idx=np.array([0,1,2,3]);G=ref['G'][np.ix_(idx,idx)];S=ref['S'][:,idx];E0=ref['E0'][idx];mref=np.array(read(P/'EXPERIMENT_CONTRACT.json')['E2']['public_reference_material'])
    for row in s['starts']:
        z=arrays(P/'results/E2'/f"START_{row['start']:03d}.npz");q=metrics(z['m'],z['j'],z['y'],z['target'],G,S,E0,mref,z['E']);check(q['qualified_recomputed']==row['qualified'],f"E2 start {row['start']}",'qualification differs');e2rows.append(dict(kind='multistart',index=row['start'],**q))
    reused=[]
    for cy in s['cycles']:
        for row in cy['nodes']:
            z=arrays(P/'results/E2'/f"CYCLE_{cy['cycle']}_{row['step']:02d}.npz");q=metrics(z['m'],z['j'],z['y'],z['target'],G,S,E0,mref);check(q['qualified_recomputed']==row['qualified'],f"E2 cycle {cy['cycle']} step{row['step']}",'qualification differs');e2rows.append(dict(kind='cycle',cycle=cy['cycle'],index=row['step'],**q))
            if row['reused_paid_first_node']:reused.append(dict(cycle=cy['cycle'],step=row['step'],new_RHS=row['new_RHS']))
    first=runs['E2_ASYMMETRIC_DISCRETE_BRANCHES']['actual_complex_RHS'];recovery=runs['E2_CYCLE_REPORTING_RECOVERY']['actual_complex_RHS'];check(first==s['paid_first_run_RHS'] and recovery==s['paid_recovery_RHS'] and first+recovery==s['total_paid_RHS'],'E2','reporting repair paid totals mismatch')
    check(reused==[dict(cycle=0,step=1,new_RHS=0)],'E2','first paid node not reused exactly once at zero additional RHS')
    e2=dict(status=s['status'],starts=len(s['starts']),qualified_starts=sum(r['qualified'] for r in s['starts']),found_root_count=s['root_count_found'],cycle_nodes=sum(len(c['nodes']) for c in s['cycles']),completed_cycles=sum(c['completed'] for c in s['cycles']),cycle_return_RMS=[c['return_RMS'] for c in s['cycles']],paid_first_run_RHS=first,paid_recovery_RHS=recovery,total_paid_RHS=first+recovery,reused_nodes=reused,uniqueness='NOT_ESTABLISHED',nontrivial_monodromy='NOT_RUN; out-and-back paths',independent_numpy_validation='33 multistart endpoints yes; cycles fresh LU forward plus saved-state audit, no independent numpy re-solve receipt')
# Certificates: static argument check and numerical terminal-output consistency only.
certificates=[]
cf=P/'audit/INTERVAL_CHART_CERTIFICATES.json'
if cf.exists():
 for c in read(cf):
    t=c['attempts'][-1];m=np.array(c['material_center']);w=np.zeros(len(m));w[c['dependent']]=t['dependent_radius'];w[c['free']]=t['halfwidth'];legal_margin=float(min(np.min(m-w-.5),np.min(4-m-w)))
    ok=t['status']=='CERTIFIED_FROZEN_COEFFICIENT_LOCAL_GRAPH' and max(t['mapping_inclusion_ratios'])<1 and max(t['contraction_ratios'])<1 and t['Neumann_norm']<1 and legal_margin>0
    check(ok,f"interval parent{c['parent']}",'terminal certificate conditions fail')
    certificates.append(dict(parent=c['parent'],halfwidth=t['halfwidth'],free=c['free'],dependent=c['dependent'],max_inclusion_ratio=max(t['mapping_inclusion_ratios']),max_contraction_ratio=max(t['contraction_ratios']),Neumann_norm=t['Neumann_norm'],legal_box_margin=legal_margin,status='STATIC_ARGUMENT_AND_TERMINAL_OUTPUT_CONSISTENT' if ok else 'FAIL',scope='Exact stored binary64 coefficients and target only; no interval replay, Hankel transcendental, PDE, full fiber, or global probability certificate'))
# E4 local saved node densities/integration. Products/SVD/determinants, no tangent solve.
e4={"status":"NOT_PRESENT"};e4rows=[]
f=P/'results/E4_local/SUMMARY.json'
if f.exists():
    s=read(f);ref=refs[0,'A'];dep,free=ref['dep'],ref['free'];orders={}
    prior=lambda m:float(np.prod(2*m/(16-.25)))
    for order,result in s['orders'].items():
        am=[];bm=[];wa=[];wb=[];ev=[]
        for row in result['rows']:
            z=arrays(P/'results/E4_local'/f"NODE_{order}_{row['index']:03d}.npz");da=prior(z['m_A'])/abs(np.linalg.det(z['J_A'][:,dep]));area=float(np.sqrt(np.linalg.det(z['tangent_B'].T@z['tangent_B'])));normal=float(np.prod(np.linalg.svd(z['J_B'],compute_uv=False)));db=prior(z['m_B'])*area/normal
            qa=metrics(z['m_A'],z['j_A_full'],z['y_A'],z['target'],ref['G'],ref['S'],ref['E0'],ref['m']);qb=metrics(z['m_B'],z['j_B_full'],z['y_B'],z['target'],ref['G'],ref['S'],ref['E0'],ref['m']);check(qa['qualified_recomputed'] and qb['qualified_recomputed'],f"E4 {order}/{row['index']}",'stored node physical state fails')
            ec=ref['E0']+ref['G']@z['j_B_chart'];shared=relative(z['j_B_chart']-z['m_B'][:,None]*ec,z['j_B_chart']);check(shared<1e-10,f"E4 {order}/{row['index']}",'shared material fails')
            check(abs(da/row['direct_density']-1)<1e-12 and abs(db/row['current_density']-1)<1e-12,f"E4 {order}/{row['index']}",'saved density differs')
            e4rows.append(dict(order=int(order),index=row['index'],material_pair_RMS=rms(z['m_A']-z['m_B']),direct_density=da,current_density=db,density_relative_gap=abs(da-db)/da,A_pair_relative=qa['pair_relative'],B_pair_relative=qb['pair_relative'],B_shared_material_relative=shared))
            am.append(z['m_A']);bm.append(z['m_B']);wa.append(row['integration_weight']*da);wb.append(row['integration_weight']*db);ev.append(float(row['offset'][0]>0))
        am,bm,wa,wb,ev=map(np.asarray,[am,bm,wa,wb,ev]);out={}
        for method,mats,weights in [('A',am,wa),('B',bm,wb)]:
            val=dict(coarea_tile_density_mass=float(sum(weights)),mean=(weights@mats/sum(weights)).tolist(),second_moment=(weights@(mats*mats)/sum(weights)).tolist(),event_probability=float(weights@ev/sum(weights)))
            out[method]=val
            check(abs(val['coarea_tile_density_mass']/result[method]['coarea_tile_density_mass']-1)<1e-12 and np.max(abs(np.array(val['mean'])-result[method]['mean']))<1e-12 and abs(val['event_probability']-result[method]['event_probability'])<1e-12,f'E4 order{order}/{method}','quadrature summary differs from nodes')
        out['nodes']=len(am);orders[order]=out
    e4=dict(status=s['status'],orders=orders,convergence_reported=s['convergence'],cross_chart_reported=s['cross_chart'],paid=runs['E4_LOCAL_CONDITIONAL_MATERIAL_INTEGRATION']['actual_complex_RHS'],planning_allocation=read(P/'EXPERIMENT_CONTRACT.json')['planning_allocation_combined']['E4_local'],global_missing_mass=s['global_missing_mass'],original_noisy_posterior=s['original_noisy_posterior'],interpretation='Local numerical quadrature reference only; no certified global normalization/missing mass')
legacy={"status":"NOT_PRESENT"}
f=P/'audit/E4_legacy_replay/REPLAY_STATUS.json'
if f.exists():
    s=read(f);r=runs['legacy_E4'];legacy=dict(gate=s['gate'],actual_complex_RHS=r['actual_complex_RHS'],closed_form_source_RHS_equivalent=r['closed_form_source_RHS_equivalent'],new_levels=s['new_levels'],reused_levels=s['reused_levels'],scope=s['certification'])
    check(r['actual_complex_RHS']==s['fresh_actual_complex_RHS'] and r['closed_form_source_RHS_equivalent']==s['fresh_closed_form_source_RHS_equivalent'],'legacy E4','ledger and receipt differs')
    saved=read(P/'audit/E4_legacy_replay/e4_two_cell_posterior.json');frozen=read(P/'inputs/Research_Z_ALG_Research_Kit/zalgo_probability/e4_two_cell_posterior.json')
    check(saved['results']['64']==frozen['results']['64'],'legacy E4','cached64 differs from immutable input')
    keys=['scaled_evidence','mean_x1','mean_x2','event_x1_gt_x2','psi_x1_squared']
    def difference(a,b):return {k:abs(a[k]-b[k])/(abs(b[k]) if k=='scaled_evidence' else 1) for k in keys}
    recomputed={m:difference(saved['results']['32'][m],saved['results']['64'][m]) for m in ['original_box','ordered_quotient']}
    cross=difference(saved['results']['32']['original_box'],saved['results']['32']['ordered_quotient'])
    check(recomputed==s['new32_vs_cached64'] and cross==s['cross_new32_direct_vs_quotient'],'legacy E4','saved comparison arithmetic differs')
    legacy.update(cached64_exactly_reused=True,new32_vs_cached64_recomputed=recomputed,new32_direct_vs_quotient_recomputed=cross)

# Final equal-trace benchmark: all saved state/J/current trial arrays, no solves.
cost={"status":"NOT_PRESENT"};costrows=[]
f=P/'results/COST_CONTROL/SUMMARY.json'
if f.exists():
    s=read(f);physical_states=0;Js=0;algebra_states=0;maxstate=0.;maxtransfer=0.;maxalg=0.;cached={};buffered_actions=[]
    for file in sorted((P/'results/COST_CONTROL').glob('*.npz')):
        z=arrays(file);groups=collections.defaultdict(dict)
        for k,v in z.items():label,key=k.split('_',1);groups[label][key]=v
        cached[file.name]=groups
        buffered_actions.extend(int(label[1:]) for label in groups if not label.startswith('salg'))
        parent=int(file.name.split('_')[1][1:]) if file.name.startswith('COMMON') else int(file.name.split('_')[0][1:]);ref=refs[parent,'A']
        for label,a in groups.items():
            if 'j' in a and 'y' in a:
                physical_states+=1;L=np.eye(6)-a['m'][:,None]*ref['G'];rhs=a['m'][:,None]*ref['E0'];err=relative(L@a['j']-rhs,rhs);maxstate=max(maxstate,err)
                check(err<=1e-10 and relative(a['E']-ref['E0']-ref['G']@a['j'],a['E'])<=1e-10 and relative(a['y']-ref['S']@a['j'],a['y'])<=1e-10,str(file),'buffered physical state inconsistent')
            elif 'receiver_transfer' in a:
                Js+=1;L=np.eye(6)-a['m'][:,None]*ref['G'];err=relative(a['receiver_transfer']@L-ref['S'],ref['S']);maxtransfer=max(maxtransfer,err)
                check(err<1e-10,str(file),'buffered transfer fails transpose linear relation')
            elif 'r' in a and 'j' in a:
                algebra_states+=1;field=ref['E0']+ref['G']@a['j'];rscale=max(norm(ref['E0'])/math.sqrt(12),1e-6);rc=(a['j']-a['m'][:,None]*a['E']).T;rr=np.r_[rc.real.ravel(),rc.imag.ravel()]/rscale;err=norm(rr-a['r']);maxalg=max(maxalg,err)
                check(relative(a['E']-field,field)<1e-10 and err<1e-10,str(file),'buffered current algebra state/residual differs')
    expected_buffered=[r['action'] for r in ledger if r['run']=='OPTIMIZED_EQUAL_TRACE_COST_CONTROL' and r['event']=='PAID_ARRAYS_BUFFERED']
    check(sorted(buffered_actions)==sorted(expected_buffered) and len(buffered_actions)==len(set(buffered_actions)),'cost trace','paid buffered snapshots missing or duplicated')
    for row in s['rows']:
        p=row['parent'];ref=refs[p,'A'];groups=cached[row['trace']];found=[]
        for label,a in groups.items():
            if label.startswith('salg'):continue
            try:begin=begins.get(('OPTIMIZED_EQUAL_TRACE_COST_CONTROL',int(label[1:])),{})
            except ValueError:continue
            if begin.get('label')=='QA_COMMON_FULL_VALIDATION':found.append(a)
        check(len(found)==1,row['trace'],'not exactly one independent endpoint state')
        if not found:continue
        a=found[0];q=metrics(a['m'],a['j'],a['y'],ref['y'],ref['G'],ref['S'],ref['E0'],ref['m'],a['E'])
        check(q['qualified_recomputed']==row['qualified'] and np.array_equal(a['m'],np.array(row['material'])),row['trace'],'cost endpoint report differs')
        check(abs(row['total_seconds']-row['solver_seconds']-row['validation_seconds']-row['persistence_seconds'])<1e-12,row['trace'],'method window arithmetic differs')
        costrows.append(dict(parent=p,repeat=row['repeat'],target=row['target'],method=row['method'],qualified=row['qualified'],RHS=row['RHS'],solver_seconds=row['solver_seconds'],validation_seconds=row['validation_seconds'],persistence_seconds=row['persistence_seconds'],setup_seconds=row['setup_seconds'],method_window_total=row['total_seconds']+row['setup_seconds'],trace=row['trace']))
    grouped=[]
    for old in s['grouped']:
        rows=[r for r in costrows if r['parent']==old['parent'] and r['method']==old['method']];v=dict(parent=old['parent'],method=old['method'],cases=len(rows),qualified=sum(r['qualified'] for r in rows),RHS=sum(r['RHS'] for r in rows),total_seconds=sum(r['method_window_total'] for r in rows),solver_plus_validation_seconds=sum(r['solver_seconds']+r['validation_seconds']+r['setup_seconds'] for r in rows));grouped.append(v)
        check(all(abs(v[k]-old[k])<=1e-12 for k in ['cases','qualified','RHS','total_seconds','solver_plus_validation_seconds']),'cost grouped','summary arithmetic differs')
    maxpair=0.
    bytarget=collections.defaultdict(list)
    for r in s['rows']:bytarget[r['parent'],r['repeat'],r['target']].append(r)
    for key,rows in bytarget.items():
        check(len(rows)==3,key,'missing method for target/repeat')
        for a in rows:
            for b in rows:maxpair=max(maxpair,rms(np.array(a['material'])-b['material']))
    pooled={m:{k:sum(r[k] for r in grouped if r['method']==m) for k in ['total_seconds','solver_plus_validation_seconds','RHS','cases','qualified']} for m in ['DYNAMIC_FULLJ','FROZEN_J','CURRENT_FULL_KERNEL']}
    saving=1-pooled['CURRENT_FULL_KERNEL']['total_seconds']/pooled['DYNAMIC_FULLJ']['total_seconds'];kernel_saving=1-pooled['CURRENT_FULL_KERNEL']['solver_plus_validation_seconds']/pooled['DYNAMIC_FULLJ']['solver_plus_validation_seconds']
    check(sum(r['RHS']+r['setup_RHS'] for r in s['rows'])+s['shared_reference_RHS']==runs['OPTIMIZED_EQUAL_TRACE_COST_CONTROL']['actual_complex_RHS'],'cost control','paid row+shared RHS mismatch')
    cost=dict(status='ALL_243_SAVED_CASES_QUALIFIED' if len(costrows)==243 and all(r['qualified'] for r in costrows) else 'FAIL',pooled=pooled,grouped=grouped,observed_method_window_saving_current_vs_dynamic=saving,solver_plus_validation_saving=kernel_saving,total_cost_gate_20percent='PASS' if saving>=.20 else 'FAIL',physical_buffered_states=physical_states,buffered_material_J=Js,buffered_current_algebra_states=algebra_states,max_state_relative=maxstate,max_transfer_equation_relative=maxtransfer,max_current_residual_absolute_difference=maxalg,max_material_pair_RMS=maxpair,controlled_whole_phase_wall_seconds=runs['OPTIMIZED_EQUAL_TRACE_COST_CONTROL']['controlled_job_wall_seconds'],shared_reference_RHS=s['shared_reference_RHS'],actual_RHS=runs['OPTIMIZED_EQUAL_TRACE_COST_CONTROL']['actual_complex_RHS'],timing_scope='Pooled method-window totals include setup, solver, validation, batch persistence. Whole phase RUN_END additionally includes shared references and checkpoint/report overhead; not assigned to individual methods. Three deterministic parents, three timing repeats each; not nine independent parents.')
# Posthoc rank diagnostics: recover transfer by scalar source consistency and verify TL=S.
rank={"status":"NOT_PRESENT"};rankrows=[]
f=P/'audit/ENDPOINT_RANK_REVIEW.json'
if f.exists():
    s=read(f);saved=[r for r in ledger if r['event']=='SAVE' and r['run']=='POSTHOC_SAVED_ENDPOINT_RANK_REVIEW'];check(len(saved)==len(s['rows'])==738,'posthoc rank','row/save count mismatch');maxtransfer=0.;maxsource=0.
    for row,save in zip(s['rows'],saved):
        ref=refs[row['parent'],row['method']];z=arrays(P/save['path']);original=arrays(P/'results/E1'/f"P{row['parent']}"/row['method']/row['key']/'ENDPOINT.npz');J=z['J'];sv=np.linalg.svd(J,compute_uv=False);sd=np.linalg.svd(J[:,ref['dep']],compute_uv=False)
        check(np.array_equal(z['m'],original['m']),'posthoc rank','saved material/row ordering differs');C=J[0::2]+1j*J[1::2];E=original['E_full'].T;den=np.sum(abs(E)**2,axis=0);T=np.sum(E.conj()*C,axis=0)/den
        transfer=relative(T[None,:]@(np.eye(6)-z['m'][:,None]*ref['G'])-ref['S'],ref['S']);source=relative(C-T[None,:]*E,C);maxtransfer=max(maxtransfer,transfer);maxsource=max(maxsource,source)
        check(transfer<1e-10 and source<1e-10,'posthoc rank','saved J inconsistent with cached fields/transpose equation')
        ranks={str(c):int(sum(sv>sv[0]*c)) for c in [1e-8,1e-10,1e-12]};condition=float(sd[0]/sd[-1]);check(ranks==row['numerical_rank'] and abs(sv[-1]-row['response_sigma_min'])<1e-12 and abs(condition-row['dependent_condition'])<1e-10,'posthoc rank','singular/rank summary differs')
        rankrows.append(dict(parent=row['parent'],method=row['method'],key=row['key'],original_qualified=row['original_qualified'],rank_1e8=ranks['1e-08'],rank_1e10=ranks['1e-10'],rank_1e12=ranks['1e-12'],response_sigma_min=float(sv[-1]),dependent_condition=condition,transfer_equation_relative=transfer,source_shared_transfer_relative=source))
    check(s['new_adjoint_RHS']==runs['POSTHOC_SAVED_ENDPOINT_RANK_REVIEW']['actual_complex_RHS']==738,'posthoc rank','adjoint paid totals mismatch')
    rank=dict(status=s['status'],stored_nodes=len(rankrows),new_paid_adjoint_RHS=738,all_rank4=all(r['rank_1e8']==r['rank_1e10']==r['rank_1e12']==4 for r in rankrows),minimum_response_sigma_min=min(r['response_sigma_min'] for r in rankrows),maximum_dependent_condition=max(r['dependent_condition'] for r in rankrows),maximum_transfer_equation_relative=maxtransfer,maximum_source_consistency_relative=maxsource,scope='POSTHOC at 738 saved legal endpoint arrays, including failures. Not original online stopping and no between-node certificate.')
# Read replay traces; do not execute CurrentChart or any algebra/state solve.
replay={"status":"NOT_PRESENT"};replayrows=[]
f=P/'raw/E1_CURRENT_TRACE_REPLAY/REPRODUCIBILITY_RECEIPT.json'
if f.exists():
    s=read(f);N=arrays(P/'audit/PROBE_ARRAYS.npz')['N'];maxm=0.;maxz=0.;maxdecode=0.;maxres=0.;maxD=0.;calls=0
    for r in s['receipts']:
        ref=refs[r['parent'],'B'];a=arrays(P/r['output']);original_path=P/'results/E1'/f"P{r['parent']}"/'B'/r['key']/'ENDPOINT.npz';dm=dz=None
        if original_path.exists():
            original=arrays(original_path);dm=float(np.max(abs(a['final_m']-original['m'])));dz=float(np.max(abs(a['final_z']-original['z'])));maxm=max(maxm,dm);maxz=max(maxz,dz);check(dm==r['comparison']['m_max_absolute_difference'] and dz==r['comparison']['z_max_absolute_difference'] and dm<=1e-12 and dz<=1e-12,r['output'],'replay final/original mismatch')
        count=len(a.get('z',[]));check(count==r['trace_calls'],r['output'],'replay call count differs');calls+=count
        ns,q=2,N.shape[1];cut=q*ns;scale=max(norm(ref['j'])/math.sqrt(12),1e-6);rscale=max(norm(ref['E0'])/math.sqrt(12),1e-6)
        z0=a['reference_z'];c0=(z0[:cut].reshape(ns,q)+1j*z0[cut:2*cut].reshape(ns,q)).T;j0=ref['j']-scale*N@c0;e0=ref['E0']+ref['G']@j0
        for i in range(count):
            z=a['z'][i];c=(z[:cut].reshape(ns,q)+1j*z[cut:2*cut].reshape(ns,q)).T;j=j0+scale*N@c;E=e0+scale*(ref['G']@N)@c;m=np.empty(6);m[ref['dep']]=z[2*cut:];m[ref['free']]=a['free'][i]
            decode=max(norm(j-a['j'][i]),norm(E-a['E'][i]),norm(m-a['m'][i]));rc=(a['j'][i]-m[:,None]*a['E'][i]).T;rr=np.r_[rc.real.ravel(),rc.imag.ravel()]/rscale;res=norm(rr-a['r'][i]);maxdecode=max(maxdecode,decode);maxres=max(maxres,res)
            check(decode<1e-10 and res<1e-10,r['output'],'replay decoded state or residual mismatch')
            if a['D_present'][i]:
                K=scale*(N-m[:,None]*(ref['G']@N));D=np.zeros((24,24))
                for source in range(2):
                    rows=slice(source*6,(source+1)*6);imrows=slice(12+source*6,12+(source+1)*6);cr=slice(source*q,(source+1)*q);ci=slice(cut+source*q,cut+(source+1)*q)
                    D[rows,cr]=K.real;D[rows,ci]=-K.imag;D[imrows,cr]=K.imag;D[imrows,ci]=K.real
                    for h,node in enumerate(ref['dep']):D[source*6+node,2*cut+h]=-a['E'][i,node,source].real;D[12+source*6+node,2*cut+h]=-a['E'][i,node,source].imag
                gap=norm(D/rscale-a['D'][i]);maxD=max(maxD,gap);check(gap<1e-10,r['output'],'replay saved current Jacobian differs')

        replayrows.append(dict(parent=r['parent'],key=r['key'],trace_calls=count,original_endpoint_present=original_path.exists(),m_max_absolute_difference=dm,z_max_absolute_difference=dz))
    check(runs['E1_CURRENT_TRACE_REPLAY']['actual_complex_RHS']==0 and runs['E1_CURRENT_TRACE_REPLAY']['closed_form_source_RHS_equivalent']==0,'replay','nonzero physical RHS charged')
    replay=dict(status='SAVED_TRACE_ARRAYS_CONSISTENT',targets=len(replayrows),current_algebra_calls=calls,existing_endpoint_comparisons=sum(r['original_endpoint_present'] for r in replayrows),maximum_m_absolute_difference=maxm,maximum_z_absolute_difference=maxz,maximum_decoded_state_absolute_difference=maxdecode,maximum_residual_absolute_difference=maxres,maximum_current_J_absolute_difference=maxD,new_physics_RHS=0,scope='Stored pure-algebra replay; original qualified flags/original z determine original warm-start policy. No new success selection, physical validation, or timing advantage.')
# Closed campaign/source guard checked statically; no Book instance.
closure={"status":"NOT_PRESENT"}
if (P/'CLOSURE.json').exists():
    closure=read(P/'CLOSURE.json');resource=read(P/'RESOURCE_SUMMARY.json');source=(P/'src/accounting.py').read_text()
    closure['source_hardcap30000_and_perrun10000_present']=bool('self.actual+self.formula+rhs+formula>30000 or self.run_used+rhs+formula>10000' in source)
    closure['source_guard_present']=bool("if (P/'CLOSURE.json').exists():raise BudgetStop" in source and source.index("if (P/'CLOSURE.json').exists():raise BudgetStop")<source.index("self.append(dict(event='RUN_START'"))
    check(closure['source_hardcap30000_and_perrun10000_present'],'closure','frozen source hard-cap checks differ')
    check(closure['source_guard_present'],'closure','closed-campaign guard missing/below first append')
    check(closure['actual_RHS']==sum(r['actual_complex_RHS'] for r in runs.values()) and closure['formula_equivalent']==sum(r['closed_form_source_RHS_equivalent'] for r in runs.values()),'closure','totals differ from BEGIN sums')
    check(resource['LU_factorizations']==sum(r['LU'] for r in runs.values()) and resource['conservative_combined']==closure['combined'],'closure','resource summary differs')
    check(closure['no_active_jobs'] and all(r['RUN_END_count']==1 and r['unclosed_BEGIN_count']==0 for r in runs.values()),'closure','active/unclosed run remains')

# Comparison by deterministic context, no inferential CI.
comparisons=[]
for p in range(3):
    a=next(x for x in context_rows if x['parent']==p and x['method']=='A');b=next(x for x in context_rows if x['parent']==p and x['method']=='B');pairs=[x for x in pair_rows if x['parent']==p and x['both_qualified']]
    comparisons.append(dict(parent=p,A_attempts=a['attempts'],B_attempts=b['attempts'],A_qualified=a['qualified'],B_qualified=b['qualified'],grid_qualified_each=a['grid_qualified'],paired_qualified=len(pairs),max_material_pair_RMS=max(x['material_pair_RMS'] for x in pairs),A_RHS=a['actual_RHS'],B_RHS=b['actual_RHS'],RHS_reduction=1-b['actual_RHS']/a['actual_RHS'],A_job_wall=a['controlled_job_wall_seconds'],B_job_wall=b['controlled_job_wall_seconds'],observed_controlled_wall_reduction=1-b['controlled_job_wall_seconds']/a['controlled_job_wall_seconds']))
totals={k:sum(r[k] for r in runs.values()) for k in ['actual_complex_RHS','closed_form_source_RHS_equivalent','LU','legacy_complex_LU_subset','legacy_complex_factorizations','legacy_real_factorizations','algebraic_real_factorizations','algebraic_complex_factorizations','algebraic_real_SVD','algebraic_real_RHS','algebraic_complex_RHS']};totals['combined_RHS']=totals['actual_complex_RHS']+totals['closed_form_source_RHS_equivalent'];totals['nonlegacy_LU']=totals['LU']-totals['legacy_complex_LU_subset'];totals['failed_runs_actual_complex_RHS']=sum(r['actual_complex_RHS'] for r in runs.values() if r.get('status','').startswith('FAILED'));totals['controlled_job_wall_seconds']=sum(r.get('controlled_job_wall_seconds',0) for r in runs.values());totals['zero_RHS_END_kernel_wall_seconds']=sum(r['zero_RHS_END_kernel_wall_seconds'] for r in runs.values())
check(totals['combined_RHS']<=30000,'campaign','combined hard cap exceeded')
check(all(r['combined_RHS']<=10000 for r in runs.values()),'campaign','per-run hard cap exceeded')
result=dict(status='PASS_SAVED_ARRAY_AUDIT' if not issues else 'FAIL_SAVED_ARRAY_AUDIT',scope='Read saved arrays and ledger only; zero new physics, zero new state solves, no Model/Book instantiation',ledger_snapshot=dict(lines=len(ledger),last_seq=ledger[-1]['seq'],last_run=ledger[-1]['run']),totals=totals,run_count=len(runs),E1_contexts=comparisons,E1_arrays=sum(x['endpoint_exists'] for x in endpoint_rows),E1_qualified=sum(x['claimed_qualified'] for x in endpoint_rows),E1_attempts=len(endpoint_rows),E1_same_grid_and_radial_targets=all(x['same_free_target'] for x in pair_rows),E2=e2,interval_certificates=certificates,E4_local=e4,E4_legacy=legacy,cost_control=cost,posthoc_rank=rank,current_trace_replay=replay,closure=closure,high_precision_link_addendum=link_addendum,issues=issues,metadata_notes=metadata_notes,protocol_scope_notes=['Atlas online radial rank stopping was not explicitly gated; later paid posthoc J diagnostics cover saved endpoint arrays only, not between-node regularity or retroactive online compliance.','E4 integration halfwidth0.05 exceeds interval-certified halfwidth0.00625 for parent0; quadrature convergence is numerical and does not certify the whole integration graph tile.'],audit_wall_seconds=time.perf_counter()-audit_t0,statistical_scope='Three selected deterministic contexts; no population confidence interval or repeated-run timing inference',timing_scope='RUN_END controlled job wall includes saved artifacts/bookkeeping within run; END operation wall is only a subtotal and can overlap for nested legacy operations. Current algebra state kernel_seconds is separate. No overall efficiency GO judgment made here.')
write('AUDIT.json',result);csvwrite('LEDGER_BY_RUN.csv',list(runs.values()));csvwrite('E1_CONTEXTS.csv',context_rows);csvwrite('E1_ENDPOINTS.csv',endpoint_rows);csvwrite('E1_PAIRED_TARGETS.csv',pair_rows);csvwrite('E2_ENDPOINTS.csv',e2rows);csvwrite('E4_LOCAL_NODES.csv',e4rows);csvwrite('COST_CONTROL_CASES.csv',costrows);csvwrite('POSTHOC_RANK.csv',rankrows);csvwrite('CURRENT_REPLAY.csv',replayrows)
print(json.dumps(dict(status=result['status'],totals=totals,issues=issues,E1_contexts=comparisons,E2={k:v for k,v in e2.items() if k!='independent_numpy_validation'},E4_local_status=e4['status']),ensure_ascii=False,indent=2))
