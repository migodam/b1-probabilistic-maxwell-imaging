"""Prepared fresh postscore energy diagnostic. Never run before root unlock.

Only24 optional GT-current forwards are allowed. All inverse/final/stage
physical currents must be reused from saved files or reported missing.
"""
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[_k]='1'
from pathlib import Path
import argparse,json,sys,time,traceback
from collections import defaultdict
import numpy as np
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from research.research_max_struct_r10.experiments.score_twod_comparison import (
    public_gate,readjson,loadnpz,resolve,stamp,dump,csvwrite,algorithm,METHODS,LABELS)
from research.research_max_struct_r10.src.energy_diagnostics import Geometry,physical_current,energy,ref

def material_valid(m):return m.shape==(576,) and not np.iscomplexobj(m) and np.isfinite(m).all()
def current_valid(j):return j is not None and j.shape==(576,12) and np.isfinite(j).all()

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--authorize-postseal',required=True);ap.add_argument('--authorization-note',required=True)
    ap.add_argument('--baseline-manifest',type=Path,default=ROOT/'configs/twod_tasks.json');ap.add_argument('--candidate-manifest',type=Path,default=ROOT/'configs/denominator_confirmation_tasks.json')
    ap.add_argument('--baseline-queue',type=Path,default=ROOT/'delegated/twod/REMOTE_QUEUE_RECEIPT.json');ap.add_argument('--candidate-queue',type=Path,default=ROOT/'delegated/twod/DENOMINATOR_QUEUE_RECEIPT.json')
    ap.add_argument('--freeze',type=Path,default=ROOT/'delegated/twod/FREEZE_SCORING.json');ap.add_argument('--truth-dir',type=Path,default=ROOT/'sealed_truth/twod')
    ap.add_argument('--scored-dir',type=Path,default=ROOT/'delegated/twod_score/scored');ap.add_argument('--output',type=Path,default=ROOT/'delegated/twod_energy/scored')
    ap.add_argument('--gt-current-cache',type=Path,default=ROOT/'delegated/twod_energy/gt_current_cache')
    ap.add_argument('--figure-dir',type=Path,default=ROOT/'FIGURES/twod_comparison');args=ap.parse_args()
    start=time.perf_counter();cpu=time.process_time();freeze,tasks,cases,entries,evidence,_=public_gate(args)
    score_receipt=readjson(args.scored_dir/'RECEIPT.json')
    if score_receipt['status'] not in ('COMPLETE','COMPLETE_WITH_SCORING_ERRORS') or score_receipt['expected_tasks']!=168:
        raise PermissionError('Fresh material scoring must complete before energy diagnostics')
    if score_receipt['truth_access']!='POST_ALL168_TERMINAL_ROOT_AUTHORIZED_ONLY':raise PermissionError('Missing valid material-scoring unlock')
    scored_seal=readjson(args.scored_dir/'ALL168_PUBLIC_SEAL.json')
    for s in scored_seal['evidence']:
        p=Path(s['path']);now=stamp(p)
        if (now['bytes'],now['mtime_ns'])!=(s['bytes'],s['mtime_ns']):raise RuntimeError('Material-scoring public seal changed: '+str(p))
    score_rows=readjson(args.scored_dir/'SCORES_FINAL.json');score_index={r['task_id']:r for r in score_rows}
    if len(score_index)!=168:raise ValueError('Material scores do not cover all168 terminal tasks')
    if args.output.exists():raise FileExistsError('Write-once fresh energy output exists')
    if any((args.figure_dir/name).exists() for name in ('ENERGY_FIXED_TAU_POINT1.png','ORACLE_CAPACITY_BY_CASE.png')):
        raise FileExistsError('Fresh energy figure outputs already exist; no GT-forward work started')
    # Public cache: one operator archive and one observation archive per case,
    # never reread/re-factor these operators for each of seven methods.
    filecache={};pubs={};sources=evidence+[stamp(args.scored_dir/'RECEIPT.json'),stamp(__file__)];requests={};joint_paths={}
    for e in entries:
        task=e['task'];case=task['case_id'];pub={}
        for pp in task['public_files']:
            p=resolve(pp)
            if p not in filecache:filecache[p]=loadnpz(p);sources.append(stamp(p))
            pub.update(filecache[p])
        if case not in pubs:pubs[case]=pub
        else:
            for key in ('G','S','E0','Y'):
                if not np.array_equal(pubs[case][key],pub[key]):raise RuntimeError('Public operator/observation mismatch within case')
        if not e['sealed']:continue
        d=e['directory'];paths=sorted(d.glob('stage_*_joint.npz'))
        if paths:joint_paths[task['task_id']]=paths[-1]
        rank=None
        if paths and paths[-1].with_suffix('.json').exists():rank=readjson(paths[-1].with_suffix('.json')).get('receiver_rank')
        if rank is None and (d/'R10_DPCS_RULE.json').exists():rank=readjson(d/'R10_DPCS_RULE.json').get('rank')
        if rank is not None and (d/'NUMERICAL_SCHEDULE.json').exists():
            schedule=readjson(d/'NUMERICAL_SCHEDULE.json')
            for stage in schedule['stages']:requests[case,int(rank),int(stage['mf'])]=None
    first=pubs[cases[0]];G,S,E=first['G'],first['S'],first['E0']
    if G.shape!=(576,576) or S.shape[1]!=576 or E.shape!=(576,12):raise RuntimeError('Declared24 LU/288 RHS budget applies only to frozen576x12 physics')
    for pub in pubs.values():
        if not all(np.array_equal(pub[k],first[k]) for k in ('G','S','E0')):raise RuntimeError('Expected one exact frozen G/S/E operator group')
    out=args.output;out.mkdir(parents=True)
    dump(out/'ENERGY_PUBLIC_GATE.json',dict(status='ALL168_TERMINAL_AND_MATERIAL_SCORING_COMPLETE',root_authorization_note=args.authorization_note,
        prior_material_scoring=str(args.scored_dir/'RECEIPT.json'),truth_opened_by_this_process=False,evidence=evidence))
    # FIRST fresh GT access in this process, after both public and score gates.
    truth={};jtruth={};GT_LU=0;GT_RHS=0;gt_cost=[];args.gt_current_cache.mkdir(parents=True,exist_ok=True)
    operator_stamps=[stamp(p) for p,data in filecache.items() if 'G' in data or 'E0' in data]
    for case in cases:
        p=args.truth_dir/(case+'.npz');z=loadnpz(p);sources.append(stamp(p));m=z['m'].ravel()
        if not material_valid(m):raise ValueError('Invalid sealed GT material schema')
        truth[case]=m;j=z.get('j',z.get('J'));cp=args.gt_current_cache/(case+'.npz');mp=cp.with_suffix('.json');used_cache=False
        if cp.exists() or mp.exists():
            if not cp.exists() or not mp.exists():raise RuntimeError('Incomplete offline GT-current cache; preserve it, do not silently rerun a forward')
            cached=loadnpz(cp);cm=readjson(mp)
            if cm['truth_source']!=stamp(p) or cm['operator_sources']!=operator_stamps or not np.array_equal(cached['m_gt'],m) or not current_valid(cached['j_gt']):
                raise RuntimeError('GT-current cache provenance mismatch; no automatic replacement/forward')
            j=cached['j_gt'];used_cache=True
        if current_valid(j):
            jtruth[case]=j;gt_cost.append(dict(case_id=case,source=str(cp if used_cache else p),mode='REUSED_OFFLINE_GT_CURRENT_CACHE' if used_cache else'REUSED_SAVED_GT_CURRENT',LU=0,complex_RHS=0))
        else:
            if GT_LU>=24 or GT_RHS+12>288:raise RuntimeError('Frozen offline GT-forward budget exceeded')
            t,c=time.perf_counter(),time.process_time();jtruth[case]=physical_current(G,E,m);GT_LU+=1;GT_RHS+=12
            gt_cost.append(dict(case_id=case,source=str(p),mode='OFFLINE_SCORING_GT_CURRENT_ONLY',LU=1,LU_dimension=576,complex_RHS=12,
                CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t))
        if not used_cache:
            np.savez_compressed(cp,m_gt=m,j_gt=jtruth[case])
            dump(mp,dict(status='OFFLINE_SCORING_ONLY_NOT_INVERSE_INPUT',truth_source=stamp(p),operator_sources=operator_stamps,
                creation_forward_LU=gt_cost[-1]['LU'],creation_complex_RHS=gt_cost[-1]['complex_RHS'],source_mode=gt_cost[-1]['mode']))
    geom_start,geom_cpu=time.perf_counter(),time.process_time();geom=Geometry(G,S)
    geometry_cost=dict(CPU_seconds=time.process_time()-geom_cpu,wall_seconds=time.perf_counter()-geom_start,operator_groups=1,
        scope='One complete geometry setup shared across all cases/arms; all fixed tau .03/.1/.3, real control, both ordered projections/intersections/remainder')
    # Operator-only coefficient-range SVDs live in geom.affines[(rank,shape,mf)].
    # No arm ID enters this cache. JD and oracle values additionally depend on case.
    oracle=[];jd_cache={};oracle_cache={}
    for case,rank,mf in sorted(requests):
        t,c=time.perf_counter(),time.process_time();Y=pubs[case]['Y'];JD,U,q=geom.affine(Y,rank,(24,24),mf)
        JD=jd_cache.setdefault((case,rank),JD);off=jtruth[case]-JD;projected=q@(q.conj().T@off);irreducible=off-projected
        row=dict(case_id=case,regime=next(r['regime'] for r in score_rows if r['case_id']==case),receiver_rank=rank,mf=mf,
            fourier_columns=U.shape[1],actual_coefficient_range_rank=q.shape[1],ambient_current_dimension=576,
            affine_GT_irreducible_squared_energy=energy(irreducible),GT_current_squared_energy=energy(jtruth[case]),
            affine_GT_irreducible_fraction_of_GT_current_energy=energy(irreducible)/energy(jtruth[case]),
            physical_final_current_is_affine_constrained=False,truth_access='OFFLINE_SCORING_ONLY',
            CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t)
        oracle.append(row);oracle_cache[case,rank,mf]=(row,JD,q,JD+projected)
    single=[];two=[];remainders=[];metrics=[];gaps=[]
    def decompose(base,quantity,x,family='complex'):
        aa,bb,cc=geom.decompose(x,family)
        for rr,target in ((aa,single),(bb,two),(cc,remainders)):
            target.extend(dict(**base,quantity=quantity,**r) for r in rr)
    for e in entries:
        task=e['task'];case=task['case_id'];method=algorithm(task);sr=score_index[task['task_id']]
        base=dict(task_id=task['task_id'],case_id=case,regime=sr['regime'],method=method,canonical_arm=task['arm'],
            truth_access='OFFLINE_SCORING_ONLY',geometry_threshold_selection='Fixed .03/.1/.3 all reported; paper primary .1')
        if not e['sealed']:
            gaps.append(dict(**base,state='final',reason='TERMINAL_PREDICTION_NOT_SEALED',producer_status=(e['receipt'] or {}).get('status')));continue
        d=e['directory'];pred=loadnpz(d/'prediction.npz');m=pred['m'];finalbase=dict(**base,state='final_physical',source=str(d/'prediction.npz'))
        if not material_valid(m):gaps.append(dict(**finalbase,reason='INVALID_FINAL_MATERIAL'));continue
        dm=m-truth[case];decompose(finalbase,'material',dm);decompose(finalbase,'material_real',dm,'real')
        row=dict(**finalbase,material_error_squared_energy=energy(dm),physical_current_source='SAVED_PREDICTION_ONLY')
        j=pred.get('j')
        if current_valid(j):
            dj=j-jtruth[case];decompose(finalbase,'physical_current_error',dj)
            row.update(physical_current_error_squared_energy=energy(dj),receiver_true_error_squared_energy=energy(S@dj),internal_field_true_error_squared_energy=energy(G@dj))
        else:gaps.append(dict(**finalbase,reason='PHYSICAL_CURRENT_NOT_SAVED_NO_RECOMPUTATION'))
        metrics.append(row)
        p=joint_paths.get(task['task_id'])
        if p is None:gaps.append(dict(**base,state='last_joint',reason='NO_SAVED_JOINT_STATE'));continue
        aux=loadnpz(p);meta=readjson(p.with_suffix('.json'));jm=aux['m'];J=aux['J'];Q=aux['Q'];rank=int(meta['receiver_rank']);mf=int(meta['basis_info']['mf'])
        ab=dict(**base,state=p.stem,receiver_rank=rank,mf=mf,source=str(p))
        if not material_valid(jm) or not current_valid(J) or not current_valid(Q):gaps.append(dict(**ab,reason='INVALID_SAVED_JOINT_SCHEMA'));continue
        auxerr=J-jtruth[case];defect=J-jm[:,None]*Q;decompose(ab,'auxiliary_current_error',auxerr);decompose(ab,'physical_consistency_defect',defect)
        ar=dict(**ab,auxiliary_current_error_squared_energy=energy(auxerr),physical_consistency_defect_squared_energy=energy(defect),
            physical_pair_status='MISSING_SAVED_MATCHED_PHYSICAL_STATE',cross_energy=None)
        pp=p.with_name(p.name.replace('_joint','_physical'))
        if pp.exists():
            phys=loadnpz(pp)
            if not np.array_equal(jm,phys.get('m')):
                ar['physical_pair_status']='MATERIAL_MISMATCH_CROSS_NOT_COMPUTED';gaps.append(dict(**ab,reason='JOINT_PHYSICAL_MATERIAL_MISMATCH_NO_FORWARD_REPAIR',physical_source=str(pp)))
            elif current_valid(phys.get('j')):
                phys_error=phys['j']-jtruth[case];gap=J-phys['j'];ar.update(physical_pair_status='EXACT_SAME_M_SAVED_PAIR',
                    physical_current_error_squared_energy=energy(phys_error),auxiliary_minus_physical_current_squared_energy=energy(gap),
                    cross_energy=2*float(np.vdot(phys_error,gap).real))
                decompose(ab,'matched_joint_physical_current_error',phys_error);decompose(ab,'auxiliary_minus_physical_current',gap)
            else:gaps.append(dict(**ab,reason='MATCHED_PHYSICAL_CURRENT_NOT_SAVED_NO_RECOMPUTATION',physical_source=str(pp)))
        else:gaps.append(dict(**ab,reason='MATCHED_PHYSICAL_SNAPSHOT_MISSING_NO_RECOMPUTATION'))
        cache=oracle_cache.get((case,rank,mf))
        if cache:
            floor,JD,q,oracle_projection=cache;ar['shared_oracle_key']=case+'|rank'+str(rank)+'|mf'+str(mf)
            ar['affine_GT_irreducible_squared_energy']=floor['affine_GT_irreducible_squared_energy']
            ar['auxiliary_in_affine_range_error_squared_energy']=energy(J-oracle_projection)
        else:gaps.append(dict(**ab,reason='NO_PREDECLARED_ORACLE_KEY'))
        metrics.append(ar)
    # Pooled fractions use sums of absolute band/total energies, never mean percentages.
    aggregate=defaultdict(lambda:dict(component_energy=0.,total_energy=0.,cases=0))
    for r in single:
        key=tuple(r[k] for k in ('regime','method','state','quantity','tau','operator','band'));v=aggregate[key]
        v['component_energy']+=r['component_energy'];v['total_energy']+=r['total_energy'];v['cases']+=1
    pooled=[]
    for key,v in aggregate.items():pooled.append(dict(zip(('regime','method','state','quantity','tau','operator','band'),key),**v,
        pooled_fraction=v['component_energy']/v['total_energy'] if v['total_energy'] else None))
    for name,rr in [('STATE_ENERGY.csv',metrics),('SINGLE_OPERATOR_ENERGY.csv',single),('TWOFOLD_ENERGY.csv',two),('TWOFOLD_REMAINDER_CROSS.csv',remainders),('ORACLE_CAPACITY_BY_CASE_RANK_MF.csv',oracle),('POOLED_ENERGY.csv',pooled)]:csvwrite(out/name,rr)
    dump(out/'GAPS.json',gaps);dump(out/'GT_CURRENT_OFFLINE_COST.json',dict(rows=gt_cost,LU=GT_LU,LU_dimension=576,complex_RHS=GT_RHS,adjoint_RHS=0,online_inverse_calls=0))
    dump(out/'GEOMETRY_METADATA.json',geom.meta);dump(out/'SOURCE_PROVENANCE.json',sources)
    figure_outputs=paper_figures(args.figure_dir,single,pooled,oracle)
    dump(out/'RECEIPT.json',dict(status='COMPLETE_WITH_EXPLICIT_GAPS' if gaps else'COMPLETE',root_authorization_note=args.authorization_note,
        truth_access='OFFLINE_SCORING_POST_MATERIAL_SCORE_ONLY',expected_tasks=168,expected_cases=24,new_inverse_runs=0,new_online_physical_calls=0,
        offline_scoring_GT_forward_LU=GT_LU,offline_scoring_GT_forward_dimension=576,offline_scoring_GT_complex_RHS=GT_RHS,offline_scoring_adjoint_RHS=0,
        new_endpoint_or_joint_forward_LU=0,geometry_setup=geometry_cost,geometry_instances=1,coefficient_range_cache_entries=len(geom.affines),
        unique_case_rank_mf_oracles=len(oracle),per_arm_oracle_SVD_repetitions=0,primary_tau=.1,reported_taus=[.03,.1,.3],GT_threshold_selection=False,
        final_energy_rows=sum(r['state']=='final_physical' for r in metrics),joint_energy_rows=sum(r['state']!='final_physical' for r in metrics),
        physical_pair_gaps=len(gaps),figures=figure_outputs,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start,
        scope='Auxiliary and matched physical currents separate; final m-derived physical current is outside auxiliary affine constraints; Green material projections are coordinate diagnostics, not material-Jacobian observability. No GT online/root/iterate choice.',
        nested_spectral_geometry_cost=ref.cost))
    print(json.dumps(dict(status='FRESH_ENERGY_COMPLETE',output=str(out),GT_forward_LU=GT_LU,GT_complex_RHS=GT_RHS)))

def paper_figures(figdir,single,pooled,oracle):
    """Only fixed tau.1, absolute-plus-fraction and same-case capacity figures."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    figdir.mkdir(parents=True,exist_ok=True);targets=[figdir/'ENERGY_FIXED_TAU_POINT1.png',figdir/'ORACLE_CAPACITY_BY_CASE.png']
    if any(p.exists() for p in targets):raise FileExistsError('Preserve existing fresh energy paper figures')
    colors=('#c54237','#e9ac3a','#3e79b0');bands=('strong','weak','complement')
    # All24 final endpoints pooled per method, with absolute energies retained.
    def values(method,quantity,operator,band):
        rr=[r for r in single if r['state']=='final_physical' and r['method']==method and r['quantity']==quantity and r['operator']==operator and r['band']==band and r['tau']==.1]
        e=sum(r['component_energy'] for r in rr);t=sum(r['total_energy'] for r in rr);return e,100*e/t if t else 0.,len(rr)
    fig,axs=plt.subplots(2,4,figsize=(17,8),constrained_layout=True);ys=np.arange(len(METHODS))
    for qi,q in enumerate(('material','physical_current_error')):
        for oi,operator in enumerate(('G_in','G_ob')):
            absolute,fraction=axs[qi,2*oi],axs[qi,2*oi+1];left=np.zeros(7)
            for bi,(band,color) in enumerate(zip(bands,colors)):
                vals=[values(m,q,operator,band) for m in METHODS];ee=np.array([x[0] for x in vals]);pct=np.array([x[1] for x in vals])
                absolute.barh(ys+(bi-1)*.23,np.where(ee>0,ee,np.nan),height=.22,color=color);fraction.barh(ys,pct,left=left,height=.7,color=color);left+=pct
            absolute.set_xscale('log');absolute.set_xlabel('sum absolute squared error energy (log)');fraction.set_xlim(0,100);fraction.set_xlabel('pooled error energy fraction (%)')
            for a in (absolute,fraction):a.set_yticks(ys,[LABELS[m]+' n='+str(values(m,q,operator,'strong')[2]) for m in METHODS]);a.invert_yaxis();a.grid(axis='x',alpha=.2)
            absolute.set_title(q+' | '+operator+' absolute',fontsize=9);fraction.set_title(q+' | '+operator+' percent',fontsize=9)
    fig.legend([Patch(color=c) for c in colors],bands,loc='outside lower center',ncol=3)
    fig.suptitle('Fresh24 postseal | fixed diagnostic tau=.1 | material coordinates vs physical current | failures/gaps remain in tables')
    fig.savefig(targets[0],dpi=170);fig.savefig(targets[0].with_suffix('.pdf'));plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(14,5),constrained_layout=True);cases=sorted({r['case_id'] for r in oracle})
    for mf,color in zip((3,6,10,12),('#a174b2','#e9ac3a','#227d83','#333333')):
        rr=[r for r in oracle if r['mf']==mf]
        for rank in sorted({r['receiver_rank'] for r in rr}):
            z=[r for r in rr if r['receiver_rank']==rank];xs=[cases.index(r['case_id']) for r in z]
            axs[0].scatter(xs,[r['actual_coefficient_range_rank'] for r in z],s=22,color=color,label='mf'+str(mf) if rank==min(r['receiver_rank'] for r in rr) else None)
            axs[1].scatter(xs,[max(r['affine_GT_irreducible_squared_energy'],1e-30) for r in z],s=22,color=color)
    for a in axs:a.set_xticks(range(len(cases)),[c.removeprefix('r10_twod_') for c in cases],rotation=90);a.set_xlabel('case; oracle shared across arms at same rank');a.grid(alpha=.2)
    axs[0].set_ylabel('actual affine coefficient-range rank');axs[1].set_ylabel('affine GT irreducible squared current energy');axs[1].set_yscale('log');axs[0].legend()
    fig.suptitle('Actual FFT/APCS capacity and offline oracle | not a bound on final physical current or material recovery')
    fig.savefig(targets[1],dpi=170);fig.savefig(targets[1].with_suffix('.pdf'));plt.close(fig)
    return [str(p) for p in targets]

if __name__=='__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):main()
