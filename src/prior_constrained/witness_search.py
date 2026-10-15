"""Bounded offline P0 pair construction, exact A30 full solver, no truth guide.

This is benchmark construction, not an online generator result. Generator
endpoints remain in the original distribution's parameter support. Full
Jacobian ranking is a proposal tool, never an acceptance certificate.
"""
from pathlib import Path
import copy,csv,json,sys,time
import numpy as np
from scipy.stats import chi2,binomtest
from . import prior
from .common import ROOT,OUT,Book,dump,plain,descriptor,separation,ell

TAU=float(np.sqrt(chi2.ppf(.95,512)))
LIMITS={'center':(-.12,.12),'rotation':(-np.pi,np.pi),'scale':(.78,1.04),'axis_ratio':(.82,1.18),
 'shell_radius':(.14,.19),'shell_width':(.035,.055),'disk_radius':(.065,.10),'disk_weight':(.45,.8),
 'second_shell_weight':(.45,.75),'rectangle_width':(.17,.25),'rectangle_height':(.19,.29),
 'gaussian_centers':(-.13,.13),'gaussian_widths':(.05,.095),'gaussian_weights':(.4,.9),'segment_offset':(-.06,.06)}

def mutate(row,rng,scale,index):
    r=copy.deepcopy(row); r['parent_id']=f"profile_{index}"
    r['parent_seed']=None; r['proposal_provenance']={'source_parent':row['parent_id'],'scale_fraction':scale,'index':index}
    for key,(lo,hi) in LIMITS.items():
        x=np.asarray(row['geometry'][key]); low=np.maximum(lo+1e-9,x-scale*(hi-lo)); high=np.minimum(hi-1e-9,x+scale*(hi-lo))
        r['geometry'][key]=np.asarray(rng.uniform(low,high)).tolist()
    if rng.uniform()<.15: r['geometry']['subtype']=1-int(r['geometry']['subtype'])
    return r

def save(path,state,row,**meta):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path,chi=state.chi,prediction=state.data,field=state.field,current=state.current,backward_residual=np.asarray(state.residual))
    dump(path.with_suffix('.json'),dict(generator=row,**meta))

def replay_noise(fa,fb,sigma,pair_id):
    records=[]
    for truth in ['a','b']:
        clean=fa if truth=='a' else fb
        for i in range(32):
            seed=[20261015,3201,i]
            rng=np.random.default_rng(np.random.SeedSequence(seed))
            y=clean+sigma/np.sqrt(2)*(rng.standard_normal(clean.shape)+1j*rng.standard_normal(clean.shape))
            ea,eb=ell(fa,y,sigma),ell(fb,y,sigma); ra,rb=np.sqrt(2*ea),np.sqrt(2*eb)
            records.append(dict(pair_id=pair_id,truth_endpoint=truth,noise_index=i,seed=json.dumps(seed),sigma_fixed=sigma,
                residual_a=float(ra),residual_b=float(rb),ell_a=ea,ell_b=eb,delta_ell=eb-ea,
                both_noise_ball=bool(max(ra,rb)<=TAU),near_likelihood=bool(abs(eb-ea)<=np.log(10)),
                competitive_primary=bool(max(ra,rb)<=TAU and abs(eb-ea)<=np.log(10)),
                competitive_log100=bool(max(ra,rb)<=TAU and abs(eb-ea)<=np.log(100))))
    return records

def freeze():
    p=OUT/'PROTOCOL_FROZEN.json'
    if p.exists(): raise RuntimeError('immutable protocol already exists')
    cal=json.loads((OUT/'delegated/prior/STRUCTURE_CALIBRATION.json').read_text())
    protocol=dict(identity='A25-A30-P0-EXACT-PCSP-20261015-v1',prior='original four-family TRAIN generator; original parameters and per-image maximum normalization; no hole-count labels',
        physics=dict(grid=24,sources=8,receivers=32,real_channels=512,gain=1,precision='complex128',packing='source-major Re32 then Im32',whitening='sqrt(2)/original complex-RMS sigma',noise='original differential RMS at20dB, fixed from anchor, never candidate rescaled'),
        structure=dict(calibration=cal,joint_rule='all three calibrated minima AND persistent disagreement >=.05 on >=3/5 thresholds AND >=6 cells of normalized rho change >=.1; DC excluded from centered shape'),
        E1=dict(purpose='offline benchmark construction, not free online discoveries',anchors='first fresh reservoir parent of each family: shell/ring, asymmetric/weak, smooth/weak, mixed/ring',
            proposal_seed=2026101531,local_material_proposals_per_anchor_per_round=128,rounds=2,
            round1='rank 128 exact local generator secants using paid exact material J; full-evaluate 3 local +3 independent global reservoir endpoints',
            round2='refresh exact J at closest structurally qualified first-round local endpoint; rank128 local parameter changes relative fixed original anchor; full-evaluate6',
            parameter_scale_fractions=[.08,.16,.32,.48],full_state_cap=64,RHS_cap=2048,
            main_search_planned_states=52,remaining12='only nonzero .02 parameter-range prior-local perturbation of an actually qualifying closest pair',
            preferred_d_pred_max=2,noise_replays_per_truth=32,noise_seed_domain=[20261015,3201],
            repeat_competition='both noise balls and |Delta ell|<=log(10) on at least16/32 fixed noises separately for each endpoint truth; log100 sensitivity only',
            E2_requires='two parent-disjoint exact-P0 pairs satisfying structural, d_pred<=2 and repeated competition; no obvious extreme-parameter or numerical defect'),
        E2=dict(status='CONDITIONAL_NOT_ACTIVE',full_state_cap=192,RHS_cap=8192,arms=['prior-respecting multistart','T-soft+C_P','QD adapted+C_P','PC-SP+C_P'],max_observations=4,max_states_per_arm_observation=10,max_proposals_per_arm=3,max_accepted_corrections=2,eta=.5),
        all_stage_caps=dict(states=256,RHS=10240,GPU_training=0),material_only_cap=2048,
        material_accounting='256 fresh generator parents +256 diagnostic transforms +up to1024 generator secants +up to12 robustness perturbations <=1548; repeated rendering of an identical rho in another known background is not a new parent',
        source='CODEX_NEXT_ROUND_ZH.md',posterior_mass=None,locked_access=False,new_noise_relaxation=False)
    dump(p,protocol)

def run():
    if (OUT/'E1_DECISION.json').exists(): raise RuntimeError('completed E1 cannot silently repeat')
    if not (OUT/'PROTOCOL_FROZEN.json').exists(): freeze()
    book=Book('E1'); book.check(); start=time.perf_counter()
    src=prior.A30/'src'; sys.path.insert(0,str(src))
    from a30.physics import FrozenMaxwellAdapter
    reservoir=json.loads((OUT/'delegated/prior/MATERIAL_ONLY_PROVENANCE.json').read_text())
    if isinstance(reservoir,dict):
        for k in ['rows','parents','samples','materials']:
            if k in reservoir: reservoir=reservoir[k]; break
    if not isinstance(reservoir,list): raise ValueError('reservoir parent list required')
    rows=[]; replays=[]; pairs=[]; models={}; emitted=0; cache={}
    configs=[(0,'ring'),(1,'weak'),(2,'weak'),(3,'ring')]
    # Resume uses the SAME authorization/ledger and preserves failed arrays.
    execution=1+sum(r['event']=='JOB_BEGIN' for r in book.rows())
    dest=OUT/'offline_only'/f'E1_execution_{execution}'; dest.mkdir(parents=True,exist_ok=True)
    if any(dest.glob('*.npz')): raise RuntimeError('completion arrays exist; explicit array-safe resume required')
    book.append(dict(event='JOB_BEGIN',name='E1_bounded_original_P0_search_completion',failure_history_preserved=True))
    def full(model,row,bg,tag):
        nonlocal emitted
        material=prior.decode(row,background_id=bg)
        key=(bg,material.tobytes())
        if key in cache:
            book.append(dict(event='CACHE_REUSE',tag=tag,source=cache[key][2],new_states=0,new_RHS=0)); return cache[key][0],cache[key][1]
        book.context=dict(tag=tag,background=bg,role='OFFLINE_BENCHMARK_CONSTRUCTION')
        s=model.full_state(material,purpose=tag); p=dest/f'{tag}.npz'; save(p,s,row,background_id=bg)
        replay=prior.verify_exact_endpoint(s.chi,row,background_id=bg)
        if replay['prior_status']!='P0_EXACT': raise RuntimeError('generated endpoint lost P0 support')
        if s.residual>1e-10: raise RuntimeError('full solver backward error')
        cache[key]=(s,p,tag)
        book.append(dict(event='CANDIDATE',tag=tag,array=str(p.relative_to(ROOT)),prior_status='P0_EXACT',backward_residual=s.residual,posterior_mass=None))
        return s,p
    def linear(model,state):
        tr=model.receiver_transfer(state)
        return np.einsum('qp,ps->qsp',tr,state.field).reshape(256,576)
    rng=np.random.default_rng(2026101531)
    for fi,bg in configs:
        anchor_row=next(r for r in reservoir if r['family']==prior.FAMILIES[fi])
        if bg not in models: models[bg]=FrozenMaxwellAdapter.from_background_npz(prior.INPUTS/'background'/f'{bg}.npz',book,purpose_prefix='PRIOR_CONSTRAINED_E1')
        model=models[bg]; anchor,path=full(model,anchor_row,bg,f'f{fi}_anchor')
        sigma=float(np.sqrt(np.mean(abs(anchor.data-model.yb)**2))*.1)
        b=prior.background(bg); da=descriptor(anchor.chi,b); jac=linear(model,anchor)
        qualified=[]; current=anchor; current_row=anchor_row
        for roundno in range(2):
            proposals=[]
            for j in range(128):
                emitted+=1; scale=[.08,.16,.32,.48][j%4]
                rr=mutate(current_row,rng,scale,emitted)
                x=prior.decode(rr,background_id=bg); dd=descriptor(x,b); sep=separation(da,dd)
                if sep['structure_qualified']:
                    pred=current.data.reshape(-1)+jac@(x-current.chi)
                    score=float(np.linalg.norm(pred-anchor.data.reshape(-1))*np.sqrt(2)/sigma)
                    proposals.append((score,rr,x,sep,'local'))
            proposals.sort(key=lambda r:r[0])
            if roundno==0:
                local=proposals[:3]; glob=[]
                for rr in reservoir:
                    x=prior.decode(rr,background_id=bg); sep=separation(da,descriptor(x,b))
                    if sep['structure_qualified']:
                        score=float(np.linalg.norm(jac@(x-anchor.chi))*np.sqrt(2)/sigma)
                        glob.append((score,rr,x,sep,'global'))
                glob.sort(key=lambda r:r[0]); chosen=local+glob[:3]
            else: chosen=proposals[:6]
            seen=set()
            for j,(score,rr,x,sep,role) in enumerate(chosen):
                if x.tobytes() in seen: continue
                seen.add(x.tobytes()); tag=f'f{fi}_r{roundno}_{role}_{j}'
                ss,pp=full(model,rr,bg,tag)
                d=float(np.linalg.norm(ss.data-anchor.data)*np.sqrt(2)/sigma)
                row=dict(pair_id=tag,anchor_id=f'f{fi}_anchor',anchor_parent=anchor_row['parent_id'],candidate_parent=rr['parent_id'],family=prior.FAMILIES[fi],background=bg,
                    role=role,round=roundno,predicted_d=score,d_pred=d,sigma_fixed=sigma,prior_status='P0_EXACT',
                    anchor_array=str(path.relative_to(ROOT)),candidate_array=str(pp.relative_to(ROOT)),**sep)
                rows.append(row); qualified.append((d,ss,rr,row))
                with (OUT/'CANDIDATE_LEDGER.jsonl').open('a') as f: f.write(json.dumps(plain(dict(stage='E1',evaluation_role='OFFLINE_WITNESS_CONSTRUCTION',posterior_mass=None,**row)))+'\n')
            if qualified:
                _,current,current_row,_=min(qualified,key=lambda t:t[0])
                if roundno==0: jac=linear(model,current)
        if qualified:
            d,ss,rr,record=min(qualified,key=lambda t:t[0]); replay=replay_noise(anchor.data,ss.data,sigma,record['pair_id']); replays.extend(replay)
            fractions={truth:sum(r['competitive_primary'] for r in replay if r['truth_endpoint']==truth)/32 for truth in ['a','b']}
            p=dict(record,noise_competition_fractions=fractions,qualified=bool(d<=2 and record['structure_qualified'] and min(fractions.values())>=.5),posterior_mass=None,
                status='FULL_VALIDATED_NUMERICAL',numerical_error_bound=None)
            pairs.append(p)
    # All candidate endpoints and pairs are frozen before opening any evaluator.
    with (OUT/'E1_SEARCH.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    with (OUT/'PAIR_NOISE_REPLAY.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(replays[0])); w.writeheader(); w.writerows(replays)
    good=[p for p in pairs if p['qualified']]
    # Local robustness is paid only if a qualifying pair actually exists.
    robust=[]
    for pair in good[:2]:
        if book.totals('E1')['states']>=64: break
        candidate_path=ROOT/pair['candidate_array']; meta=json.loads(candidate_path.with_suffix('.json').read_text()); rr=meta['generator']; bg=pair['background']; model=models[bg]
        with np.load(ROOT/pair['anchor_array']) as z: fa=z['prediction'].copy()
        for k in range(6):
            if book.totals('E1')['states']>=64: break
            emitted+=1; pp=mutate(rr,rng,.02,emitted); ss,where=full(model,pp,bg,f"robust_{pair['pair_id']}_{k}")
            robust.append(dict(pair_id=pair['pair_id'],scale_fraction=.02,d_pred=float(np.linalg.norm(ss.data-fa)*np.sqrt(2)/pair['sigma_fixed']),array=str(where.relative_to(ROOT))))
    dump(OUT/'WITNESS_PAIRS.json',dict(visibility='OFFLINE_ONLY; no E2 online reader authorized',posterior_mass=None,pairs=pairs,local_robustness=robust))
    # The present run has zero qualifying pairs. A future positive branch must
    # check full parent ancestry, local robustness and prior-tail diagnostics;
    # raw pair count alone must never open E2. No such positive result is
    # implemented/claimed by this closed negative pilot.
    decision=dict(status='PROVISIONAL_PAIRS_REQUIRE_REMAINING_PROTOCOL_CHECKS' if good else 'NO_GO_FOR_GENERATOR_COMPARISON',
        material_state='ONE_SUPPORTED_CANDIDATE_AMBIGUITY_UNRESOLVED',qualified_parent_disjoint_pairs=0,
        provisional_pair_count=len(good),E2_authorized_by_evidence=False,
        finite_search_not_global_uniqueness=True,new_generator_secants=emitted,material_only_conservative_total=512+emitted,
        costs=book.totals('E1'),controlled_wall_seconds=time.perf_counter()-start,
        minimum_d_pred=min(r['d_pred'] for r in rows),all_pairs=pairs)
    dump(OUT/'E1_DECISION.json',decision)
    book.append(dict(event='JOB_END',name='E1_bounded_original_P0_search',wall_seconds=decision['controlled_wall_seconds'],counts=dict(material_only_generator_secants=emitted)))
    print(json.dumps(plain({k:v for k,v in decision.items() if k!='all_pairs'})))
if __name__=='__main__': run()
