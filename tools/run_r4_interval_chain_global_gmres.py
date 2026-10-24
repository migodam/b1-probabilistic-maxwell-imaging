#!/usr/bin/env python3
"""Prospective four-arm same-cache global-Krylov control. No autorun.

Reuses the fast driver's opt-in hooks; its closed/default two-arm behavior is
unchanged. Independent Full/Richardson/global-GMRES/fixed-accuracy processes.
"""
from pathlib import Path
import argparse,json,sys,time,subprocess,os
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'tools')]
import numpy as np
import run_r4_interval_chain_fast as engine
from p3r4.decision_refresh_global_gmres import DecisionRefreshGlobalGMRES,resolve_global
from p3r4.evaluation import checkpoint_identity

EXPERIMENT='interval_chain_global_gmres_v1';ARMS=('full','richardson','global_gmres','global_gmres_accuracy')
engine.EXPERIMENT=EXPERIMENT;engine.PARENTS=tuple(f'R4DEV_{i:05d}' for i in range(48,56));engine.SEED=20261025161
engine.RESULT=engine.OUT/'results/II'/EXPERIMENT
engine.CAP={'full':65536,'richardson':78154,'global_gmres':78155,'global_gmres_accuracy':78155}
engine.CACHE_CLASS_BY_ARM={a:DecisionRefreshGlobalGMRES for a in ARMS if a.startswith('global_gmres')}


def configure(cache,book,costs):
    def small_solve(H,rhs):
        before=costs['RHS'];book.reserve_rhs(1,'samecache_global_GMRES_reduced_LS')
        with book.span('global_gmres_reduced_complex_least_squares',reduced_LS_RHS=1,reduced_LS_rows=H.shape[0],reduced_LS_columns=H.shape[1]):
            tick=time.perf_counter();result=np.linalg.lstsq(H,rhs,rcond=None)[0];numeric=time.perf_counter()-tick
        cache._last_reduced_numeric_seconds=numeric
        cache.counts['global_gmres_reduced_LS_numeric_nanoseconds']+=int(numeric*1e9)
        return result,costs['RHS']-before
    cache.configure_reduced_solver(small_solve)


engine.CACHE_CONFIGURE_BY_ARM={a:configure for a in engine.CACHE_CLASS_BY_ARM}
engine.RESOLVER_BY_ARM={'global_gmres':resolve_global,'global_gmres_accuracy':lambda *a,**kw:resolve_global(*a,**kw,fixed_accuracy=True)}
engine.FROZEN_DESIGN_EXTENSION=dict(comparison='same-paid576columnR0/full8sourceBLASrights',arms=ARMS,
  source_pCN_beta=.002,max_global_GMRES_iterations=4,global_vector='one C^(576*8) vector, Frobenius complex innerproduct',
  small_LS='one actual complex reduced RHS periteration conservatively charged book.reserve_rhs1; smallshape/wall separate, never FullLUchargeforR0matvec',
  basis_cache='bounded current+candidate endpoint bases only, max2, originalR0 unchanged',
  fixed_accuracy_arm='BOTH current/candidate error_white_upper<=.05 AND Phiwidth<=.05 beforedecision; then sameU MH interval refinement/fallback',
  success='allfourarms flags/paths/ALLproofintervals match independent Full; wall/RHS/fallbacks descriptive, parentownsbenefitclaim',
  no_A17_or_current_space_advantage_claim=True,no_new_physics_until_parentrelease=True)


class GlobalProofWriter(engine.EndpointWriter):
    def proof(self,state,role,event,step):
        super().proof(state,role,event,step)
        profile=state.point.get('global_gmres_profile')
        if profile is not None:engine.append(self.folder/'GLOBAL_GMRES_REDUCED_PROOFS.jsonl',dict(step=step,endpoint_id=state.endpoint_id,role=role,event=event,**profile))


engine.EndpointWriter=GlobalProofWriter


def compare():
    root=engine.RESULT
    if not all((root/a/'ARM_COMPLETE.json').exists() for a in ARMS):raise ValueError('allfourarms closed before otherarm reads')
    gates={a:json.loads((root/a/'ARM_COMPLETE.json').read_text()) for a in ARMS}
    if any(g['status']!='COMPLETE' or g['array_codec']!='stored' for g in gates.values()) or len({g['pid'] for g in gates.values()})!=4:raise ValueError('four independent closedCOMPLETE storedarms')
    for arm,gate in gates.items():
        if gate['input_identities']!=gates['full']['input_identities'] or json.loads(Path(gate['receipt']).read_text())['status']!='COMPLETE':raise ValueError('commoninputs andclosedBook required')
        for identity in gate['all_array_seals']:
            if checkpoint_identity(identity['path'])!=identity:raise ValueError('sealedarray changed')
    timing=json.loads((root/'ARM_PROCESS_WALL.json').read_text());rows=[];total_misses=0
    for parent in engine.PARENTS:
        with np.load(root/'full'/parent/'CHAIN.npz') as z:f={k:z[k] for k in z.files}
        gold={'initial':float(f['current_phi_before'][0])};gold.update({f'c{i:04d}':float(v) for i,v in enumerate(f['candidate_phi'])})
        for arm in ARMS[1:]:
            with np.load(root/arm/parent/'CHAIN.npz') as z:q={k:z[k] for k in z.files}
            flags=np.flatnonzero(q['accept']!=f['accept']);paths=np.flatnonzero(np.any(q['x_path']!=f['x_path'],axis=1));proof=[]
            for line in (root/arm/parent/'ENDPOINT_PROOFS.jsonl').read_text().splitlines():
                item=json.loads(line);lo,hi=item['Phi_interval'];lo=0. if lo is None else lo;hi=float('inf') if hi is None else hi
                if not lo<=gold[item['endpoint_id']]<=hi:proof.append(dict(step=item['step'],endpoint_id=item['endpoint_id'],refinement=item['accuracy_iteration']))
            candidate=np.flatnonzero((f['candidate_phi']<q['candidate_interval'][:,0])|(f['candidate_phi']>q['candidate_interval'][:,1]));current=np.flatnonzero((f['current_phi_before']<q['current_interval'][:,0])|(f['current_phi_before']>q['current_interval'][:,1]))
            accuracy=[]
            if arm=='global_gmres_accuracy':
                for line in (root/arm/parent/'CHAIN_REFERENCES.jsonl').read_text().splitlines():
                    step=json.loads(line);final=step['history'][-1]
                    if not final.get('current_accuracy') or not final.get('candidate_accuracy'):accuracy.append(step['step'])
            count=len(flags)+len(paths)+len(proof)+len(candidate)+len(current)+len(accuracy);total_misses+=count
            rows.append(dict(parent=parent,arm=arm,flag_mismatches=flags.tolist(),path_mismatches=paths.tolist(),all_proof_misses=proof,
                             candidate_interval_misses=candidate.tolist(),current_interval_misses=current.tolist(),fixed_accuracy_violations=accuracy))
    combined=sum(g['costs']['RHS'] for g in gates.values())
    if combined>300000:raise AssertionError('combined300kRHS cap')
    table={arm:dict(actual_RHS=gates[arm]['costs']['RHS'],whole_process_seconds=timing[arm]['seconds'],relative_wall_saving_vsFull=1-timing[arm]['seconds']/timing['full']['seconds']) for arm in ARMS}
    result=dict(status='PASS_SAME_KERNEL_FOUR_ARM_NUMERICS' if total_misses==0 else 'NO_GO_NUMERICAL_OR_PATH_MISMATCH',total_misses=total_misses,rows=rows,costs=table,combined_RHS=combined,
                scope='ordinarysamecache numerical control; no forcedRichardson advantage, no A17/current-space claim, no posteriorreference; float64notcertified')
    from p3r4.core import write
    if (root/'COMPARISON.json').exists():raise FileExistsError('preserve closedcomparison')
    write(root/'COMPARISON.json',result);print(json.dumps(result),flush=True)
    if total_misses:raise SystemExit(2)


def pair(release):
    check=json.loads(release.read_text())
    if check.get('released') is not True or check.get('experiment')!=EXPERIMENT or check.get('campaign_id')!='P3LR4-20261025':raise ValueError('matchingparentrelease required')
    marker=engine.OUT/'results/II/interval_chain_global_gmres_pair'
    if marker.exists() or engine.RESULT.exists():raise FileExistsError('no paidreplay')
    marker.mkdir(parents=True);records=[];timings={};environment=dict(os.environ,R4_ARRAY_CODEC='stored')
    for arm in (*ARMS,'compare'):
        start=time.perf_counter();p=subprocess.run([sys.executable,'-u',str(Path(__file__).resolve()),arm,'--release',str(release.resolve())],cwd=ROOT,env=environment)
        timings[arm]=dict(seconds=time.perf_counter()-start,returncode=p.returncode);records.append(dict(arm=arm,returncode=p.returncode))
        if arm==ARMS[-1] and p.returncode==0:(engine.RESULT/'ARM_PROCESS_WALL.json').write_text(json.dumps(timings,indent=2)+'\n')
        if p.returncode:(marker/'FAILED.json').write_text(json.dumps(dict(status='FAILED',records=records),indent=2)+'\n');raise SystemExit(p.returncode)
    (marker/'COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',records=records,comparison=str(engine.RESULT/'COMPARISON.json')),indent=2)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('arm',nargs='?',default='pair',choices=('pair',*ARMS,'compare'));p.add_argument('--release',type=Path);a=p.parse_args()
    if a.arm=='compare':compare()
    elif a.release is None:p.error('parentrelease required; noautorun')
    elif a.arm=='pair':pair(a.release)
    else:engine.run(a.arm,a.release)

if __name__=='__main__':main()
