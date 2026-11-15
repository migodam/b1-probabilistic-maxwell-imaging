"""Freeze six fixed rank19 controls on the existing exposed24 public objects."""
import argparse,json
from copy import deepcopy
from datetime import datetime,timezone
from pathlib import Path

ARMS=('g19_internal_galerkin','g19_observation_galerkin','g19_internal_right','g19_observation_right','g19_anchor_observation_galerkin','g19_anchor_observation_right')

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,required=True);a=ap.parse_args();s=a.root.resolve();R=s.parent
    source=R/'configs/independent24_v1.json';original=json.loads(source.read_text());templates={t['case_id']:t for t in original['tasks'] if t['method']=='g19'}
    cases=[f'r14_new_{i:04d}' for i in range(24)]
    if set(templates)!=set(cases):raise ValueError('Complete original24 G19 public task headers required')
    campaign='observable19_exposed24_v1';tasks=[]
    for cid in cases:
        for arm in ARMS:
            t=deepcopy(templates[cid]);t.update(task_id=cid+'__'+arm,campaign=campaign,method=arm,label=arm,algorithm='g19',basis_rank=19,
                cohort='EXPOSED24_FOLLOWUP',initialization_scope='Original G19 paid physical-current warm start; unchanged uniform public m0=.1',
                runner_family='ORIGINAL_C2_WITH_PRIVATE_FIXED19_FACTORY',basis_family='anchor_observation' if '_anchor_observation_' in arm else 'observation' if '_observation_' in arm else 'internal',
                denominator_projection='right' if arm.endswith('_right') else 'galerkin',basis_scope='predeclared_resolvent_visible_control' if '_anchor_observation_' in arm else 'primary_four_basis_projection_control')
            path=Path(t['input_file'])
            if path.parts[:2]!=('data_public','independent24'):raise ValueError('Unexpected original public input path')
            t['input_file']='../'+path.as_posix();t['observation_contract']['input_file']=t['input_file'];tasks.append(t)
    value=dict(campaign=campaign,status='FROZEN_BEFORE_FOLLOWUP_RUNS',exposure_status='EXPOSED24_FOLLOWUP',original_cohort_was_independent=True,new_independent_or_blind_claim=False,
        source_manifest='../configs/independent24_v1.json',source_campaign=original['campaign'],protocol=deepcopy(original['protocol']),protocol_deepcopy_exact=True,
        methods=list(ARMS),primary_four=list(ARMS[:4]),predeclared_resolvent_visible_controls=list(ARMS[4:]),basis_rank=19,material_real_DOF=576,
        basis_contract={'internal':'first19 G right singular vectors; direct original prepare_full when its .1 strong rank equals19; otherwise fail, never lower rank',
                        'observation':'first19 unscaled S right singular vectors; same right subspace as common scalar W S; no inverse-singular weighting; source metric untouched',
                        'anchor_observation':'first19 unscaled O0 right singular vectors (same common-scalar whitened subspace); O0=S(I-D_m0 G)^-1, m0=original public uniform .1; cached initialLU only, one nr-adjoint-RHS solve, zero additionalLU, no full inverse, NOT a material Jacobian',
                        'galerkin':'K=P G P; B=U(U* G U)', 'right':'K=G P; B=G U',
                        'rank_policy':'fixed19; insufficient requested basis produces a retained failure, never rank adjustment'},
        fixed_readout_snapshots=dict(calls=[1,128,512,1500],fields=['J','Q','m_before','m_after_if_returned','beta','rho','Os','Ls','U','B'],kind='MATERIAL_READOUT_INPUT_OUTPUT_NOT_FULL_PR_CHECKPOINT',all_six_arms_same=True,CPU_and_IO_paid=True,interrupted='Input-only INTERRUPTED; never fabricate m_after'),
        selector='Same original paid minPhi/secondary/prefix rules; no GT use or parameter changes',budget='Original C2 same-process absolute120 and proposal<=72 with remaining original corrector; imports/setup/snapshots and failures paid',
        physics='Original full scalar TM Maxwell model; no substitution of S or O0 for physical G; 576 real material DOF',
        fixed_expansion_policy='Exactly three bases x two projection modes, no further expansion/sweep after freeze',
        basis_spectral_audit='Save complete already-computed G/S/O0 spectra, sigma19, sigma20, (sigma19-sigma20)/sigma1. No extra SVD and no rank adjustment from gap.',
        portable_runtime_source_map='audit/RUNTIME_SOURCE_MAP.json; exact copied R13 full_current/field_current/material_events; remote legacy need not contain R13 package',tasks=tasks)
    target=s/'configs'/f'{campaign}.json'
    with target.open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')
    receipt=dict(created_utc=datetime.now(timezone.utc).isoformat(),source=str(source),source_bytes=source.stat().st_size,manifest=str(target),cases=24,tasks=len(tasks),methods=list(ARMS),protocol_exactly_equal=value['protocol']==original['protocol'],public_input_references_only=True,truth_access='NONE',data_generated=False,launches=0,scientific_scope='EXPOSED24_FOLLOWUP, not fresh independent confirmation')
    with (s/'audit/MANIFEST_FREEZE.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
    print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
