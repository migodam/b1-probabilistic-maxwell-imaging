"""Prepare frozen public-only dual-metric manifests; never run or score tasks."""
from __future__ import annotations
import argparse
from copy import deepcopy
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil

PILOT = (2, 5, 7, 16, 17, 24, 29, 37, 38, 42, 44, 45)
PAIRS = ('00', '01', '10', '11')
READOUT = dict(aplus_mode='mobius', aplus_arc_candidates=True,
               aplus_inherited_gtol=True, aplus_sweeps=2, aplus_step_tol=1e-8)


def write_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def prepare(root, history):
    root, history = Path(root).resolve(), Path(history).resolve()
    if root == history or root.is_relative_to(history):
        raise ValueError('R14 output must be separate from the historical R13 tree')
    source = history / 'data_public/fresh48'
    historical = json.loads((history / 'configs/confirmation_core48.json').read_text())
    contracts = {}
    for task in historical['tasks']:
        if task.get('method') == 'g19':
            contracts[task['case_id']] = deepcopy(task.get('observation_contract', {}))
    ids = [f'r13_fresh_{i:04d}' for i in range(48)]
    sources = [source / (cid + '.npz') for cid in ids]
    if any(not path.is_file() for path in sources):
        raise FileNotFoundError('All 48 original public files are required')
    if set(contracts) != set(ids):
        raise ValueError('Original 48 public observation contracts must all be present')
    target = root / 'data_public/exposed48'
    campaigns = [('dualmetric48_v1', ids),
                 ('dualmetric_pilot12_v1', [ids[i] for i in PILOT]),
                 ('dualmetric_remaining36_v1', [cid for i, cid in enumerate(ids) if i not in PILOT])]
    outputs = [root / 'configs' / (name + '.json') for name, _ in campaigns]
    receipt_path = root / 'audit/CAMPAIGN_PREPARATION.json'
    if target.exists() or any(p.exists() for p in outputs + [receipt_path]):
        raise FileExistsError('Prepared public data/manifests are never overwritten')
    protocol = deepcopy(historical['protocol'])
    protocol.update(CPU_cap=120, proposal_fraction=.6, steps=3000,
                    audit_every=50, physical_audit_CPU_interval=2., BLAS_threads=1,
                    P_cores=list(range(0, 16, 2)), sourcewise_cg=False,
                    current_solver='pr', rho=1., b0=1.,
                    stages=[dict(mf=mf, maxiter=its, beta_current=1.,
                                 beta_material=b, b0=1., b1=b)
                            for mf, its, b in [(3, 500, 1.), (6, 500, .2),
                                               (10, 500, .2), (12, 1500, .2)]],
                    coefficient_tol=0., miniter=50,
                    material_inner_maxiter=60, material_inner_ftol=1e-12,
                    material_inner_gtol=1e-8, initialize_physical_current=False,
                    readout_options=READOUT,
                    material_M0_beta='b0=1 in every stage; not historical Dual beta_material',
                    material_M1_beta='b1=1,.2,.2,.2 from the original A+ schedule',
                    scope='EXISTING_EXPOSED48_MECHANISM_STUDY_NOT_FRESH_BLIND',
                    truth_policy='NONE online; seal new predictions before independent offline scoring',
                    source_scope='R14 dual-metric implementation over inherited public R13 inputs',
                    lane_policy='Case group pinned to one core; four sequential pairs then next case',
                    priority='NORMAL')
    # Preserve provenance but do not mislabel the new code with R13 source_commit.
    if 'source_commit' in protocol:
        protocol['historical_source_commit'] = protocol.pop('source_commit')
    target.mkdir(parents=True, exist_ok=False)
    copies = []
    for src in sources:
        dst = target / src.name
        shutil.copy2(src, dst)
        copies.append(dict(source=str(src), destination=str(dst), bytes=dst.stat().st_size))
    for (name, members), dest in zip(campaigns, outputs):
        tasks = []
        for cid in members:
            for pair in PAIRS:
                tasks.append(dict(task_id=cid + '__pair' + pair, case_id=cid,
                    method='pair' + pair, algorithm='pair' + pair, pair=pair,
                    current_metric='M' + pair[0], material_metric='M' + pair[1],
                    campaign=name, cohort='EXISTING_EXPOSED48', physics_kind='scalar_tm_2d',
                    input_file='data_public/exposed48/' + cid + '.npz', seed=13013,
                    CPU_cap=120, sourcewise_cg=False, initialize_physical_current=False,
                    readout_version='R7_FINAL_READOUT' if pair[1] == '1' else 'ORDINARY_SCALAR_B0',
                    readout_options=deepcopy(READOUT) if pair[1] == '1' else {},
                    observation_contract=deepcopy(contracts[cid])))
        write_new(dest, dict(campaign=name, status='FROZEN_BEFORE_ANY_R14_RUN',
            scope='EXISTING_EXPOSED48_MECHANISM_STUDY_NOT_FRESH_BLIND',
            protocol=deepcopy(protocol), frozen_readout_options=deepcopy(READOUT),
            fixed_pilot_case_indices=list(PILOT), tasks=tasks))
    receipt = dict(status='PREPARED_NOT_EXECUTED', created_utc=datetime.now(timezone.utc).isoformat(),
        history_root=str(history), root=str(root), public_copies=copies,
        manifests=[str(p) for p in outputs], task_counts=[192, 48, 144],
        truth_access='NONE', scoring='NOT_RUN', optimization='NOT_RUN',
        copy_policy='Original public NPZ bytes copied; no array transformation or truth read',
        partition_policy='Pilot12 plus remaining36 equals main48; never run main in addition to partition unless separately requested')
    write_new(receipt_path, receipt)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--history-root', type=Path, required=True,
                        help='Historical R13 study root, not the repository root')
    args = parser.parse_args()
    receipt = prepare(args.root, args.history_root)
    print(json.dumps(dict(status=receipt['status'], manifests=receipt['manifests'],
                          truth_access='NONE')))


if __name__ == '__main__':
    main()
