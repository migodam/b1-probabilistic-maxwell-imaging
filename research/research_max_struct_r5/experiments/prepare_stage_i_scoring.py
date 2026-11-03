"""Prepare fresh isolated Stage I scoring inputs; stdlib only, no GT reads."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import json

ROOT=Path(__file__).resolve().parents[1]

def prepare(root):
    root=Path(root).resolve()
    audit=root/'audit'
    seal_path=audit/'STAGE_I_SEAL.json'
    receipt_path=audit/'OFFLINE_SCORING_RECEIPT.json'
    if seal_path.exists() or receipt_path.exists():
        raise FileExistsError('Use a fresh isolated copy: existing seal/receipt must never be overwritten')
    paths=set()
    def add(path):
        if not path.is_file():raise FileNotFoundError(path)
        paths.add(path)
    def completed(folder,status):
        rec=folder/'RECEIPT.json'
        if json.loads(rec.read_text())['status']!=status:raise ValueError('Incomplete output: '+str(folder))
        for path in folder.rglob('*'):
            if path.is_file() and '__pycache__' not in path.parts:add(path)
    replay_config=root/'configs/trajectory_replay_frozen.json'
    factorial_config=root/'configs/factorial_frozen.json'
    replay=json.loads(replay_config.read_text())
    factorial=json.loads(factorial_config.read_text())
    for run in replay['runs']:
        completed(root/'RESULTS/trajectory_replay'/run['run_id'],'COMPLETED_AUX_ONLY')
        completed(root/'RESULTS/trajectory_diagnostics'/run['run_id'],'COMPLETED')
    for arm in ('self_absorption','onehop_phase_pi'):
        completed(root/'RESULTS/trajectory_diagnostics'/('legacy_'+arm),'COMPLETED_FIXED_STATES')
    for arm in factorial['arms']:
        completed(root/'RESULTS/factorial'/(factorial['case_id']+'__'+arm),'COMPLETED_DEVELOPMENT_FACTORIAL')
    # Freeze the exact present source and configuration, not a stale original snapshot.
    for path in (replay_config,factorial_config,root/'STAGE_I_PROTOCOL.md',root/'THEORY_DERIVATIONS.md'):add(path)
    for folder in (root/'src',root/'experiments'):
        for path in folder.glob('*.py'):add(path)
    for name in ('TRAJECTORY_REPLAY_SUMMARY.json','FACTORIAL_SUMMARY.json'):
        path=root/'RESULTS'/name
        if path.exists():add(path)
    files=[]
    for path in sorted(paths):
        stat=path.stat()
        files.append(dict(path=str(path.relative_to(root)),bytes=stat.st_size,mtime_ns=stat.st_mtime_ns))
    now=datetime.now(timezone.utc).isoformat()
    seal=dict(status='STAGE_I_DIAGNOSTIC_AND_FIRST_FACTORIAL_SEALED',sealed_utc=now,
        truth_arrays_opened=False,classification='All cases previously disclosed; offline diagnostic, not fresh confirmation',
        new_physical_LU=0,new_complex_RHS=0,hashes='not computed',files=files,
        scientific_parameter_changes_after_score='No changes to this sealed batch; any further algorithm is a separate version',
        prepared_by='experiments/prepare_stage_i_scoring.py; fresh isolated reproduction only')
    receipt=dict(status='PREPARED_NOT_SCORED',classification='OFFLINE_DIAGNOSTIC_ALL_CASES_PUBLIC_DEVELOPMENT_NOT_BLIND',
        model_assignment=None,scoring_executor='standalone reproduction script; no agent model inferred',seal_path=str(seal_path),seal_utc=now,
        prepared_utc=now,first_GT_read_utc=None,truth_reads=[],physical_LU=0,physical_complex_RHS=0,
        online_use_of_truth='FORBIDDEN; no parameter/selection changes to sealed batch',
        preparation_truth_access='NONE',preparation_physics_imported=False)
    audit.mkdir(parents=True,exist_ok=True)
    for path,value in ((seal_path,seal),(receipt_path,receipt)):
        with path.open('x') as stream:json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')
    return dict(status='PREPARED_NOT_SCORED',sealed_files=len(files),seal_path=str(seal_path),receipt_path=str(receipt_path),truth_access='NONE')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT,help='R5 directory in a fresh isolated copy; defaults to this script parent')
    print(json.dumps(prepare(parser.parse_args().root)))
