"""Freeze V7 guarded portfolios before four entirely new confirm targets."""
import datetime
import json
from pathlib import Path
import shutil
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from src.portfolio import PORTFOLIO_DEFAULTS


def main():
    target=ROOT/'experiments'/'CONFIRM_V7_CONTRACT.json'
    destination=ROOT/'audit'/'CONFIRM_V7_SOURCE_SNAPSHOT'
    if target.exists() or destination.exists():raise FileExistsError('Preserving V7 freeze')
    old=json.loads((ROOT/'audit'/'HOLDOUT_SOURCE_FREEZE.json').read_text())
    for item in old['files']:
        if (ROOT/item['source']).read_bytes()!=(ROOT/item['snapshot']).read_bytes():
            raise RuntimeError('V6 core or original CLI changed: '+item['source'])
    frozen=datetime.datetime.now(datetime.timezone.utc).isoformat()
    methods=[dict(name='GuardedGlobalDual',method='guarded_globaldual',config=dict(PORTFOLIO_DEFAULTS)),
             dict(name='GuardedReciprocal',method='guarded_reciprocal',config=dict(PORTFOLIO_DEFAULTS))]
    contract=json.loads((ROOT/'experiments'/'DEV12_CONTRACT.json').read_text())
    contract.update(schema='Z-ALG-R3-CONFIRM-V7-V1',frozen_utc=frozen,split='confirm',case_prefix='NEWCONFIRM12',
                    confirm_count=4,method_configs=methods,methods=[m['name'] for m in methods],
                    solver_revision='R3_PUBLIC_V7_GUARDED_PORTFOLIO',runner='run_portfolio.py',
                    source_manifest='audit/CONFIRM_V7_SOURCE_FREEZE.json',source_snapshot='audit/CONFIRM_V7_SOURCE_SNAPSHOT',
                    role='Fresh confirmatory targets after V7 revision freeze. Exposed V6 holdout results are development for V7, preserved unchanged.',
                    prior_development_evidence='results/HOLDOUT_V6_SUMMARY.json',
                    common_guard='One native-material FullJ low-constant public start, max_nfev80, cap0.20s, residual goal1e-10; retain every guard candidate',
                    continuation='Only remaining measured total3s budget; V6GlobalDual exactDEVconfig versus V6reciprocalFullJ1024',
                    total_wall_boundary='Outer pipeline timer includes public load, common guard, repeated setup, continuation and wrapper bookkeeping; actual overrun reported',
                    candidate_selection='Minimum original all-observation residual across all guard and continuation candidates; no GT',
                    source_checks='All eight V6 source/CLI files byte-identical to their existing immutable freeze; new portfolio wrapper and CLI frozen separately',
                    primary_comparison='GuardedGlobalDual versus GuardedReciprocal after same common guard',
                    data_accepted_relative_residual=1e-6,independent_validation='all returned candidates',
                    metrics=['Residual-selected top1 material RMS','All-data relative residual','Offline oracle best@K RMS','Full/reduced/inner-LS costs','Actual outer wall and overrun','Guard versus continuation phase costs'],
                    original_geometry_and_prior_changed=False,
                    prohibitions=['No V7 tuning after new-confirm GT scoring','No reuse of exposed V6 cases as V7 blind evidence','No budget extension or silent retries','No N24 expansion in this contract'])
    target.write_text(json.dumps(contract,indent=2)+'\n');target.chmod(0o444)
    sources=sorted((ROOT/'src').glob('*.py'))+[ROOT/'run_solver.py',ROOT/'run_portfolio.py',ROOT/'experiments'/'prepare_data.py',target]
    files=[]
    for source in sources:
        copy=destination/source.relative_to(ROOT);copy.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(source,copy);copy.chmod(0o444)
        if source.read_bytes()!=copy.read_bytes():raise RuntimeError('Source snapshot mismatch')
        stat=copy.stat();files.append(dict(source=str(source.relative_to(ROOT)),snapshot=str(copy.relative_to(ROOT)),bytes=stat.st_size,mtime_ns=stat.st_mtime_ns))
    receipt=ROOT/'audit'/'CONFIRM_V7_SOURCE_FREEZE.json'
    receipt.write_text(json.dumps(dict(frozen_utc=frozen,files=files,method_contract=str(target.relative_to(ROOT)),
                                       V6_sources_unchanged=True,solver_confirmed_final=True,
                                       new_confirm_targets_generated=False,truth_read=False,
                                       identity='Immutable source/config snapshots byte-for-byte compared; no SHA checks'),indent=2)+'\n')
    receipt.chmod(0o444)
    print(json.dumps(dict(frozen_utc=frozen,source_files=len(files),methods=[x['name'] for x in methods],new_confirm_targets_generated=False)))


if __name__=='__main__':main()
