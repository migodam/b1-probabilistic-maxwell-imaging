"""Freeze five preselected DEV method configs and source before new targets."""
import datetime
import json
import shutil
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def main():
    destination=ROOT/'audit'/'HOLDOUT_SOURCE_SNAPSHOT'
    contract_path=ROOT/'experiments'/'HOLDOUT12_CONTRACT.json'
    if destination.exists() or contract_path.exists():raise FileExistsError('Preserving method/source freeze')
    frozen=datetime.datetime.now(datetime.timezone.utc).isoformat()
    sources=sorted((ROOT/'src').glob('*.py'))+[ROOT/'run_solver.py']
    source_rows=[]
    for source in sources:
        target=destination/source.relative_to(ROOT)
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(source,target);target.chmod(0o444)
        if target.read_bytes()!=source.read_bytes():raise RuntimeError('Snapshot mismatch')
        stat=target.stat()
        source_rows.append(dict(source=str(source.relative_to(ROOT)),snapshot=str(target.relative_to(ROOT)),bytes=stat.st_size,mtime_ns=stat.st_mtime_ns))
    mapping=[('FullJ1024','fullj','dev_mechanism_DEV12_00_FullJ1024.json'),
             ('FullJReciprocal1024','fullj','dev_mechanism_DEV12_00_FullJReciprocal1024.json'),
             ('Coordinate','coordinate','dev_reduced_strong_DEV12_00_Coordinate.json'),
             ('ABAR4_coordpre','algebraic','dev_mechanism_DEV12_00_ABAR4_coordpre.json'),
             ('GlobalDual','globaldual','dev_mechanism_DEV12_00_GlobalDual.json')]
    defaults=json.loads((ROOT/'results'/mapping[-1][2]).read_text())['config']
    methods=[]
    for name,method,receipt in mapping:
        config=dict(defaults)
        config.update(json.loads((ROOT/'results'/receipt).read_text())['config'])
        methods.append(dict(name=name,method=method,config=config,DEV_config_source='results/'+receipt))
    contract=json.loads((ROOT/'experiments'/'DEV12_CONTRACT.json').read_text())
    contract.update(schema='Z-ALG-R3-HOLDOUT12-V1',frozen_utc=frozen,
                    role='Four fresh confirmatory targets after method/source freeze; no target-dependent retuning',
                    holdout_count=4,method_configs=methods,methods=[x['name'] for x in methods],
                    solver_revision='R3_PUBLIC_V6_GLOBALDUAL_RECIPROCAL_CONTROL',
                    source_snapshot='audit/HOLDOUT_SOURCE_SNAPSHOT',source_manifest='audit/HOLDOUT_SOURCE_FREEZE.json',
                    primary_method='GlobalDual',primary_baseline='FullJ1024',
                    data_accepted_relative_residual=1e-6,independent_validation='all returned candidates',
                    original_geometry_and_prior_changed=False,
                    budget_boundary='Maximum 3 s numerical wall per method/case; terminal atomic overshoot is measured; no retry or retune',
                    reporting='Residual-selected top1 RMS primary; offline oracle best@K separately labeled; all costs retained')
    contract_path.write_text(json.dumps(contract,indent=2)+'\n');contract_path.chmod(0o444)
    manifest=dict(frozen_utc=frozen,solver_confirmed_final=True,files=source_rows,
                  identity='Read-only source snapshots byte-for-byte compared at freeze; no SHA checks',
                  target_truth_read=False,holdout_targets_generated=False,
                  method_contract=str(contract_path.relative_to(ROOT)))
    path=ROOT/'audit'/'HOLDOUT_SOURCE_FREEZE.json'
    path.write_text(json.dumps(manifest,indent=2)+'\n');path.chmod(0o444)
    print(json.dumps(dict(frozen_utc=frozen,methods=[x['name'] for x in methods],source_files=len(sources),holdout_targets_generated=False)))


if __name__=='__main__':main()
