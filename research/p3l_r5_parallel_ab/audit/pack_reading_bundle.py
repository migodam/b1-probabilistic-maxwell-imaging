"""Compact Pro handoff. Raw datasets/models remain in independent archives."""
from pathlib import Path
import json,zipfile
P=Path(__file__).resolve().parents[1]
GAUSSIAN=P.parents[3]
OUT=GAUSSIAN.parent/'research/delegated/r5_delivery_20261026'

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    dest=OUT/'Research_P3LR5_Codex_Results_and_Pro_Handoff.zip'
    files=[]
    for f in sorted(P.rglob('*')):
        if not f.is_file() or '__pycache__' in f.parts:continue
        rel=f.relative_to(P)
        keep=(f.suffix in ('.md','.py','.png') or
              (len(rel.parts)==1 and f.suffix=='.json') or
              (rel.parts[0]=='audit' and f.suffix=='.json' and f.stat().st_size<200000) or
              f.name in ('OBSERVATION_SUMMARY.csv','RESULTS.json','METRICS.csv','BOOTSTRAP_EXPLORATORY.json',
                         'OPTIMIZED_COST_BENCHMARK.json','MOMENT_VALIDATION.json','JOINT_CDF_VALIDATION.json',
                         'SELECTION_FROZEN.json','B_STATUS.json','A_FINAL_QA.json','COST_LEDGER_FINAL.json'))
        if keep:files.append(f)
    note='R5 Codex actual results and next-Pro handoff. Start with START_HERE.md and GPT_PRO_HANDOFF/START_HERE.md. This compact ZIP excludes complete raw NPZ/PT/training data. Independent raw evidence archives and the repository contain all original endpoints, checkpoints and failure records. No claim of FULL posterior or new locked generalization.\n'
    with zipfile.ZipFile(dest,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        z.writestr('R5_Codex_Results/READ_ME_FIRST.txt',note)
        for f in files:z.write(f,'R5_Codex_Results/'+str(f.relative_to(P)))
        z.writestr('R5_Codex_Results/READING_MANIFEST.json',json.dumps(dict(
            raw_arrays_included=False,files=[dict(path=str(f.relative_to(P)),bytes=f.stat().st_size) for f in files],
            scientific_gate='A REDIRECT; B CONDITIONAL DEVELOPMENT GO; independent confirmation NOT_RUN'),indent=2))
    with zipfile.ZipFile(dest) as z:
        assert z.testzip() is None
    print(json.dumps(dict(path=str(dest),files=len(files),bytes=dest.stat().st_size,raw_arrays=False)))
if __name__=='__main__':main()
