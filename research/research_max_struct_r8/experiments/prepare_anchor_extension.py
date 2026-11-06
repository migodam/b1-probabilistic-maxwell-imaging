"""Freeze the affine-current translation extension without reading test truth."""
from pathlib import Path
from datetime import datetime,timezone
import json,copy,zipfile,shutil
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]

def main():
    stamp=datetime.now(timezone.utc).isoformat();cfg=json.loads((ROOT/'configs/warm_anchor_development.json').read_text())
    cfg.update(run_id='r8_affine_anchor_extension',freeze_utc=stamp,status='FROZEN_BEFORE_CONFIRMATION_TRUTH',development_selection='All42 unscored objects; no geometry or outcome exclusion; unchanged original network and solver parameters')
    cpath=ROOT/'configs/warm_anchor_extension.json';cpath.write_text(json.dumps(cfg,indent=2))
    source=json.loads((ROOT/'configs/warm_all_manifest.json').read_text())['tasks'];tasks=[]
    for t in source:
        if t['initialization']=='nn' and t['arm'] in ('cie_dual','a_plus','a_plus_sourcewise'):
            n=copy.deepcopy(t);n['initialization']='nn_affine_anchor';n['task_id']=n['case_id']+'__'+n['arm']+'__nn_affine_anchor'
            n['output_dir']='research/research_max_struct_r8/RESULTS/warm_anchor_confirmation/'+n['task_id'];n['initialize_physical_current']=True;n['config']=str(cpath.relative_to(REPO));tasks.append(n)
    remote=set(json.loads((ROOT/'audit/HARDWARE_ASSIGNMENT.json').read_text())['remote_objects'])
    for name,rows in [('confirmation',tasks),('local',[t for t in tasks if t['case_id'] not in remote]),('xinan',[t for t in tasks if t['case_id'] in remote])]:
        (ROOT/f'configs/warm_anchor_{name}_manifest.json').write_text(json.dumps(dict(config=str(cpath.relative_to(REPO)),tasks=rows),indent=2))
    dev=json.loads((ROOT/'configs/warm_development_manifest.json').read_text())['tasks'];rows=[]
    for t in dev:
        if int(t['case_id'].split('_')[-1]) in (0,3):continue
        if (t['initialization']=='nn' and t['arm'] in ('cie_dual','a_plus')) or (t['initialization']=='uniform' and t['arm']=='a_plus'):
            n=copy.deepcopy(t);n['initialization']+='_affine_anchor';n['task_id']=n['case_id']+'__'+n['arm']+'__'+n['initialization'];n['initialize_physical_current']=True;n['config']=str(cpath.relative_to(REPO))
            n['output_dir']='research/research_max_struct_r8/RESULTS/warm_anchor_development_extended/'+n['task_id'];rows.append(n)
            if t['initialization']=='nn' and t['arm']=='a_plus':
                n=copy.deepcopy(n);n['arm']='a_plus_sourcewise';n['task_id']=n['case_id']+'__'+n['arm']+'__'+n['initialization'];n['output_dir']='research/research_max_struct_r8/RESULTS/warm_anchor_development_extended/'+n['task_id'];rows.append(n)
    (ROOT/'configs/warm_anchor_development_extended_manifest.json').write_text(json.dumps(dict(config=str(cpath.relative_to(REPO)),tasks=rows),indent=2))
    freeze=ROOT/'audit/ANCHOR_EXTENSION_FREEZE';freeze.mkdir(exist_ok=False)
    for p in [ROOT/'src/cie_fft.py',ROOT/'experiments/run_case.py',cpath,ROOT/'configs/warm_anchor_confirmation_manifest.json',ROOT/'configs/warm_anchor_development_extended_manifest.json']:shutil.copy2(p,freeze/p.name)
    (freeze/'PROTOCOL.json').write_text(json.dumps(dict(frozen_utc=stamp,confirmation_tasks=len(tasks),confirmation_objects=42,development_extension_tasks=len(rows),truth_access='NONE',rationale='Pilot validation0 preserves recovery; validation3 recovers under all three anchored CIE arms while uniform anchor fails. Test transfer across all42 before any testGT is read.',new_mathematics='None claimed: fixed affine APCS translation retaining initial physical current; changes early truncated family, no added optimized DOF.',comparisons='Same-object same-hardware unanchored CIE/Dual/A+ and NN direct physical; all original750 remain unchanged.',same_current_budget=True,additional_physics_LU_for_anchor=0,additional_fixed_complex_storage='576 times source count; initial projection charged',no_test_adaptation=True),indent=2))
    # Self-contained public bundle in a separate remote directory; never change
    # the running primary checkout, and never transfer validation/test labels.
    paths={p for p in REPO.glob('research/**/*.py') if not {'__pycache__','audit','delegated'}.intersection(p.relative_to(REPO).parts)}
    for t in tasks:
        if t['case_id'] in remote:
            paths.update(REPO/p for p in t['public_files']);paths.add(REPO/t['initialization_file'])
    paths.update([cpath,ROOT/'configs/solver.json',ROOT/'configs/warm_anchor_xinan_manifest.json'])
    archive=ROOT/'audit/transfer_anchor_xinan.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=3) as z:
        for p in sorted(paths):z.write(p,p.relative_to(REPO))
    print(json.dumps(dict(confirmation=len(tasks),remote=sum(t['case_id'] in remote for t in tasks),development=len(rows),archive_bytes=archive.stat().st_size)))

if __name__=='__main__':main()
