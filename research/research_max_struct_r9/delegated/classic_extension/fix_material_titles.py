"""Title-only revision from stored arrays; no scoring or scientific calls."""
import time
START,CPU_START=time.perf_counter(),time.process_time()
from pathlib import Path
import json
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
REPO=Path(__file__).resolve().parents[4]
ROOT=REPO/'research/research_max_struct_r9'
OUT=ROOT/'RESULTS/classic_extension_xinan'
FIG=ROOT/'FIGURES/classic_extension'
ARMS=['cie_single','cie_dual','global_dual','a_material','a_plus','physical_lbfgs','paper_nie','paper_lsie']


def main():
    receipt_path=OUT/'PLOT_TITLE_FIX_RECEIPT.json'
    if receipt_path.exists():raise FileExistsError('Retain the previous title revision')
    manifest=json.loads((ROOT/'configs/classic_extension_main.json').read_text())
    paths=[FIG/f'{c["case_id"]}_{scale}.png' for c in manifest['cases'] for scale in ['full_range','morphology']]
    archive=FIG/'attempt1';archive.mkdir(exist_ok=True)
    for path in paths:
        backup=archive/path.name
        if backup.exists():raise FileExistsError('Retain archived original '+str(backup))
    for path in paths:shutil.copy2(path,archive/path.name)
    for case in manifest['cases']:
        truth=np.load(ROOT/'sealed_truth/classic_extension'/(case['case_id']+'.npz'))['m']
        panels=[('Truth',truth)]+[(arm,np.load(OUT/f'{case["case_id"]}__{arm}'/'PREDICTION.npz')['m']) for arm in ARMS]
        for scale,top in [('full_range',20.),('morphology',case['epsilon_r']-1)]:
            fig,axes=plt.subplots(3,3,figsize=(12,12),constrained_layout=True)
            for ax,(name,m) in zip(axes.flat,panels):
                im=ax.imshow(1+m.reshape(64,64),origin='lower',extent=(-1,1,-1,1),
                             vmin=1,vmax=1+top,cmap='viridis')
                ax.set_title(f'{name}\nmax εr={1+m.max():.3g}')
                ax.set_xlabel('x (m)');ax.set_ylabel('y (m)')
            fig.colorbar(im,ax=list(axes.flat),label='εr')
            if scale=='full_range':
                fig.suptitle(f'{case["case_id"]}: full material box, common εr scale 1–21\n'
                    'Complete declared range; max labels show amplitudes')
            else:
                fig.suptitle(f'{case["case_id"]}: morphology, CLIPPED at truth amplitude (εr≤{case["epsilon_r"]:.3g})\n'
                    'Shape detail only; amplitudes shown by max labels and full-range companion')
            fig.savefig(FIG/f'{case["case_id"]}_{scale}.png',dpi=150)
            plt.close(fig)
    cpu=time.process_time()-CPU_START
    receipt=dict(status='COMPLETED_TITLE_ONLY_REVISION_2',cpu_seconds=cpu,
        wall_seconds=time.perf_counter()-START,physical_solves=0,scoring_runs=0,
        source='existing sealed prediction arrays and existing offline truth arrays only',
        redrawn=[str(p.relative_to(ROOT)) for p in paths],
        archived_attempt1=[str((archive/p.name).relative_to(ROOT)) for p in paths],
        remaining_five_figures='UNCHANGED; not redrawn',
        scope='includes imports/archive/read/render through this point; final serialization tail unharvested')
    receipt_path.write_text(json.dumps(receipt,indent=2))
    ledger_path=OUT/'PHASE2_COST_LEDGER.json';d=json.loads(ledger_path.read_text())
    old_total=d['actual_new_known_CPU_seconds'];old_mac=d['by_platform']['Mac']['known_CPU_seconds']
    d['actual_new_known_CPU_seconds']+=cpu
    d['by_platform']['Mac']['known_CPU_seconds']+=cpu
    d['by_platform']['Mac']['analysis_preparation_plotting_CPU_seconds']+=cpu
    d['primary_cohort_plus_shared_analysis_known_CPU_seconds']+=cpu
    d['overhead'].append(dict(path=str(receipt_path),**receipt))
    d['plot_title_revision']=2;d['finalized_for_parent_compiler']=True
    ledger_path.write_text(json.dumps(d,indent=2,allow_nan=False))
    report=ROOT/'CLASSIC_BENCHMARKS.md';text=report.read_text()
    text=text.replace(f'{old_total:.6f}',f'{d["actual_new_known_CPU_seconds"]:.6f}')
    text=text.replace(f'{old_mac:.6f}',f'{d["by_platform"]["Mac"]["known_CPU_seconds"]:.6f}')
    text+='\n图题说明已作title-only revision2：两例fullrange明确完整εr=1–21 box；morphology明确按truth幅度截断，仅供形态细节，实际幅度见max标签和fullrange图。首轮4张材料图保留在FIGURES/classic_extension/attempt1；其余5图未重画，预测、评分、预算全部未改。本次实际绘图CPU已追加canonical。\n'
    report.write_text(text)
    print(json.dumps(dict(status=receipt['status'],CPU_seconds=cpu,
        actual_new_known_CPU_seconds=d['actual_new_known_CPU_seconds'],Mac_known_CPU_seconds=d['by_platform']['Mac']['known_CPU_seconds'])),flush=True)


if __name__=='__main__':main()
