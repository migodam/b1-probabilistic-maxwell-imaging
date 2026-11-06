"""Render-only correction: fresh Twofold actual rank30, historical DEV28 unchanged."""
from pathlib import Path
import sys,json,csv,time
import numpy as np
P=Path(__file__).resolve().parent;ROOT=P.parents[2];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r10.experiments.score_twod_comparison import figures
from research.research_max_struct_r10.experiments.analyze_twod_energy import paper_figures
def main():
    start=time.perf_counter();cpu=time.process_time();tmp=ROOT/'delegated/twod_score/label_refresh_rank30';tmp.mkdir(parents=True,exist_ok=False)
    rows=json.loads((P/'SCORES_FINAL.json').read_text());truth={};pics={}
    for r in rows:
        case=r['case_id']
        if case not in truth:
            with np.load(ROOT/'delegated/twod_energy/gt_current_cache'/(case+'.npz'),allow_pickle=False) as z:truth[case]=z['m_gt'].copy()
        with np.load(ROOT/'RESULTS/twod'/r['task_id']/'prediction.npz',allow_pickle=False) as z:pics[case,r['method']]=z['m'].copy()
    figures(tmp,tmp/'material_figures',rows,truth,pics)
    ep=ROOT/'delegated/twod_energy/scored'
    with (ep/'SINGLE_OPERATOR_ENERGY.csv').open() as f:single=list(csv.DictReader(f))
    for r in single:
        for k in ('tau','component_energy','total_energy'):r[k]=float(r[k])
    with (ep/'ORACLE_CAPACITY_BY_CASE_RANK_MF.csv').open() as f:oracle=list(csv.DictReader(f))
    for r in oracle:
        for k in ('mf','receiver_rank','actual_coefficient_range_rank'):r[k]=int(r[k])
        r['affine_GT_irreducible_squared_energy']=float(r['affine_GT_irreducible_squared_energy'])
    paper_figures(tmp/'energy_figures',single,[],oracle)
    dest=ROOT/'FIGURES/twod_comparison';archive=tmp/'previous_rank28_label_archive';archive.mkdir()
    changed=[]
    for folder in ('material_figures','energy_figures'):
        for src in sorted((tmp/folder).iterdir()):
            dst=dest/src.name
            if dst.exists():dst.replace(archive/dst.name)
            src.replace(dst);changed.append(str(dst))
    (tmp/'RECEIPT.json').write_text(json.dumps(dict(status='LABEL_ONLY_RENDER_CORRECTION',reason='Fresh actualTwofold rank30; previous label inherited historicalDEV28',
        rank_source='Fresh DENOMINATOR_SETUP and counted smallLU dimension30',historical_DEV_labels_changed=False,
        scoring_values_changed=False,methods_thresholds_budgets_changed=False,new_GT_forward_LU=0,new_inverse_runs=0,
        truth_images_source='Previously authorized offline GT-current cache m_gt, no new sealed_truth reads',files=changed,
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start),indent=2)+'\n')
    print(json.dumps(dict(status='COMPLETE_LABEL_ONLY',files=len(changed))))
if __name__=='__main__':main()
