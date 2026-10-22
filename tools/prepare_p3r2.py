"""Freeze an independent authorization; never runs physics or edits history."""
from pathlib import Path
import sys,json,shutil,subprocess
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r2.core import OUT,write,CAPS

def main():
    if (OUT/'configs/AUTHORIZATION.json').exists():raise FileExistsError('authorization already frozen')
    write(OUT/'configs/AUTHORIZATION.json',dict(campaign_id='P3LR2-20261022',explicit_user_authorization=True,authority='User Research P3L-R2 full campaign prompt, 2026-10-22 Asia/Singapore',RHS_caps=CAPS,total_RHS=20000000,controlled_CPU_seconds=86400,GPU_seconds=86400,paid_cloud=False,new_3D=False,new_Mie=False,old_balances_used=False))
    write(OUT/'configs/PROTOCOL.json',dict(target_id='P1-CORRELATED-HOST-MIXTURE-TM24-v1-20261018',master_seed=2026102202,independent_tracks=['A GT-supervised material MSE; no posterior teacher required','B exact-corrected nonGaussian inference; SLICE separately identified','C material-native trajectory and supervised risk policy; no reference prerequisite'],training_sizes=[4096,16384],optional_training_size=32768,validation_parents=2048,fresh_locked_parents=256,old_locked_status='EXPOSED; never called fresh',primary_baseline='native96 plus incremental probit64',primary_metric='parent mean full-grid complex material MSE',primary_comparison='validation-selected risk estimator versus strong probit point',bootstrap=dict(draws=2000,seed=2026102292),positive_research_threshold=dict(relative_risk_reduction=.05,paired_CI_excludes_zero=True),GT_isolation='train supervision permitted; online physics never sees GT; fresh locked truth opened only by separate offline evaluator after allmethods frozen',memory=dict(paid_workers=1,transfer_entries=4,RSS_stop_MiB=2000,numerical_threads=1),target_labels=['FULL','SLICE','LOCAL-C','PROFILE','SUPERVISED_RISK','CANDIDATE_ONLY'],unqualified_full_does_not_block_risk_or_geometry=True))
    old=ROOT.parent/'b1-p3l-local-posterior/research/p3l_local_posterior_imaging'
    for i in range(32):
        parent=f'DEV_{i:03d}'
        for src,dest in [(old/f'inputs/online/{parent}.npz',OUT/f'inputs/online/{parent}.npz'),(old/f'A/anchors/{parent}.npz',OUT/f'A/anchors/{parent}.npz')]:
            dest.parent.mkdir(parents=True,exist_ok=True);dest.symlink_to(src)
    write(OUT/'AUDIT/HISTORY_REUSE.json',dict(source_commit='c1171ab9c9006f5061f10293c7eadd5b977dd608',old_sealed_RHS=8374472,new_RHS=0,old_anchor_observation_cache_parent_ids=[f'DEV_{i:03d}' for i in range(32)],source_root=str(old),GT_used=False,reused_background='paid P1 G/S/incident, copied weak.npz',physics_setup_from_zero_not_free=True))
    (OUT/'START_HERE.md').write_text('# P3L-R2：材料风险、非 Gaussian 推断与几何轨迹\n\n状态：RUNNING。新独立授权 20,000,000 complex RHS / 24h CPU / 24h GPU 上限；原 P1 不变，旧预算封账。A/B/C 独立推进，未合格参考不阻断监督风险或轨迹。\n\n当前不是完成报告；结果按阶段追加。\n')
    write(OUT/'SOURCE_INVENTORY.json',dict(base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),old_data_root=str(old),frozen_A30_commit='717d3d26a290248bd2b9aeb55f3d03694b407b84',historical_ledgers_modified=False))
    print('R2 authorization and protocol frozen; zero new physics')

if __name__=='__main__':main()
