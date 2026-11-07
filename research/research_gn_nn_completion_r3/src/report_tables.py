"""Generate report tables from frozen TEST scores without changing selections."""
from pathlib import Path
import argparse,json
import numpy as np

ARMS=['v_image','v_physics','j_full','j_null','j_null_current']
LABELS={'continue_trf':'继续TRF（10秒）','weakreg45_anchor':'独立弱正则TRF（45秒）',
        'v_image':'V-image','v_physics':'V-physics','j_full':'J-full',
        'j_null':'J-null','j_null_current':'J-null-current'}

def avg(x):
    x=[float(v)for v in x if v is not None]
    return float(np.mean(x))if x else None
def f(x,d=4):return '—'if x is None else f'{x:.{d}f}'

def main():
    p=argparse.ArgumentParser();p.add_argument('--scores',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    scores=json.loads(a.scores.read_text());rows=scores['rows']
    runs=['continue_trf','weakreg45_anchor']+[f'confirm1024_{x}_seed7'for x in ARMS]
    names=[LABELS[x]for x in ['continue_trf','weakreg45_anchor']+ARMS]
    lines=['# 冻结确认结果表','','本表只汇总已经封存并离线评分的预测。材料均值基于有在线top1的对象，缺失数量及联合门另列；不将best@K当作top1。NRMSE=RMSE/材料增量。','']
    for level,noise in [('low',0.),('low',.001),('stress',0.),('stress',.001)]:
        title=('低对比'if level=='low'else'A=2.5外推')+(' clean'if noise==0 else' noisy')
        lines+=['## '+title,'','| 方法 | 在线合格/总数 | 联合成功/总数 | top1 NRMSE | IoU | 原NN NRMSE | 合格best@K NRMSE | 旧排序 NRMSE | 后处理秒 |','|---|---:|---:|---:|---:|---:|---:|---:|---:|']
        rr=[r for r in rows if r['run']==runs[2]and r.get('level','low')==level and r['noise']==noise]
        anchors=[next(c for c in r['candidates']if c['kind']=='anchor')for r in rr]
        lines.append(f"| GN输入 | {sum(c['physical_data_feasible']for c in anchors)}/{len(anchors)} | {sum(c['joint_success']for c in anchors)}/{len(anchors)} | {f(avg(c['NRMSE']for c in anchors))} | {f(avg(c['IoU_Otsu']for c in anchors),3)} | — | — | — | 前置费用单列 |")
        for run,name in zip(runs,names):
            rr=[r for r in rows if r['run']==run and r.get('level','low')==level and r['noise']==noise]
            top=[r['top1']for r in rr if r.get('top1')]
            raw=[next((c for c in r.get('candidates',[])if c['label']=='raw_1'),{})for r in rr]
            value=lambda key:avg((r.get(key)or{}).get('NRMSE')for r in rr)
            lines.append(f"| {name} | {len(top)}/{len(rr)} | {sum(r.get('joint_success',False)for r in rr)}/{len(rr)} | {f(avg(c['NRMSE']for c in top))} | {f(avg(c['IoU_Otsu']for c in top),3)} | {f(avg(c.get('NRMSE')for c in raw))} | {f(value('qualified_bestK_OFFLINE'))} | {f(value('old_score_top1'))} | {f(avg(r.get('wall_seconds')for r in rr),3)} |")
        lines+=['','独立弱正则TRF的秒数包含自己的45秒上限求解及验收，其余列为共享GN输入之后的增量费用，不能把两者混作同一冷启动延迟。弱正则行的best@K包含用于比较的公共GN数组，仅为离线参考池，不代表该冷启动方法实际获得了一个免费GN候选；其在线top1只取自己的端点。','']
    lines+=['## 三种训练种子','','按每个种子完整报告；没有按TEST选择最优种子。每格为clean/noisy平均NRMSE，均为64个母体；联合成功另见结构化表。','','| 方法 | seed7 | seed17 | seed27 |','|---|---:|---:|---:|']
    for arm in ARMS:
        vals=[]
        for seed in [7,17,27]:
            ns=[]
            for noise in [0.,.001]:
                rr=[r for r in rows if r['run']==f'confirm1024_{arm}_seed{seed}'and r.get('level','low')=='low'and r['noise']==noise]
                ns.append(f(avg(r['top1']['NRMSE']for r in rr if r.get('top1'))))
            vals.append(' / '.join(ns))
        lines.append('| '+LABELS[arm]+' | '+' | '.join(vals)+' |')
    a.output.write_text('\n'.join(lines)+'\n')
    print(str(a.output))

if __name__=='__main__':main()
