"""Render the publication tables from sealed paired 2D outputs."""
from paired_stats import ROOT,load,paired
import numpy as np
import json

def f(x):return '—' if x is None else f'{float(x):.3f}'
def pct(x):return f'{100*x:+.1f}%'
def tab(head,rows):
    return '\n'.join(['| '+' | '.join(head)+' |','|'+'|'.join(['---']*len(head))+'|']+['| '+' | '.join(map(str,r))+' |' for r in rows])
def summary(cohort, arms=('dual','srap')):
    r=load(cohort+'_metrics.csv');lines=[]
    for T in ('15','30'):
        for v in ('raw','light'):
            out=[]
            for array in sorted({q['array'] for q in r if q['suite']=='LA_SYNTH'}):
                ss=[z for z in r if z['budget']==T and z['view']==v and z['array']==array]
                z=paired(ss,*arms)
                if z.get('eligible',0)==0:continue
                ci=' / '.join(pct(x) for x in z['ci95'])
                out.append([array,z['eligible'],f(z['baseline_median']),f(z['srap_median']),pct(z['median_gain']),ci,
                            f"{z['wins']}/{z['ties']}/{z['losses']}",f"{z['dual_better_air']}/{z['srap_better_air']}",
                            f(z['iou_dual'])+' / '+f(z['iou_srap']),f(z['background_dual'])+' / '+f(z['background_srap'])])
            lines.extend([f'### {T} 秒 · {v}','',tab(['阵列','配对数','基线 L2','SRAP L2','配对中位降幅','95%材料 bootstrap','胜/平/负','优于空气 基线/SRAP','IoU 基线/SRAP','背景 L2 基线/SRAP'],out),''])
    return lines

def real_table(fname):
    r=load(fname)
    mains=[q for q in r if q['budget']=='30' and q['suite']=='LA_FRESNEL']
    methods=sorted({q['method'] for q in mains})
    rows=[]
    for case in sorted({q['case'] for q in mains}):
        by={(q['method'],q['view']):q for q in mains if q['case']==case}
        for key in ('held_tx','held_rx'):
            if not any(q.get(key) for q in by.values()):continue
            rows.append([case, key]+[' / '.join(f((by.get((m,v)) or {}).get(key) or None) for v in ('raw','light')) for m in methods])
    return tab(['设置','留出通道']+methods,rows)

def main():
    out=['# 逐阵列统计与配对表','','表中“胜”按材料相对 L2 下降计数。相对空气参照采用 L2 < 1；它仅衡量连续材料误差。每行同一批对象在同一阵列的两个算法配对。95% 区间用材料对象成对重采样，未校正多阵列探索。不把 7 次布阵当作 7 件材料。','']
    for nm,label in [('confirmation','24个冻结新材料'),('followup','8个连续材料与窄桥补集')]:
        if (ROOT/'results'/(nm+'_metrics.csv')).is_file():
            out+=['## '+label,'']+summary(nm)
            if nm=='confirmation':out+=['## 直线阵及同角弧阵的共享 DPCS 截断','']+summary(nm,('dual_cap','srap_cap'))
            if nm=='followup':out+=['## 两个120°角窗的共享DPCS截断','']+summary(nm,('dual_cap','srap_cap'))
    for nm,title in [('legacy_metrics.csv','XINAN 半环实测重放'),('measured_extra_metrics.csv','XINAN 通道数配平实测'),('fresnel_lazy_metrics.csv','Mac 同预算端口延迟准备')]:
        if (ROOT/'results'/nm).is_file():out+=['## '+title,'','每格依次为 raw / 轻校正的留出复场相对误差，越低越好。','',real_table(nm),'']
    (ROOT/'PAIRED_TABLES_ZH.md').write_text('\n'.join(out)+'\n')
if __name__=='__main__':main()
