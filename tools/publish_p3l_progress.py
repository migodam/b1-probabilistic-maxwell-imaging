"""Render saved P3-L progress locally in the existing B1 publication path.

No physics, truth arrays, training, ledger-journal scan, archive, git or push.
Only explicitly named saved summaries and small production-Phi trace members
are read. Existing research prose and parent-authored decisions are preserved.
The page is a RUNNING checkpoint, never an automatic final scientific verdict.
"""
from pathlib import Path
from datetime import datetime,timezone
import ast
import html
import json
import math
import re
import zipfile

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'research/p3l_local_posterior_imaging'
REL='research/p3l_local_posterior_imaging/'


def read_json(path):
    path=Path(path)
    if path.stat().st_size>24*1024*1024:
        raise RuntimeError(f'unexpectedly large summary; manual inspection required: {path.name}')
    return json.loads(path.read_text())


def fmt(value,digits=5):
    return 'UNASSESSED' if value is None or not math.isfinite(float(value)) else f'{value:.{digits}g}'


def link(path,label=None):
    return '<a href="'+html.escape(str(path),quote=True)+'">'+html.escape(label or str(path))+'</a>'


def inline(text):
    text=html.escape(text)
    text=re.sub(r'\[([^\]]+)\]\(([^)]+)\)',r'<a href="\2">\1</a>',text)
    text=re.sub(r'\*\*([^*]+)\*\*',r'<strong>\1</strong>',text)
    return re.sub(r'`([^`]+)`',r'<code>\1</code>',text)


def render_markdown(text):
    # Match the existing B1 renderer's plain prose/tables and add real lists.
    lines=text.splitlines();output=[];i=0
    while i<len(lines):
        line=lines[i]
        if not line.strip():i+=1;continue
        if line.startswith('|'):
            rows=[]
            while i<len(lines) and lines[i].startswith('|'):
                cells=lines[i].strip('|').split('|');i+=1
                if not all(re.fullmatch(r'\s*:?-+:?\s*',cell) for cell in cells):rows.append(cells)
            output.append('<div class="scroll"><table>'+''.join('<tr>'+''.join('<'+('th' if n==0 else 'td')+'>'+inline(cell.strip())+'</'+('th' if n==0 else 'td')+'>' for cell in row)+'</tr>' for n,row in enumerate(rows))+'</table></div>')
            continue
        if line.startswith('#'):
            level=min(6,len(line)-len(line.lstrip('#')))
            output.append(f'<h{level}>'+inline(line.lstrip('#').strip())+f'</h{level}>');i+=1;continue
        if line.startswith('- '):
            items=[]
            while i<len(lines) and lines[i].startswith('- '):items.append('<li>'+inline(lines[i][2:])+'</li>');i+=1
            output.append('<ul>'+''.join(items)+'</ul>');continue
        output.append('<p>'+inline(line)+'</p>');i+=1
    return '\n'.join(output)


def existing_style():
    # Read the existing B1 style literal without importing old publication code.
    tree=ast.parse((ROOT/'tools/publish_lxz.py').read_text())
    for node in tree.body:
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='STYLE' for t in node.targets):
            return ast.literal_eval(node.value)
    raise RuntimeError('existing B1 STYLE literal missing')


def save_figure(fig,name):
    folder=OUT/'figures';folder.mkdir(exist_ok=True)
    fig.savefig(folder/f'{name}.png',dpi=170,bbox_inches='tight')
    fig.savefig(folder/f'{name}.pdf',bbox_inches='tight')
    plt.close(fig)


def references():
    rows=[]
    for stage in ('A','C'):
        for path in sorted((OUT/stage/'references').glob('*/*/REFERENCE_STATS.json')):
            stats=read_json(path)
            if stats.get('target_family')!='FULL':continue
            material=stats.get('summaries',{}).get('material_functions',{})
            rows.append(dict(path=path.relative_to(OUT).as_posix(),label=path.parent.parent.name+'/'+path.parent.name,
                stats=stats,rhat=material.get('maximum_rhat'),bulk=material.get('minimum_bulk_ess'),
                tail=material.get('minimum_tail_ess'),mcse=stats.get('full_grid_mean_rms_mcse')))
    return rows


def figures(refs,risk,geometry):
    names=[]
    if refs:
        labels=[row['label'] for row in refs]
        fig,axes=plt.subplots(1,3,figsize=(17,max(4,len(refs)*.55)),constrained_layout=True)
        for ax,key,threshold,title in zip(axes,['rhat','bulk','mcse'],[1.01,1000,.001],
                ['Material-function maximum Rhat','Minimum material bulk ESS','Full-grid mean RMS MCSE']):
            values=[np.nan if row[key] is None else row[key] for row in refs]
            ax.barh(range(len(refs)),values,color='#25667d');ax.axvline(threshold,color='#b64c3a',ls='--')
            ax.set_yticks(range(len(refs)),labels);ax.invert_yaxis();ax.set_title(title)
            if key!='rhat':ax.set_xscale('log')
            for i,value in enumerate(values):
                if not np.isfinite(value):ax.text(.05,i,'UNASSESSED',transform=ax.get_yaxis_transform(),color='#b64c3a')
        fig.suptitle('Saved FULL-reference diagnostics; thresholds alone do not prove global coverage')
        save_figure(fig,'p3l_progress_reference_diagnostics');names.append('p3l_progress_reference_diagnostics')
        fig,axes=plt.subplots(math.ceil(len(refs)/2),2,figsize=(14,3*math.ceil(len(refs)/2)),squeeze=False,constrained_layout=True)
        for ax,row in zip(axes.flat,refs):
            array_path=OUT/row['path'].replace('REFERENCE_STATS.json','REFERENCE_ARRAYS.npz')
            if array_path.exists():
                with zipfile.ZipFile(array_path) as archive:
                    if archive.getinfo('production_Phi.npy').file_size>32*1024*1024:
                        raise RuntimeError('unexpectedly large saved Phi member; no automatic trace loading')
                with np.load(array_path,allow_pickle=False) as data:
                    # NPZ is lazy: only Phi is decompressed, never material/q.
                    phi=data['production_Phi']
                if phi.ndim!=2 or phi.size>4*1000000:
                    raise RuntimeError('unexpected trace layout; no automatic large trace rendering')
                stride=max(1,phi.shape[1]//1500)
                for chain,values in enumerate(phi):ax.plot(np.arange(0,len(values),stride),values[::stride],lw=.7,label=f'chain {chain}')
                ax.set_title(row['label']);ax.set_xlabel('Production iteration (decimated display)');ax.set_ylabel('Full Phi');ax.legend(fontsize=7)
            else:ax.text(.1,.5,'Saved trace member unavailable',transform=ax.transAxes)
        for ax in list(axes.flat)[len(refs):]:ax.axis('off')
        fig.suptitle('Unqualified chains: saved likelihood traces, warmup already removed; not posterior samples for distillation')
        save_figure(fig,'p3l_progress_likelihood_traces');names.append('p3l_progress_likelihood_traces')
    risk_rows=risk['rows']
    methods=['optimized_anchor','gaussian_material_mean','gaussian_material_geometric_median_mc']
    labels=['Optimized anchor','GN Gaussian material mean','GN Gaussian geometric median']
    fig,axes=plt.subplots(1,2,figsize=(14,4.5),constrained_layout=True)
    for method,label in zip(methods,labels):
        values=[row['risks'][method]['full_grid']['complex_MSE'] for row in risk_rows]
        axes[0].plot(range(len(values)),values,'.-',lw=.7,label=label)
    axes[0].set(xlabel='Frozen development parent index',ylabel='Complex material MSE');axes[0].legend()
    comparison=risk['aggregation']['all_development']['comparisons']['gaussian_material_mean_minus_optimized_anchor/full_grid/complex_MSE']
    mean=comparison['difference_mean'];lo,hi=comparison['paired_bootstrap_95_interval']
    axes[1].errorbar([mean],[0],xerr=[[mean-lo],[hi-mean]],fmt='o',capsize=6,color='#b64c3a')
    axes[1].axvline(0,color='black',lw=1);axes[1].set_yticks([])
    axes[1].set(xlabel='Gaussian material mean minus anchor MSE (positive is worse)',title='Paired development bootstrap 95% interval')
    fig.suptitle('Specific GN reference worsens development risk; no FULL Bayes-risk comparison')
    save_figure(fig,'p3l_progress_gaussian_risk');names.append('p3l_progress_gaussian_risk')
    if geometry:
        fig,axes=plt.subplots(1,len(geometry['parents']),figsize=(12,4.5),squeeze=False,constrained_layout=True)
        for ax,parent in zip(axes.flat,geometry['parents']):
            for key,label in [('actual','Full-J geometry'),('frozen_Maxwell_changing_decoder','Frozen Maxwell / changing decoder')]:
                ax.plot([row['h'] for row in parent['rows']],[row[key]['normal_bracket_norm'] for row in parent['rows']],'o-',label=label)
            ax.set(xlabel='Finite-difference h',ylabel='Coordinate-normal bracket norm',title=parent['parent']);ax.legend(fontsize=8)
        fig.suptitle('Two-anchor projector-field diagnostic; no general integrability or posterior claim')
        save_figure(fig,'p3l_progress_geometry');names.append('p3l_progress_geometry')
    return names


def main():
    refs=references()
    ledger=read_json(OUT/'RESOURCE_LEDGER.json') # never open JSONL
    risk=read_json(OUT/'B/GAUSSIAN_IMAGE_RISK.json')
    departure=read_json(OUT/'B/GN_DEPARTURE_AUDIT.json')
    slice3=read_json(OUT/'B/slice3/RESULTS.json')
    geometry_path=OUT/'C/jobs/geometry_projector_DEV000_DEV028/RESULTS.json'
    geometry=read_json(geometry_path) if geometry_path.exists() else None
    names=figures(refs,risk,geometry)
    stamp=datetime.now(timezone.utc).isoformat()
    totals={key:sum(stage.get(key,0) for stage in ledger['stages'].values()) for key in ('RHS','states','LU','controlled_wall_seconds','physics_wall_seconds','GPU_seconds')}
    unqualified=sum(not row['stats'].get('reference_qualified',False) for row in refs)
    reference_table='<div class="scroll"><table><tr><th>机制 / 对象</th><th>材料 Rhat 最大</th><th>bulk ESS 最小</th><th>tail ESS 最小</th><th>材料均值 RMS MCSE</th><th>资格</th></tr>'
    for row in refs:
        reference_table+='<tr><td>'+link(row['path'],row['label'])+'</td>'+''.join('<td>'+fmt(row[k])+'</td>' for k in ('rhat','bulk','tail','mcse'))+'<td>'+html.escape(row['stats']['status'])+'</td></tr>'
    reference_table+='</table></div>'
    comparison=risk['aggregation']['all_development']['comparisons']['gaussian_material_mean_minus_optimized_anchor/full_grid/complex_MSE']
    last=slice3['rows'][-1]['refinement_from_previous_order']
    evidence=['START_HERE.md','RESEARCH_REPORT_ZH.md','THEORY_AND_ALGORITHM.md','POSTERIOR_REFERENCE_AUDIT.md','GO_NO_GO.md',
        'configs/PROTOCOL.json','RESOURCE_LEDGER.json','B/GN_DEPARTURE_AUDIT.json','B/GAUSSIAN_IMAGE_RISK.json','B/slice3/RESULTS.json',
        'C/jobs/geometry_projector_DEV000_DEV028/RESULTS.json','AUDIT/MEMORY_INTERRUPTION_RECOVERY.json','AUDIT/MEMORY_RECOVERY_REVIEW.md']
    nav='<nav>'+link('../../','B1 首页')+link('index.html','P3-L 进展')+link('report.html','中文报告')+link('START_HERE.md','匿名原始 Markdown')+'</nav>'
    top='<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>P3-L · 进行中的概率成像证据</title><style>'+existing_style()+'.running{border-left:5px solid #b57334}.negative{color:#923b2e}figcaption{font-size:14px;color:#516672}</style></head><body><main>'+nav
    body='<h1>P3-L：基础重建已建立，可信 FULL posterior 仍未建立</h1><p>RUNNING · 研究检查点 · '+html.escape(stamp)+'</p>'
    body+='<section class="running"><h2>当前边界</h2><p>本页只呈现已经落盘的证据。'+str(len(refs))+' 组 FULL 参考诊断中，'+str(unqualified)+' 组未合格；当前没有可用于 posterior distillation 的 teacher。后续作业与最终科学审查仍在进行，本页不自动生成最终 GO / NO-GO。</p><p>保持原 P1 的 64 Gaussian 场系数、2 uniform 属性与二维 scalar-TM/VIE。SLICE、PROFILE、Gaussian reference 和 FULL posterior 的含义分别保留。</p></section>'
    body+='<section><h2>费用来自资源快照</h2><div class="metrics"><span><b>'+f"{totals['RHS']:,}"+'</b>complex RHS</span><span><b>'+f"{totals['states']:,}"+'</b>paid state attempts</span><span><b>'+fmt(totals['controlled_wall_seconds']/3600)+'</b>受控 wall 小时</span><span><b>'+fmt(totals['GPU_seconds']/3600)+'</b>GPU 小时</span></div><p>快照 sequence '+str(ledger['sequence'])+'；活动 wall reservation '+html.escape(json.dumps(ledger.get('active_wall',{}),ensure_ascii=False))+'。拒绝、失败与内存中断费用不删除。快照不等同于进程存活检查。</p>'+link('RESOURCE_LEDGER.json','小型账本快照')+'</section>'
    body+='<section><h2>FULL 参考：接受率不能替代混合</h2><p>门槛为 Rhat&lt;1.01，材料函数 bulk/tail ESS≥1000，完整网格材料均值 RMS MCSE≤0.001；即使通过同机制门槛，仍需独立机制核验。UNASSESSED 常数分量保留，不能填作 ESS 无限大。</p>'+reference_table+'<img src="figures/p3l_progress_reference_diagnostics.png" alt="保存的参考诊断和原门槛"><img src="figures/p3l_progress_likelihood_traces.png" alt="不合格链的实际生产Phi轨迹"></section>'
    body+='<section><h2 class="negative">Gaussian 材料均值的实际负结果</h2><p>32 个 development parent 全部保留。GN Gaussian 材料均值相对优化 anchor 的 complex MSE 差为 +'+fmt(comparison['difference_mean'],7)+'，配对 bootstrap 95% 区间 '+html.escape(str(comparison['paired_bootstrap_95_interval']))+'。这是该具体近似的恶化证据，不能据此宣称 FULL Bayes 成像无益。</p><img src="figures/p3l_progress_gaussian_risk.png" alt="开发集材料风险与配对区间"><p>这里只读取此前独立离线评估的风险 JSON；渲染器不读取 GT。development 数据已暴露，不能称 locked 泛化。</p></section>'
    body+='<section><h2>强条件 GN reference 仍有明显非二次误差</h2><p>'+str(len(departure['rows']))+' 组中 '+str(sum(row['importance_ESS']<2 for row in departure['rows']))+' 组重要性 ESS&lt;2。小梯度不能修复整个 Gaussian 范围上的非二次偏差；decoder 和非 GN 曲率也有贡献，不能全部归因于 Maxwell feedback。</p>'+link('B/GN_DEPARTURE_AUDIT.json','保留线性梯度项后的审计')+'</section>'
    body+='<section><h2>SLICE3_BOX 有数值收敛证据，目标仍是条件切片</h2><p>最后两阶材料均值 RMS 差 '+fmt(last['material_mean_RMS_difference'])+'，logZ 差 '+fmt(last['logZ_difference'])+'。固定其他 61 个 Gaussian 坐标和两个属性，并条件化于冻结 box。它没有计算原 FULL 弱 marginal，也没有给出 box 外质量或经过认证的 quadrature 误差。</p>'+link('B/slice3/RESULTS.json','原始阶数与限制')+'</section>'
    if geometry:body+='<section><h2>几何结果是局部诊断</h2><img src="figures/p3l_progress_geometry.png" alt="两份固定anchor的有限差分投影场诊断"><p>冻结 Maxwell、只改变 decoder 已重现大部分 bracket。这里没有证明一般不可积性、exact fiber 或 posterior 概率质量。</p></section>'
    body+='<section><h2>原始说明与证据</h2><ul>'+''.join('<li>'+link(path)+'</li>' for path in evidence if (OUT/path).is_file())+'</ul><p>原始 chains、端点和 LU 留在原独立目录；没有在本渲染步骤制作新 Release、上传数组或执行推送。网页可达性与科学资格分别审查。</p></section>'
    (OUT/'index.html').write_text(top+body+'</main></body></html>')
    # Render current parent prose as read, without rewriting any Markdown.
    report='<p class="running">RUNNING · 此网页渲染已有研究记录，最终科学裁决尚未提交。生成时间 '+html.escape(stamp)+'</p>'
    for name in ('RESEARCH_REPORT_ZH.md','THEORY_AND_ALGORITHM.md','POSTERIOR_REFERENCE_AUDIT.md','GO_NO_GO.md'):
        path=OUT/name
        if path.exists():report+='<article>'+render_markdown(path.read_text())+'</article>'
    (OUT/'report.html').write_text(top+report+'</main></body></html>')
    receipt=dict(status='LOCAL_RUNNING_PROGRESS_RENDER_NOT_PUBLICATION_CONFIRMATION',rendered_at_utc=stamp,
        repository='migodam/b1-probabilistic-maxwell-imaging',existing_path=REL,new_repository=False,
        reference_groups=len(refs),qualified_full_references=sum(row['stats'].get('reference_qualified',False) for row in refs),
        ledger_sequence=ledger['sequence'],ledger_snapshot_only=True,ledger_journal_read=False,
        figures=names,physics_RHS=0,GT_arrays_read=False,training_runs=0,
        markdown_unmodified=True,push_executed=False,final_scientific_verdict='PENDING_PARENT_REVIEW')
    (OUT/'PROGRESS_RENDER.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(receipt,ensure_ascii=False))


if __name__=='__main__':main()
