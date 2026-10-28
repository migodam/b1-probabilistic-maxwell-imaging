# -*- coding: utf-8 -*-
"""Presentation metadata and an adapted stdlib scientific reading renderer."""
from pathlib import Path
import json
P=Path(__file__).resolve().parents[1]
REPO=P.parents[1]
RELEASE='https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p3lr7-material-fibers-20261028'

def main():
    (P/'EVIDENCE_PACKAGE.md').write_text('''# R7阅读与无损证据

[现有B1的R7 Release]('''+RELEASE+''')提供完整原始证据、manifest、独立校验收据及GPT Pro阅读包。

- `Research_P3LR7_Results_and_GPT_Pro_Handoff.zip`：中文报告、源码、图、compact八对象数值与必要复核入口；不是全部过程数据。
- `R7_Complete_Evidence.tar.zst`：所有新增科学数组、封存预测、逐trial current/field/J、账本、失败与必要原始复用输入。恢复后才能按原相对路径np.load。
- `R7_Evidence_Manifest.json`：逐文件尺寸、BLAKE2b、时间与来源。
- `R7_Delivery_Verification.json`：无损逐字节恢复核验；最终匿名读回另有收据。

范围只包含本轮R7及明确必要的约10.2MB历史复用输入，不重复上传R6九GB完整包、历史训练数据或虚拟环境。43个必要复用文件复制后逐字节相等，ASSET_MANIFEST记录原始路径与恢复层次；没有保存LU矩阵，保留current/field/J/T和可重建状态。

原绝对路径是来源证据；在另一台机器恢复到新的目录层次即可，不覆盖不同字节的已有文件。纯数组审计支持新输出目录与路径重定位。归档保留原文件逻辑字节与科学字段，打包、发布和读回0新增Maxwell RHS，不提高任何科学资格。没有SHA-256附加检查，没有删除原始科学数据。
''')
    (P/'RESOURCE_ACCOUNTING.md').write_text('''# R7实际费用

总新增37,456 complex RHS /50,000上限；未用12,544。forward10,600、derivative26,856；state/LU各1,325；19个物理job、物理动作失败0。拒绝试探、恢复及修正重放全部计费。

|用途|RHS|
|---|---:|
|首次实际probe|200|
|最终DEV4|7200|
|保留的DEV修正前/重放|2952|
|八份新truth|64|
|八对象anchor、探索、路径|26336|
|普通DEV20态profile|576|
|同材料与交叉copula物理控制|128|
|合计|37456|

物理动作wall14.7676秒；包含持久化/setup的物理job wall72.2100秒、CPU67.4749秒；物理进程峰值189,972,480 bytes。E1另7.18秒与490MB峰值。编辑、文献、绘图、打包和交付不伪称都被物理job钟覆盖。没有训练、GPU或云服务。

本轮fresh anchor是direct-Probit64，目标是几何；不能与R6不同优化流程的1.73秒换算成成像系统加速比。R6的6,860,984 RHS不重计；九次合法T/J复用的cold-equivalent288 derivative RHS是历史cache权利，不是新支付。

中央追加账本在完整Release；snapshot及独立逐行replay在仓库。C的24 RHS setup按四档整体计一次，缺失端点不得构成加速主张。
''')
    source=(REPO/'research/p3l_r6_task_optimal_maxwell/build_site.py').read_text().replace('R6','R7')
    a=source.index('IDENTITIES={');b=source.index('\nFIXED_APPENDICES',a)
    source=source[:a]+'''IDENTITIES={
 'decoder_quotient':('参数gauge','不同参数可表示同一材料；分层秩和有限同材料见证分开报告。'),
 'material_geometry':('真实材料弱响应','完整Maxwell与实际材料RMS；有限端点不自动是连续exact fiber。'),
 'task_scope':('任务与概率资格','局部有限见证不等于posterior质量、可信区间或Bayesian风险。')}
'''+source[b:]
    source=source.replace('先分清物理问题与三类输出','先分清参数、材料、物理与概率')
    source=source.replace('R7任务最优Maxwell材料计算：材料风险、低维条件概率、真实物理成本分别呈现','R7材料商空间与真实Maxwell有限响应：gauge、弱方向、路径及任务资格')
    mathconfig={'tex':{'inlineMath':[['\\(','\\)']],'displayMath':[['\\[','\\]'],['$$','$$']]},'svg':{'fontCache':'global'}}
    script='<script>window.MathJax='+json.dumps(mathconfig,separators=(',',':'))+';</script><script defer src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-svg.js"></script>'
    source=source.replace('class Builder:','MATH_SCRIPT='+repr(script)+'\n\nclass Builder:')
    source=source.replace("+STYLE+'</style></head><body>","+STYLE+'</style>'+MATH_SCRIPT+'</head><body>")
    source=source.replace('JS_bytes=0,remote_assets_added=0','JS_bytes=len(MATH_SCRIPT.encode()),remote_assets_added=1')
    # Preserve display math as a single text block. Rendering has a readable fallback.
    old="""   elif stripped=='$$':
    flush();mathlines=[];i+=1
    while i<len(lines) and lines[i].strip()!='$$':mathlines.append(lines[i]);i+=1
    out.append('<pre class="math" aria-label="数学公式">'+html.escape('\\n'.join(mathlines))+'</pre>')"""
    new="""   elif stripped in ('$$',r'\\['):
    closer='$$' if stripped=='$$' else r'\\]'
    flush();mathlines=[];i+=1
    while i<len(lines) and lines[i].strip()!=closer:mathlines.append(lines[i]);i+=1
    out.append('<div class="math" aria-label="数学公式">'+html.escape(r'\\['+'\\n'+'\\n'.join(mathlines)+'\\n'+r'\\]')+'</div>')"""
    if old in source:source=source.replace(old,new)
    else:
        # Record rather than silently claiming an unimplemented renderer change.
        raise RuntimeError('display-math renderer source pattern changed')
    (P/'build_site.py').write_text(source)
    summary=json.loads((P/'RESULTS_SUMMARY.json').read_text())
    appendices=[{'path':str(f.relative_to(P))} for f in sorted(P.glob('*.md')) if f.name not in ['START_HERE.md','RESEARCH_REPORT_ZH.md']]
    appendices += [{'path':f} for f in ['audit/results/REVIEW.md','audit/final_review/REVIEW.md','audit/literature/FOCUSED_PRIOR_ART_AUDIT_ZH.md','audit/E1_LOCAL_IFT_QUALIFICATION.md']]
    summary['site']=dict(title='P3L-R7｜Maxwell材料解纤维与任务可辨识性',subtitle='E0/E1零物理 · 四开发＋八新对象 · 37,456 RHS封账 · 高阶独立优势未建立',
      intro_markdown='不同潜变量不一定是不同材料。不同材料在有限噪声下难以区分，也不自动构成精确解纤维，更不自动带有posterior质量。R7用真实材料RMS和完整Maxwell将这些对象分开验证。',
      highlights=[{'label':'新增complex RHS','value':'37,456 / 50,000','note':'包含拒绝和修正重放','identity':'material_geometry'},
       {'label':'.02 RMS Full-J竞争材料','value':'8 / 8','note':'另7/8满足pair≤1；不是credible mass','identity':'material_geometry'},
       {'label':'高阶独立增量','value':'NO-GO','note':'C6/8，B8/8；Full66 posterior未认证','identity':'task_scope'}],
      appendices=appendices,section_sources={'cost':'RESOURCE_ACCOUNTING.md'},downloads=[{'path':RELEASE,'label':'完整原始证据与GPT Pro包'}])
    (P/'RESULTS_SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
    print('R7 presentation prepared; zero physical calls')

if __name__=='__main__':main()
