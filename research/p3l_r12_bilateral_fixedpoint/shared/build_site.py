"""Render only R12 documents; no experiments, formula reset or numerical replay."""
from pathlib import Path
import importlib.util,html,os,shutil
P=Path(__file__).resolve().parents[1];REPO=P.parents[1]
spec=importlib.util.spec_from_file_location('r12_reader_base',REPO/'research/p3l_r9_stabilized_subspace/src/build_reading_site.py');base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
class Reader(base.Reader):
    def shell(self,title,body,page,subtitle='精确双边模式 · 普通一阶延拓 · 实际材料覆盖 · 数字拓扑'):
        links=[(P/'index.html','R12完整报告'),(P/'DIGITS_VISUAL_REVIEW_ZH.html','数字视觉结果'),(P/'GO_NO_GO.html','科研裁决'),(P/'RESEARCH_REPORT_ZH.md','原始Markdown'),(REPO/'research/maxwell_formula_atlas/index.html','167条公式')]
        nav='<nav>'+''.join('<a href="'+os.path.relpath(a,page.parent)+'">'+t+'</a>' for a,t in links)+'<a href="https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p3lr12-bilateral-fixedpoint-20261030">完整证据Release</a></nav>'
        note='<div class="notice"><strong>开发结论：分层扩大 NO-GO。</strong>四个暴露P1难例主覆盖0/4；数字8未恢复双孔。163,994 RHS / 6,974 LU，在本轮上限内。没有独立新parent或Full66 posterior声明。</div>'
        return '<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+html.escape(title)+'</title><style>'+base.old.STYLE+' .math,mjx-container{max-width:100%;overflow-x:auto;overflow-y:hidden}img{max-width:100%;height:auto}table{display:block;overflow-x:auto}</style>'+base.MATH+'</head><body><header>'+nav+'<h1>'+html.escape(title)+'</h1><p>'+subtitle+'</p></header><main>'+note+body+'</main><footer>完整非线性Maxwell；原P1与N原生材料合同分开。候选频率不是posterior mass。</footer></body></html>'
def main():
    docs=list(P.glob('*.md'))+list((P/'GPT_PRO_HANDOFF').glob('*.md'))
    for d in ('agent_a','agent_b','agent_c','agent_a_max','agent_c_max'):docs+=list((P/d).glob('*.md'))
    Reader(P).render(docs)
    shutil.copy2(P/'RESEARCH_REPORT_ZH.html',P/'index.html');shutil.copy2(P/'RESEARCH_REPORT_ZH.html',P/'report.html')
    print('R12_DOCUMENTS_RENDERED',len(docs))
if __name__=='__main__':main()
