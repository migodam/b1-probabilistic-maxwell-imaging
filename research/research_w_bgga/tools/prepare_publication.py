"""Reconcile readiness metadata and add one entry to the existing B1 site."""
from pathlib import Path
import ast,json,re,subprocess
from bgga.accounting import BudgetBook
ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
def main():
    for folder in ['src','tools','tests']:
        for p in (ROOT/folder).rglob('*.py'):ast.parse(p.read_text())
    bad=[]
    for p in ROOT.glob('*.md'):
        for link in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',p.read_text()):
            if link.startswith(('http','#')):continue
            if not (p.parent/link.split('#')[0]).exists():bad.append((p.name,link))
    if bad:raise ValueError(str(bad))
    book=BudgetBook();events=book.records()
    c=json.loads((ROOT/'COST_SUMMARY.json').read_text())
    c.update(measured_wrapper_process_cpu_seconds=sum(r.get('process_cpu_seconds',0) for r in events),
        wrapper_wall_seconds=sum(r.get('wall_seconds',0) for r in events if r.get('label')=='job_completion'),
        peak_measured_job_RSS_bytes=max([r.get('children_peak_RSS',0) for r in events]+[0]),
        snapshot_label='completed pre-publication jobs; final receipt for this metadata job is in local ledger',
        accounting_gap='Early source/PDF extraction, worker static reads and publication/network not fully timed. Not zero; no new physical work.',
        pending_current_report_job='prepare_publication current receipt not yet complete; no end-to-end speed claim')
    (ROOT/'COST_SUMMARY.json').write_text(json.dumps(c,ensure_ascii=False,indent=2)+'\n')
    index=ROOT.parent/'index.html';text=index.read_text()
    start,end='<!-- BGGA_START -->','<!-- BGGA_END -->'
    card=start+'<section><h2><a href="research_w_bgga/">Research W-BGGA</a></h2><p>固定公共网格与区间概率代码、32项测试完成；新增全波0，G0额度待授权。公开实际端口成本，不宣称加速。</p></section>'+end
    if start in text:text=re.sub(re.escape(start)+r'.*?'+re.escape(end),lambda m:card,text,flags=re.S)
    else:text=text.replace('</main>',card+'</main>',1)
    index.write_text(text)
    readme=REPO/'README.md';text=readme.read_text()
    entry='\n<!-- BGGA_START -->\n## Research W-BGGA\n\n[静态Green吸收实现与资格前检查](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_w_bgga/)：固定公共目标网格和概率区间核，32项测试通过；新物理待独立资源授权，不写成G0或加速PASS。\n<!-- BGGA_END -->\n'
    if start in text:text=re.sub(re.escape(start)+r'.*?'+re.escape(end),lambda m:entry.strip(),text,flags=re.S)
    else:text+=entry
    readme.write_text(text)
    (ROOT/'PUBLICATION_SOURCE_BASE.json').write_text(json.dumps(dict(
        publication_parent=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
        science_sources=json.loads((ROOT/'CONTRACT.json').read_text())['source_commits'],
        note='Publication follows current main without replacing pinned scientific source/results; old directories untouched.'),indent=2)+'\n')
    print('Own static site entry and pre-publication costs prepared; local report links valid. No PDE.')
if __name__=='__main__':main()
