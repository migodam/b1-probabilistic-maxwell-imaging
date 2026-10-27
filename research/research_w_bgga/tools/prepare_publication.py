"""Reconcile readiness metadata and add one entry to the existing B1 site."""
from pathlib import Path
import ast,json,re,subprocess
from finalize_results import main as reconcile
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
    reconcile()
    index=ROOT.parent/'index.html';text=index.read_text()
    start,end='<!-- BGGA_START -->','<!-- BGGA_END -->'
    card=start+'<section><h2><a href="research_w_bgga/">Research W-BGGA</a></h2><p>3934 RHS / 482 FE真实全波：公共空间精确吸收通过；完整端口与默认区间没有净节约，当前加速路线NO-GO。有限joint概率与单点化风险独立报告，G2未启动。</p></section>'+end
    if start in text:text=re.sub(re.escape(start)+r'.*?'+re.escape(end),lambda m:card,text,flags=re.S)
    else:text=text.replace('</main>',card+'</main>',1)
    index.write_text(text)
    readme=REPO/'README.md';text=readme.read_text()
    entry='\n<!-- BGGA_START -->\n## Research W-BGGA\n\n[真实全波与完整联合概率实验报告](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_w_bgga/)：3934 RHS / 482 FE；精确公共空间吸收通过，当前完整端口/默认区间加速路线NO-GO。G1为旧bank有限prior重放，G2未启动。[原始数据](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-w-bgga-20261028)。\n<!-- BGGA_END -->\n'
    if start in text:text=re.sub(re.escape(start)+r'.*?'+re.escape(end),lambda m:entry.strip(),text,flags=re.S)
    else:text+=entry
    readme.write_text(text)
    (ROOT/'PUBLICATION_SOURCE_BASE.json').write_text(json.dumps(dict(
        publication_parent=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
        science_sources=json.loads((ROOT/'CONTRACT.json').read_text())['source_commits'],
        note='Publication follows current main without replacing pinned scientific source/results; old directories untouched.'),indent=2)+'\n')
    print('Own result-site entry, source syntax, links and actual cost metadata prepared. No PDE.')
if __name__=='__main__':main()
