"""Prepare a compact reading kit and identities; no experiments or historical edits."""
from pathlib import Path
import csv,datetime,hashlib,json,re,subprocess,zipfile
P=Path(__file__).resolve().parents[1];R=P.parents[1];AT=R/'research/maxwell_formula_atlas';OUT=P/'packaging'
tag='p3lr10-bilateral-20261030';url='https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/'+tag+'/'
root_commit=subprocess.check_output(['git','rev-parse','b900f1c0'],cwd=R,text=True).strip()
# Append precise formula-source locations; old formula entries remain unchanged.
idx=AT/'SOURCE_FORMULA_INDEX.csv'
with idx.open() as f: rr=csv.DictReader(f);fields=rr.fieldnames;rows=list(rr)
rows=[r for r in rows if not r['formula_id'].startswith('F-BI-')]
chapter=AT/'09_BILATERAL_SHAPE_AMPLITUDE.md';lines=chapter.read_text().splitlines();heads=[]
for i,l in enumerate(lines,1):
 m=re.match(r'## (F-BI-\d+) · (.+)',l)
 if m:heads.append((i,m[1],m[2]))
for j,(line,ident,title) in enumerate(heads):
 rows.append(dict(formula_id=ident,title=title,source_path=str(chapter.relative_to(R)),source_function_or_section=ident,line_start=line,line_end=heads[j+1][0]-1 if j+1<len(heads) else len(lines),historical_identity='R10 scientific snapshot '+root_commit,status='DERIVATION_WITH_EXPLICIT_CONDITIONS',coordinate_and_claim_caveat='Standard identities, conditional bounds, empirical costs and unrun mechanisms separated'))
with idx.open('w') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

card='''<!-- P3LR10_CURRENT --><section class="card"><h2>P3L-R10-BI · 双边传递与形状—幅度物理复用</h2><p>固定形状320查询全部通过数值误差门；32查询冷形状费用相对直接LU节省51.82%。新10对象中双边提案覆盖2/10，Full66随机目标3/10，未出现材料覆盖优势。完整失败、公式与论文草稿分别交付。</p><p><a href="{prefix}p3l_r10_bilateral_structure/">中文报告与实际材料画廊 →</a> · <a href="{prefix}maxwell_formula_atlas/">140条持续推导公式 →</a></p></section><!-- /P3LR10_CURRENT -->'''
for path,prefix in [(R/'index.html','research/'),(R/'research/index.html','')]:
 s=path.read_text();s=re.sub(r'<!-- P3LR10_CURRENT -->.*?<!-- /P3LR10_CURRENT -->','',s,flags=re.S)
 s=s.replace('<main>','<main>'+card.format(prefix=prefix),1);path.write_text(s)
block='''<!-- P3LR10_CURRENT -->
## P3L-R10-BI：双边传递、条件物理复用与材料覆盖边界

[阅读完整中文报告](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/p3l_r10_bilateral_structure/) · [140条公式集](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/maxwell_formula_atlas/) · [无损证据Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p3lr10-bilateral-20261030)

固定形状的320次新幅度查询全部通过误差门；实际32节点冷形状费用节省51.82%。新材料对象的覆盖优势未成立，Full66 posterior仍未认证，强TAP贡献HOLD。完整源码、失败、独立QA与论文草稿已保存。
<!-- /P3LR10_CURRENT -->

'''
for path in [R/'README.md',R/'START_HERE.md',R/'research/README.md']:
 s=path.read_text();s=re.sub(r'<!-- P3LR10_CURRENT -->.*?<!-- /P3LR10_CURRENT -->\s*','',s,flags=re.S);path.write_text(block+s)
f=P/'PAPER_NARRATIVE_V2_EN.md';f.write_text(f.read_text().replace('cold kernel time is 0.1150 seconds','aggregate cold kernel time across all ten shapes is 0.1150 seconds'))

def bh(p):
 h=hashlib.blake2b(digest_size=32)
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
items=[]
for p in sorted(OUT.glob('*.tar.zst')):items.append(dict(name=p.name,bytes=p.stat().st_size,blake2b_256=bh(p),url=url+p.name))
manifest=dict(release_tag=tag,scientific_commit=root_commit,created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),assets=items,restore='tar --zstd -xf ASSET; R10 new evidence uses repository-relative paths. Dependency archive uses B1 sibling layout; read NOTES first. Never overwrite historical assets.',truth_scope='synthetic GT published only after sealed independent scoring',history_read_only=True)
(P/'RELEASE_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
(OUT/'RELEASE_MANIFEST.json').write_bytes((P/'RELEASE_MANIFEST.json').read_bytes())
print(json.dumps(dict(archives=len(items),total_archive_bytes=sum(x['bytes'] for x in items),scientific_commit=root_commit)))
