"""Anonymous publication readback; no physics and no hashes."""
from pathlib import Path
import concurrent.futures,datetime,json,urllib.request,urllib.error
ROOT=Path(__file__).resolve().parents[1]
COMMIT='012cdfa95ae39449a867ae3dbf6864a974d69ead'
PREFIX='https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_w_gq_r3_joint/'
RAW='https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/'+COMMIT+'/research/research_w_gq_r3_joint/'
RELEASE='https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-w-gq-r3-ji-20261030'
ASSET=RELEASE.replace('/tag/','/download/')+'/RESEARCH_W_GQ_R3_OFFLINE_AUDIT.zip'
targets=[]
for name in ['index.html','report.html','formulas.html','gallery.html','digits.html','handoff.html','run.html','evidence.html']:
    targets.append((PREFIX+name,'GET',name))
for name in ['REPORT_ZH.md','FORMULA_COMPENDIUM_ZH.md','GO_NO_GO.json']:
    targets.append((RAW+name,'GET',name))
for name in ['figures/R0_TARGET_MATERIAL_OVERVIEW.png','figures/R0_FULLGN_CONTOUR_ZOOM.png','exploration/digits/truth_vs_reconstruction.png','research/delegated/sources/CHART_MECHANISM_RUN/FIGURES/public_2x2_mechanism.png']:
    targets.append((PREFIX+name,'GET',name))
targets.extend([(RELEASE,'GET',None),(ASSET,'HEAD',None)])
def check(item):
    url,method,local=item
    out={'url':url,'method':method}
    try:
        req=urllib.request.Request(url,method=method,headers={'User-Agent':'Research-W-R3-public-readback'})
        with urllib.request.urlopen(req,timeout=45) as r:
            out.update(status=r.status,content_length=r.headers.get('Content-Length'),content_type=r.headers.get('Content-Type'))
            if method=='GET':
                b=r.read();out['bytes']=len(b)
                if local:out['equals_local_bytes']=b==(ROOT/local).read_bytes()
            else:out['expected_asset_bytes_match']=int(out['content_length'] or -1)==193037632
    except Exception as e:out['error']=str(e)
    return out
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex: results=list(ex.map(check,targets))
out={'checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_commit':COMMIT,'authenticated':False,'new_physics':0,'checks':results,
     'all_pass':all(x.get('status')==200 and x.get('equals_local_bytes',True) and x.get('expected_asset_bytes_match',True) for x in results)}
with (ROOT/'publication/READBACK_ATTEMPTS.jsonl').open('a') as f:f.write(json.dumps(out,ensure_ascii=False)+'\n')
(ROOT/'publication/PUBLICATION_READBACK.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
print(json.dumps({'all_pass':out['all_pass'],'checks':len(results),'failures':[x for x in results if x.get('status')!=200 or not x.get('equals_local_bytes',True) or not x.get('expected_asset_bytes_match',True)]},ensure_ascii=False))

