"""Anonymous HTTP/raw/asset readback. No physics, credentials or SHA checks."""
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.parse import urljoin,urlparse
import argparse,datetime,json,time,resource
from zoneinfo import ZoneInfo

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--commit',required=True);a=ap.parse_args()
    root=Path(__file__).resolve().parents[1];note=json.loads((root/'publication/PUBLICATION_NOTE.json').read_text());start=time.perf_counter()
    base='https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_w_gsda_r4/'
    rawbase='https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/'+a.commit+'/research/research_w_gsda_r4/'
    headers={'User-Agent':'R4-anonymous-publication-readback','Cache-Control':'no-cache'}
    rows=[];errors=[]
    def small(url,local=None):
        try:
            with urlopen(Request(url,headers=headers),timeout=60) as f:data=f.read();status=f.status;ctype=f.headers.get('Content-Type')
            row=dict(url=url,status=status,bytes=len(data),content_type=ctype,authentication=False)
            if local is not None:row['bytes_equal_to_local']=data==local.read_bytes()
            if status!=200 or (local is not None and not row['bytes_equal_to_local']):errors.append(row)
            rows.append(row);return data
        except Exception as e:row=dict(url=url,error=str(e),authentication=False);rows.append(row);errors.append(row);return None
    for name in ['index.html','report.html','gallery.html','methods.html','physics.html','handoff.html','evidence.html','failures.html']:
        small(base+name,root/name)
    for name in ['REPORT_ZH.md','WEB_CHAT_HANDOFF.md','GO_NO_GO.json','DIGITS_1_TO_9_METRICS.csv','EXPERIMENT_RESULTS.csv']:
        small(rawbase+name,root/name)
    for name in ['DIGITS_1_TO_9_CONTRAST.png','DIGITS_1_TO_9_LOSS.png','figures/tv_ablation/DIGIT8_TV0_VS_TV02.png','figures/paired_controls_with_B0/NATIVE_CONTROLS_CONTRAST.png','figures/physics/GD_MODES_AND_FEEDBACK.png']:
        small(base+name,root/name)
    api='https://api.github.com/repos/migodam/b1-probabilistic-maxwell-imaging/releases/tags/research-w-gsda-r4-20261030'
    payload=small(api);assets=[]
    if payload:
        release=json.loads(payload);assets=release.get('assets',[])
        if release.get('draft'):errors.append(dict(error='release is still draft'))
    assetrows=[]
    for filename,expected in [('RESEARCH_W_GSDA_R4_20261030.zip',note['full_bytes']),('RESEARCH_W_GSDA_R4_WEB_CHAT.zip',note['web_bytes'])]:
        found=[x for x in assets if x['name']==filename]
        if len(found)!=1 or found[0]['size']!=expected or found[0]['state']!='uploaded':
            errors.append(dict(asset=filename,error='anonymous API asset size/state mismatch'));continue
        asset=found[0];url=asset['browser_download_url'];local=a.source/'delivery'/filename
        record=dict(name=filename,url=url,expected_bytes=expected,api_size=asset['size'],state=asset['state'],authentication=False)
        try:
            with urlopen(Request(url,headers=headers,method='HEAD'),timeout=60) as f:
                record['HEAD_status']=f.status;record['HEAD_bytes']=int(f.headers.get('Content-Length',-1));record['download_host']=urlparse(f.url).netloc
            if record['HEAD_status']!=200 or record['HEAD_bytes']!=expected:raise RuntimeError('anonymous HEAD/size mismatch')
            if 'WEB_CHAT' in filename:
                count=0;equal=True
                with urlopen(Request(url,headers=headers),timeout=60) as remote,local.open('rb') as own:
                    record['GET_status']=remote.status
                    while True:
                        block=remote.read(1024*1024)
                        if not block:break
                        count+=len(block)
                        if block!=own.read(len(block)):equal=False;break
                    if own.read(1):equal=False
                record.update(verification='complete anonymous streamed byte comparison, no second local download or hash',downloaded_bytes=count,all_bytes_equal=equal and count==expected)
                if not record['all_bytes_equal']:raise RuntimeError('web archive byte readback mismatch')
            else:
                record['verification']='anonymous HEAD size plus first/last 1024-byte range comparisons; not full-byte readback'
                samples=[]
                for lo,hi in [(0,1023),(expected-1024,expected-1)]:
                    h={**headers,'Range':f'bytes={lo}-{hi}'}
                    with urlopen(Request(url,headers=h),timeout=60) as f:
                        block=f.read(1024);status=f.status;content_range=f.headers.get('Content-Range')
                    if status!=206:raise RuntimeError('range server did not honor partial request')
                    with local.open('rb') as f:f.seek(lo);equal=block==f.read(hi-lo+1)
                    samples.append(dict(start=lo,end=hi,status=status,content_range=content_range,bytes_equal=equal))
                    if not equal:raise RuntimeError('full archive range mismatch')
                record['range_samples']=samples
        except Exception as e:record['error']=str(e);errors.append(dict(asset=filename,error=str(e)))
        assetrows.append(record)
    out=dict(status='PASS' if not errors else 'READBACK_INCOMPLETE',checked_at_singapore=datetime.datetime.now(ZoneInfo('Asia/Singapore')).isoformat(),content_commit=a.commit,public_page=base,raw_report=rawbase+'REPORT_ZH.md',release=note['release'],http_rows=rows,assets=assetrows,errors=errors,new_FE=0,new_RHS=0,verification_process_CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-start,publication_does_not_change_scientific_gates=True,sha256_not_used=True)
    (root/'publication/PUBLICATION_READBACK.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(status=out['status'],http_rows=len(rows),assets=assetrows,errors=errors),ensure_ascii=False),flush=True)
    if errors:raise SystemExit(1)

if __name__=='__main__':main()
