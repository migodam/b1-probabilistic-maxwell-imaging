"""Read selected files from an authorized public split GitHub release ZIP.
Uses normal HTTP Range reads; no authentication/captcha bypass.
"""
from pathlib import Path
import os
import io, zipfile, urllib.request, struct, zlib, binascii, json, time
from concurrent.futures import ThreadPoolExecutor

ROOT=Path(os.environ.get('MAX_STRUCT_RUN_ROOT', str(Path(__file__).resolve().parent.parent)))/'work'
DEST=ROOT/'paired_selected_candidates'
DEST.mkdir(exist_ok=True)
SIZES=[134217728,42366822]
PREFIX='https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/research-max-struct-r12-paired-20261112/paired_results_full.zip.part'
TOTAL=sum(SIZES)
CACHE=[]
tail=ROOT/'paired_release_tail.bin'
if tail.exists() and tail.stat().st_size==65536: CACHE.append((TOTAL-65536,tail.read_bytes()))
TRANSFER=[]

def read_at(offset,n):
    if n<=0:return b''
    for start,b in CACHE:
        if start<=offset and offset+n<=start+len(b):return b[offset-start:offset-start+n]
    parts=[];current=offset;stop=min(TOTAL,offset+n)
    while current<stop:
        i=0 if current<SIZES[0] else 1
        base=0 if i==0 else SIZES[0]
        a=current-base;b=min(stop,base+SIZES[i])-base-1
        req=urllib.request.Request(PREFIX+f'{i+1:03d}',headers={'Range':f'bytes={a}-{b}'})
        with urllib.request.urlopen(req,timeout=50) as r:
            if r.status!=206:raise RuntimeError(f'Expected206, received{r.status}; stop rather than silently read full file')
            data=r.read()
            if len(data)!=b-a+1:raise RuntimeError('Unexpected Range length')
        parts.append(data);TRANSFER.append(dict(part=i+1,start=a,end=b,bytes=len(data)))
        current+=len(data)
    return b''.join(parts)

class Remote(io.RawIOBase):
    def __init__(self):self.pos=0
    def readable(self):return True
    def seekable(self):return True
    def tell(self):return self.pos
    def seek(self,p,whence=0):
        self.pos=p if whence==0 else self.pos+p if whence==1 else TOTAL+p
        return self.pos
    def read(self,n=-1):
        n=TOTAL-self.pos if n<0 else n
        b=read_at(self.pos,n);self.pos+=len(b);return b

WANTED=[f'{cohort}__{case}__{method}' for cohort,case in [
    ('A','r10_twod_0004'),('B','r12_confirm24_0001'),('B','r12_confirm24_0003'),
    ('B','r12_confirm24_0022'),('C','r8_test_0001')]
    for method in ['ap','dr']]

def extract(info):
    name=info.filename
    # Local extra fields are small, but their actual size is read from the header.
    header=read_at(info.header_offset,30)
    v=struct.unpack('<4s5H3L2H',header)
    if v[0]!=b'PK\x03\x04':raise RuntimeError('Bad local signature')
    name_len,extra_len=v[-2:]
    compressed=read_at(info.header_offset+30+name_len+extra_len,info.compress_size)
    if info.compress_type==zipfile.ZIP_STORED:data=compressed
    elif info.compress_type==zipfile.ZIP_DEFLATED:data=zlib.decompress(compressed,-15)
    else:raise RuntimeError(f'Unsupported compression{info.compress_type}')
    if len(data)!=info.file_size or binascii.crc32(data)&0xffffffff != info.CRC:raise RuntimeError('CRC/length mismatch')
    task=next(t for t in WANTED if '/'+t+'/' in '/'+name)
    dest=DEST/(task+'__CANDIDATES.npz');dest.write_bytes(data)
    return dict(task_id=task,archive_name=name,local_path=str(dest),size=len(data),compressed_size=info.compress_size,CRC=info.CRC)

if __name__=='__main__':
    start=time.perf_counter()
    with zipfile.ZipFile(Remote()) as z:
        info=[i for i in z.infolist() if i.filename.endswith('/CANDIDATES.npz') and any('/'+t+'/' in '/'+i.filename for t in WANTED)]
        print('archive_files',len(z.infolist()),'selected',len(info),flush=True)
        (DEST/'archive_inventory.json').write_text(json.dumps([dict(name=i.filename,size=i.file_size,compressed=i.compress_size) for i in z.infolist()],indent=2))
    if len(info)!=len(WANTED):raise RuntimeError('Not all ten requested candidate files found')
    with ThreadPoolExecutor(max_workers=4) as ex: results=list(ex.map(extract,info))
    receipt=dict(source=PREFIX,archive_bytes=TOTAL,selected_files=results,transferred_bytes=sum(t['bytes'] for t in TRANSFER)+65536,
                 transfers=TRANSFER,wall_seconds=time.perf_counter()-start)
    (DEST/'retrieval_receipt.json').write_text(json.dumps(receipt,indent=2))
    print(json.dumps({k:receipt[k] for k in ['transferred_bytes','wall_seconds']}),flush=True)
