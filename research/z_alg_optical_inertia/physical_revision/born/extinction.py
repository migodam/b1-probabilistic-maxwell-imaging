"""A fixed receive weighting tests first Born without estimating material."""
import os
for k in ('OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'):
    os.environ[k]='1'
from pathlib import Path
import csv,json
from time import perf_counter,process_time
import numpy as np
from scipy.linalg import svd
from flint import arb
from run import mat
HERE=Path(__file__).resolve().parent
STUDY=HERE.parents[1]
def weights(s,b,cut=.01):
    v,sigma,vh=svd(s,full_matrices=False,check_finite=False)
    keep=sigma>=cut*sigma[0]
    c=(v[:,keep]/sigma[keep])@(vh[keep]@b)
    q=np.sum((s.conj().T@c)*b.conj(),axis=1)
    if np.min(q.real)<=0:
        return dict(rank=int(keep.sum()),min_re=float(np.min(q.real)),status='one phase cannot be fixed with positive real part'),None
    eta=max(0.,float(np.max(-q.imag/q.real)))+2.e-8
    z=(1j-eta)*c
    return dict(rank=int(keep.sum()),min_re=float(np.min(q.real)),eta=eta,
            relative_field_defect=float(np.linalg.norm(s.conj().T@c-b)/np.linalg.norm(b))),z
def geometry_check(s,b,z):
    zw=mat(z);norm=sum((abs(zw[i,j])**2 for i in range(zw.nrows()) for j in range(zw.ncols())),arb(0)).sqrt()
    peak=None
    for i in range(0,len(b),384):
        zz=mat(s[:,i:i+384].conj().T)*zw
        bb=mat(b[i:i+384].conj())
        for j in range(bb.nrows()):
            v=sum(((zz[j,l]*bb[j,l]).real for l in range(bb.ncols())),arb(0)).upper()
            peak=v if peak is None else peak.max(v)
    return zw,norm,peak
def readout(zw,norm,y):
    v=zw.conjugate().transpose()*mat(y)
    p=sum((v[j,j].real for j in range(v.nrows())),arb(0))
    return max(0.,np.nextafter(float((p/norm).lower()),-np.inf))
def main():
    root=STUDY/'physical_revision/array/observations'
    manifest=json.loads((root/'manifest.json').read_text())['observations']
    arr=np.load(root/'arrays.npz'); results=[];geos=[];wv={}
    for ns in (10,16,24,32):
      for nr in (20,32):
        t=perf_counter();cpu=process_time()
        s,b=arr[f'S{nr}'],arr[f'B{ns}'];row,z=weights(s,b)
        row.update(s=ns,p=nr)
        if z is None:
            geos.append(row);continue
        zw,norm,g=geometry_check(s,b,z)
        row.update(norm_up=float(norm.upper()),largest_born_cell_upper=np.nextafter(float(g.upper()),np.inf),valid=bool(g<0),wall=perf_counter()-t,cpu=process_time()-cpu)
        geos.append(row)
        if g>=0:continue
        wv[f'w_{ns}_{nr}']=z
        for obs in manifest:
            if obs['sources']!=ns or obs['receivers']!=nr:continue
            y=np.load(root/obs['file'])['Y']; radius=obs['padded_radius'];t=perf_counter();cpu=process_time()
            d=readout(zw,norm,y)
            results.append(dict(id=obs['id'],sources=ns,receivers=nr,noise_pct=obs['noise_pct'],d_lower=d,
                   d_relative=d/np.linalg.norm(y),radius=radius,ratio=d/radius,excludes=d>radius,
                   wall_seconds=perf_counter()-t,cpu_seconds=process_time()-cpu))
        zrows=[r for r in results if r['sources']==ns and r['receivers']==nr]
        print(ns,nr,row, 'extinction fractions',min(r['d_relative'] for r in zrows),max(r['d_relative'] for r in zrows),flush=True)
    (HERE/'extinction_geometry.json').write_text(json.dumps(geos,indent=2)+'\n')
    if results:
      with (HERE/'extinction.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=results[0]);w.writeheader();w.writerows(results)
      np.savez_compressed(HERE/'fixed_port_weights.npz',**wv)
if __name__=='__main__':main()
