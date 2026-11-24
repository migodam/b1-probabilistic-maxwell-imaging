"""Resolve angular interpolation error of the measured response at equal total fluence."""
from pathlib import Path
import csv
import numpy as np
HERE=Path(__file__).resolve().parent
INPUT=HERE.parent/'array/observations'
def clean(s,p,oid):
    read=lambda q: np.load(INPUT/f'optfresh_{oid:02d}__pw360_{s}x{p}__pair{q}pct.npz')['Y']
    # Same frozen noise direction and exactly related amplitudes, up to input rounding.
    return 2*read(1)-read(2)
def interpolate(x,to):
    size=x.shape[1]
    coeff=np.fft.fft(x,axis=1)/size
    orders=np.arange(-(size//2)+1,size//2)
    return coeff[:,orders%size]@np.exp(2j*np.pi*orders[:,None]*np.arange(to)/to)*np.sqrt(size/to)
def ratio(x,y): return float(np.linalg.norm(x-y)/np.linalg.norm(x))
def main():
    arrays=np.load(INPUT/'arrays.npz')
    rows=[]
    for p in (20,32):
      for oid in range(1,7):
        with np.load(HERE.parents[1]/'synthetic/evaluation'/f'optfresh_{oid:02d}.npz') as f:
            active=f['m_fine']>0
        axis=(np.arange(96)+.5)*2/96-1; yy,xx=np.meshgrid(axis,axis,indexing='ij')
        aradius=float(np.max(np.hypot(xx.ravel()[active],yy.ravel()[active])))+np.sqrt(2)/96
        ref=clean(32,p,oid);coeff=np.fft.fft(ref,axis=1,norm='ortho')
        for s in (10,16,24):
            b=arrays[f'B{s}'];b32=arrays['B32'];x=clean(s,p,oid)
            rows.append(dict(shape=oid,receivers=p,from_sources=s,to_sources=32,
                      angular_orders=2*(s//2)-1,active_enclosing_radius_m=aradius,active_two_ka=4*np.pi*aradius/.75,
                      active_B_rel_error=ratio(b32[active],interpolate(b[active],32)),B_rel_error=ratio(b32,interpolate(b,32)),
                      clean_ports_rel_error=ratio(ref,interpolate(x,32)),
                      measured_32_tail_over_order_11=float(np.linalg.norm(coeff[:,12:21])/np.linalg.norm(coeff)),
                      measured_32_tail_over_order_12=float(np.linalg.norm(coeff[:,13:20])/np.linalg.norm(coeff))))
    with (HERE/'angular_interpolation.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    for s in (10,16,24):
        z=[r for r in rows if r['from_sources']==s]; print(s,z[0]['B_rel_error'],min(r['clean_ports_rel_error'] for r in z),max(r['clean_ports_rel_error'] for r in z))
    print('tail',[(min(r[k] for r in rows),max(r[k] for r in rows)) for k in ('measured_32_tail_over_order_11','measured_32_tail_over_order_12')])
if __name__=='__main__':main()
