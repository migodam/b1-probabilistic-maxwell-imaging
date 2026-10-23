"""ORACLE representation diagnostics, not online reconstruction."""
import sys,json,gc
import numpy as np
from wq.inverse_support import PixelGrid,pixel_mesh
from p0.geometry import wall_values
from .accounting import ROOT,Book,write
from .physics import Forward,PixelModel
from .oracle_reference import reference

def projection(seed,g,subsamples):
    x=(g.x_edges[:-1]+g.x_edges[1:])/2;y=(g.y_edges[:-1]+g.y_edges[1:])/2
    xx,yy=np.meshgrid(x,y);pts=np.c_[xx.ravel(),yy.ravel()];re=0.;loss=0.
    offsets=(np.arange(subsamples)+.5)/subsamples-.5
    for a in offsets:
        for b in offsets:
            er,si=wall_values(pts+np.array([4/g.nx*a,3/g.ny*b]),seed);re=re+er-1;loss=loss+si
    return np.r_[re/subsamples**2/8,loss/subsamples**2/.15]

def main():
    if json.loads((ROOT/'results/R0/RESULTS.json').read_text())['status']!='PASS':raise RuntimeError('R0 failed')
    book=Book('R');rows=[];cfg=json.loads((ROOT/'configs/ORIGINAL_P0.json').read_text())
    for nx,ny,ns,order in [(48,32,1,3),(48,32,8,3),(64,48,8,3)]:
        g=PixelGrid(nx,ny); high=nx>48
        try:mesh=pixel_mesh(g,order=order,max_dofs=110000 if high else 60000,max_factor_bytes=3500000000 if high else 1500000000)
        except MemoryError as e:
            rows.append(dict(grid=[ny,nx],subsamples=ns,status='NOT_RUN_MEMORY',reason=str(e)));write(ROOT/'results/R1/RESULTS.json',rows);continue
        for seed in cfg['truth_seeds']:
            x=projection(seed,g,ns);model=PixelModel(g)
            fs=sorted(set(sum(cfg['frequencies'].values(),[])))
            fw=Forward(model,mesh,fs,book,f'R1/{seed}_{nx}x{ny}_ss{ns}',phase='R1_oracle_projection',identity='ORACLE_DIAGNOSTIC',order=order,max_dofs=110000 if high else 60000,max_factor_bytes=3500000000 if high else 1500000000)
            d,_=fw.evaluate(x)
            for band,frequencies in cfg['frequencies'].items():
                a=np.stack([d[fs.index(f)] for f in frequencies]);r=reference(seed,band)
                row=dict(seed=seed,band=band,grid=[ny,nx],subsamples=ns,d_repr=float(np.linalg.norm(a-r)/.01),occupied_area=float(np.count_nonzero(x[:g.n_pixels]+x[g.n_pixels:])*g.cell_area),minimum_feature='NOT_CERTIFIED_BY_PIXEL_PROJECTION',identity='ORACLE_DIAGNOSTIC')
                rows.append(row);write(ROOT/'results/R1/RESULTS.json',rows);print(json.dumps(row),flush=True)
            del fw;gc.collect()
    write(ROOT/'results/R1/COMPLETE.json',dict(totals=book.totals(),raw_projection_is_not_legal_witness=True))
if __name__=='__main__':main()
