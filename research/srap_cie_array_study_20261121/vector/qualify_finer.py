"""Check the fine-grid simulation against a third grid on the preselected first shell."""
from generate import *
def run():
    st=time.perf_counter();cpu=time.process_time();records=[]
    for geom in ['plane_grid','hemi']:
        fields=[];rows=[]
        for n in [12,24,36]:
            mod=Model(n,geom);m=t(frac(mod,0)).real;d=m.repeat_interleave(3)[:,None];pieces=[];errs=[]
            for k in range(0,16,4):
                x,er=mod.solve(m,d*mod.E[:,k:k+4]);pieces.append(npy(mod.S@x));errs.append(er)
            fields.append(np.concatenate(pieces,axis=1))
            rows.append(dict(n=n,max_linear_residual=max(errs),counts=dict(mod.count)))
        rel=lambda a,b:float(np.linalg.norm(fields[a]-fields[b])/np.linalg.norm(fields[b]))
        records.append(dict(array=geom,case=0,relative_12_24=rel(0,1),relative_24_36=rel(1,2),relative_12_36=rel(0,2),solves=rows))
    write(ROOT/'data'/'shell_finer.json',dict(rows=records,wall=time.perf_counter()-st,cpu=time.process_time()-cpu))
if __name__=='__main__':run()
