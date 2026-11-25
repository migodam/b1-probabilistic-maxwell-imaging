"""Held-out receiving arcs for measured two-specimen discrimination."""
import csv, json, sys, warnings
from pathlib import Path
from time import perf_counter,process_time
from itertools import combinations
from zipfile import ZipFile
import numpy as np
HERE=Path(__file__).resolve().parent
STUDY=HERE.parents[1]
sys.path.insert(0,str(STUDY/'validation/manitoba'))
import physical as mb

def score(delta,halves,labels):
    for design in range(2):
        train=delta[halves[design]]; test=delta[halves[1-design]]
        _, _, vh=np.linalg.svd(train,full_matrices=False)
        beam=vh[0].conj()
        port=int(np.argmax(np.linalg.norm(train,axis=0)))
        coherent=float(np.linalg.norm(test@beam)); single=float(np.linalg.norm(test[:,port]))
        yield dict(**labels,design_half=design,train_receivers=len(train),test_receivers=len(test),
                   transmitters=train.shape[1],selected_column=port,design_coherent=float(np.linalg.norm(train@beam)),
                   design_single=float(np.linalg.norm(train[:,port])),held_coherent=coherent,
                   held_single=single,held_ratio=coherent/single)

def main():
    t=perf_counter(); c=process_time();rows=[]
    for freq in (4,6,8):
        observations={}
        for name in ('FoamDielIntTM','FoamDielExtTM','FoamTwinDielTM'):
            with np.load(STUDY/f'results/fresnel/{name}_{freq}GHz_ports.npz') as f:
                cols=[np.flatnonzero(f['source_degrees']==j)[0] for j in (0,180)]
                observations[name]=(f['raw_scattered'][:,cols],f['raw_empty'][:,cols],f['mask'][:,cols].all(axis=1))
        for i,j in combinations(observations,2):
            x,ex,m=observations[i]; y,ey,n=observations[j]
            valid=m&n&np.isfinite(x+y+ex+ey).all(axis=1)
            halves=[valid&(np.arange(360)<180),valid&(np.arange(360)>=180)]
            for cal,delta in [('own_empty',x-y),('common_empty',(x+ex)-(y+ey))]:
                rows.extend(score(delta,halves,dict(facility='Fresnel',frequency_GHz=freq,
                                                    specimen_a=i,specimen_b=j,calibration=cal)))
    with ZipFile(mb.RAW) as archive:
        measured={k:mb.load_file(archive,v)[0] for k,v in mb.FILE_MAP.items()}
    ix=np.ix_(np.arange(1,24,2),np.arange(0,24,2))
    halves=[np.arange(6),np.arange(6,12)]
    for fi,freq in enumerate(mb.FREQUENCIES):
        k=2*np.pi*freq*1.e9/299792458.
        for metal,rad in [('cal_metal_3p5in',.04445),('held_metal_2in',.0254)]:
            exact=mb.pec_scattered(k,rad)[ix]
            for blank in ('empty_first','empty_last'):
                gain=exact/(measured[metal][fi][ix]-measured[blank][fi][ix])
                for i,j in combinations(mb.TARGETS,2):
                    delta=gain*(measured[i][fi][ix]-measured[j][fi][ix])
                    rows.extend(score(delta,halves,dict(facility='Manitoba',frequency_GHz=float(freq),
                                                        specimen_a=i,specimen_b=j,calibration=metal+'_'+blank)))
    with (HERE/'held_ports.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summary={}
    for group in sorted({(r['facility'],r['calibration']) for r in rows}):
        ratios=np.array([r['held_ratio'] for r in rows if (r['facility'],r['calibration'])==group])
        summary[' '.join(group)]=dict(conditions=len(ratios),win=int(sum(ratios>1)),
                                       ratio_quantiles=np.quantile(ratios,[0,.25,.5,.75,1]).tolist())
    (HERE/'results.json').write_text(json.dumps(dict(independent_specimens={'Fresnel':3,'Manitoba':4},summary=summary,
                         wall_seconds=perf_counter()-t,cpu_seconds=process_time()-c),indent=2)+'\n')
    print(summary)
if __name__=='__main__': main()
