"""Select weights by class distance using old readings, then freeze the proposed beam."""
import argparse,json,sys, shutil
from pathlib import Path
from datetime import datetime,timezone
from time import process_time,perf_counter
import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits
from distance import translate,optimize,direction_support,herm
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'beams'))
from quadratic import ResidualForm,old_border
from beam_common import geometry,normalize,save,digest
from importlib.util import spec_from_file_location,module_from_spec
spec=spec_from_file_location('registered_beam_design',HERE.parent/'beams/design.py')
registered=module_from_spec(spec);spec.loader.exec_module(registered)

def library(data,widths,fg,old,ro,origin,a,w,side_sign):
    y,e,d,t,eta,*_=data; m=len(a)
    if eta>=1:raise ValueError('sector base inertia requires a contraction for these specific public boxes')
    form=ResidualForm(widths,e,y,t,d)
    pool=np.flatnonzero(fg)
    power=np.linalg.norm(d,axis=0)*np.linalg.norm(e,axis=1)
    pool=pool[np.argsort(power[pool])]
    cells=np.unique(pool[np.round(np.linspace(0,len(pool)-1,7)).astype(int)])
    buckets={'indefinite':[], 'positive':[]}
    counts=dict(tested=0,saturated_indefinite=0,indefinite_receiver=0,unsaturated=0)
    for mu in 2.0**np.array([-15,-11,-8,-5,-2,1,4,7,10,12]):
        prep=form.prepare(cells,mu,state=True)
        qpsd=herm(prep[3]-prep[2].conj().T@prep[5])
        options=[(0,qpsd,-1,0.)]
        for pos,cell in enumerate(cells):
            for ratio_factor in (1.8,2.05,2.25,2.75,4.,8.):
                ratio=float(ratio_factor/form.row_norm[cell])
                n,q=form.update(prep,pos,ratio)
                if n==1: options.append((n,q,int(cell),ratio))
        for nu,q,cell,ratio in options:
            counts['tested']+=1
            pivot,k,Q=old_border(q,m,old,ro)
            if k!=nu:
                counts['unsaturated']+=1
                continue
            rec_min=float(eigh(Q[m:,m:],eigvals_only=True,check_finite=False,subset_by_index=(0,0))[0])
            if nu:
                counts['saturated_indefinite']+=1
                counts['indefinite_receiver']+=int(rec_min<=0)
            lo,hi=direction_support(Q,a,w)
            signed=(-hi if side_sign==0 else lo)
            bucket=buckets['indefinite' if nu else 'positive']
            metadata=dict(nu=nu,mu=float(mu),cell=cell,ratio=ratio,old_pivot=pivot,
                          receiver_min=rec_min,support_end=None if not np.isfinite(signed) else float(signed))
            bucket.append((signed,Q,metadata))
    for tag,values in buckets.items():
        # A complementary receiver-indefinite pair can work despite two unbounded projections.
        values.sort(key=lambda x:x[0],reverse=True)
        finite=[v for v in values if np.isfinite(v[0])]
        selected=[]
        for v in finite:
            if len(selected)>=2: break
            if not selected or v[2]['mu']!=selected[0][2]['mu']: selected.append(v)
        negative=sorted([v for v in values if v[2]['receiver_min']<0],key=lambda z:abs(z[2]['receiver_min']))
        if negative:selected.append(negative[0])
        buckets[tag]=selected
    return buckets,counts

def propose(source,out,indices):
    p0=perf_counter();c0=process_time()
    src=HERE.parent/'beams'/source; out.mkdir(exist_ok=True)
    for name in ('registration.json','public_classes.npz','old_readings.npz'):
        shutil.copyfile(src/name,out/name)
    reg=json.loads((src/'registration.json').read_text())
    old=np.load(src/'old_readings.npz');bank=np.load(src/'public_classes.npz')
    geometry0=geometry(reg['side'])
    gram=bank['power'];m=len(gram)
    ev,vec=eigh(gram[1:,1:]);basis=np.vstack((np.zeros((1,m-1)),vec/np.sqrt(ev)))
    rows={};arrays={}
    for record in reg['cases']:
        idx=record['index']
        if not record['selected'] or (indices and idx not in indices):continue
        tick=perf_counter();ct=process_time()
        pref=str(idx)+'__'
        cells,centers,widths,fg=[bank[pref+k] for k in ('cells','centers','widths','fg')]
        g,s,b=geometry0.restricted(cells)
        pair=[registered.for_class(g,s,b,ctr,widths) for ctr in centers]
        diff=pair[1][0]-pair[0][0]
        direct,eigen=registered.exact_best_angle(diff,gram)
        radius=record['candidate_new_noise_radius']
        energy=record.get('fixed_fluence',(2.4*radius/np.sqrt(eigen))**2)
        x=normalize(direct,gram,energy); w=diff@x;w/=np.linalg.norm(w)
        origin=(pair[1][0]+pair[0][0])/2
        old_y=old[pref+'Y'].ravel()
        lib=[library(p,widths,fg,old_y,record['old_radius'],origin,x,w,j) for j,p in enumerate(pair)]
        candidates={};aux={}
        for tag in ('positive','indefinite'):
            best=None
            for a1,Q1,weights1 in lib[0][0][tag]:
                for a2,Q2,weights2 in lib[1][0][tag]:
                    c=[translate(q,basis,origin) for q in (Q1,Q2)]
                    est,vector,lam=optimize(*c,m-1)
                    if vector is None:continue
                    if best is None or est>best[0]:
                        best=(est,vector,lam,c,(weights1,weights2),(Q1,Q2))
            if best is None:
                candidates[tag]=dict(found=False);continue
            est,vector,dual,cond,weights,orig=best
            action=normalize(basis@vector,gram,energy)
            # Matrix-valued common linear completions, tested only as diagnostics in binary64 here.
            graph=np.vstack((np.eye(m-1),np.zeros((m,m-1))))
            common=[float(eigh(herm(graph.conj().T@q@graph),eigvals_only=True)[-1]) for q in cond]
            candidates[tag]=dict(found=True,numerical_distance_over_r=float(np.sqrt(max(0.,est)*energy)/radius),
                                 spectral_dual=float(est),dual_weights=dual.tolist(),forms=list(weights),
                                 halfway_graph_max=common)
            aux[tag]=(action,orig)
        actions=[x]
        policies=['known_pair_coherent']
        for tag in ('positive','indefinite'):
            if tag in aux:
                actions.append(aux[tag][0]);policies.append(tag+'_quadratic')
                arrays[pref+tag+'_forms']=np.stack(aux[tag][1])
        arrays[pref+'actions']=np.column_stack(actions)
        rows[str(idx)]=dict(policies=policies,cases=candidates,counts=[r[1] for r in lib],fluence=float(energy),
                            new_radius=radius,old_radius=record['old_radius'],physical_upper_for_centres_over_r=float(np.sqrt(energy*eigen)/radius),
                            elapsed=perf_counter()-tick,cpu=process_time()-ct)
        print(source,idx,[(k,v.get('numerical_distance_over_r')) for k,v in candidates.items()],flush=True)
    save(out/'proposal_arrays.npz',**arrays)
    meta=dict(rows=rows,physical_source=source,design_utc=datetime.now(timezone.utc).isoformat(),
              old_file_digest=digest(np.frombuffer((out/'old_readings.npz').read_bytes(),dtype=np.uint8)),
              design_binary_digest=digest(np.frombuffer((out/'proposal_arrays.npz').read_bytes(),dtype=np.uint8)),
              wall_seconds=perf_counter()-p0,cpu_seconds=process_time()-c0)
    (out/'proposal.json').write_text(json.dumps(meta,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('source');parser.add_argument('indices',type=int,nargs='*')
    args=parser.parse_args()
    with threadpool_limits(limits=1):
        propose(args.source,HERE/args.source,set(args.indices))
