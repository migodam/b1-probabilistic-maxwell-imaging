"""Choose an additional illumination from a location cover and update that cover."""
from pathlib import Path
from itertools import combinations
from time import perf_counter,process_time
import argparse, json, sys
import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits
HERE=Path(__file__).resolve().parent;MAPS=HERE.parent/'maps'
sys.path.insert(0,str(MAPS))
from common import Geometry,geom_path,pattern,place,save_json,clock,elapsed
from sector_map import SectorGeometry,local_certificate,response_reduction,margin_distances
from certified_map import local_models
NAMES=['uniform','random','best_single','coherent']

def dictionary(g):
    xy=(g.xy+.5)*2/g.side-1.
    angle=np.deg2rad(np.arange(-78.75,79,22.5))
    return np.sqrt(10/32)*np.exp(1j*(2*np.pi*400e6/3e8)*(xy@np.vstack((np.cos(angle),np.sin(angle)))))

def propose(side,oid,shape):
    folder=HERE/oid;folder.mkdir(exist_ok=True)
    path=folder/'proposal.json'
    if path.exists():return
    start=clock(); g=Geometry(geom_path(side))
    result=MAPS/'results'/f'{oid}__b0.03125__primary_sector.npz'
    margin=np.load(result)['margin_sector'];ids=np.flatnonzero(~(margin>0))
    if len(ids)<2:
        save_json(path,dict(stopped=True,survivors=ids.tolist(),**elapsed(start)));return
    reps=ids[np.unique(np.round(np.linspace(0,len(ids)-1,min(20,len(ids)))).astype(int))]
    D=dictionary(g);W=D.conj().T@D/np.vdot(g.B[:,0],g.B[:,0]).real;W=(W+W.conj().T)/2
    g.B=D
    offset,profile=pattern(side,shape); transfers=[]
    for idx in reps:
        c,p=place(side,g.xy[idx],offset,profile)
        m=np.zeros(g.n);m[c]=4*p
        T,_,info=g.forward(m)
        if info:raise RuntimeError('nominal material failed to converge')
        transfers.append(T)
    diffs=np.stack([transfers[i]-transfers[j] for i,j in combinations(range(len(reps)),2)])
    options=list(np.eye(W.shape[0],dtype=complex));norms=[]
    for d in diffs:
        A=d.conj().T@d
        vals,vec=eigh((A+A.conj().T)/2,W,subset_by_index=(len(W)-1,len(W)-1))
        options.append(vec[:,0])
    options=np.column_stack(options)
    options/=np.sqrt(np.einsum('ij,ij->j',options.conj(),W@options).real)
    scores=np.linalg.norm(diffs@options,axis=1).min(axis=0)
    random=np.random.default_rng(np.frombuffer(oid.encode(),dtype=np.uint8)).integers(0,len(W))
    choice=[len(W)//2,random,np.argmax(scores[:len(W)]),np.argmax(scores)]
    V=options[:,choice]; B=D@V
    Yfile=np.load(MAPS/'public'/f'{oid}.npz')
    np.savez_compressed(folder/'proposal.npz',ids=ids,representatives=reps,V=V,B=B,score=scores[choice],
                        fluence=np.einsum('ij,ij->j',V.conj(),W@V).real,
                        center_transfer=np.stack(transfers),directions_deg=np.arange(-78.75,79,22.5))
    save_json(path,dict(stopped=False,policies=NAMES,side=side,shape=shape,oid=oid,old_survivors=len(ids),
                        design_representatives=len(reps),old_radius=float(Yfile['r']),
                        new_radius=float(Yfile['r']/np.sqrt(7)),**elapsed(start)))

def acquire(folder):
    rec=json.loads((folder/'proposal.json').read_text())
    if rec['stopped'] or (folder/'future.npz').is_file():return
    g=Geometry(geom_path(rec['side']));plan=np.load(folder/'proposal.npz');g.B=plan['B']
    actual=np.load(MAPS/'offline'/f"{rec['oid']}.npz")['truth']
    t=clock();clean,current,info=g.forward(actual)
    if info:raise RuntimeError('full-wave new acquisition failed to converge')
    defect=current-actual[:,None]*(g.B+g.apply(current))
    rng=np.random.default_rng([917799,*rec['oid'].encode()])
    noise=rng.normal(size=clean.shape[0])+1j*rng.normal(size=clean.shape[0])
    noise*=.75*rec['new_radius']/np.linalg.norm(noise)
    np.savez_compressed(folder/'future.npz',Y=clean+noise[:,None])
    save_json(folder/'instrument.json',dict(state_defect_fro=float(np.linalg.norm(defect)),right_hand_sides=len(NAMES),**elapsed(t)))

def update(folder):
    rec=json.loads((folder/'proposal.json').read_text())
    if rec['stopped'] or (folder/'update.json').exists():return
    start=clock();side=rec['side'];plan=np.load(folder/'proposal.npz');ids=plan['ids']
    Yold=np.load(MAPS/'public'/f"{rec['oid']}.npz")['Y'];read=np.load(folder/'future.npz')['Y']
    g=Geometry(geom_path(side));g.B=np.column_stack((g.B,plan['B']*np.sqrt(7)))
    sector=SectorGeometry(g,1/32)
    rad=np.nextafter(np.hypot(rec['old_radius'],np.sqrt(7)*rec['new_radius']),np.inf)
    times={};margins=np.full((len(ids),len(NAMES)),-np.inf);ranks=np.full_like(margins,-1,dtype=int)
    subset=np.full(g.n,-1);subset[ids]=np.arange(len(ids))
    scales=2.**np.arange(-12,5)
    src=[np.r_[np.arange(7),7+j] for j in range(len(NAMES))]
    for group,cells,profile in local_models(g,rec['shape']):
        included=np.flatnonzero(subset[group]>=0)
        if not len(included):continue
        for weight in (.125,.5,1.,2.,8.):
            local=local_certificate(sector,cells[0],profile,3.75,4.25,weight)
            if local is None: continue
            for block in np.array_split(included,max(1,int(np.ceil(len(included)/64)))):
                ix=subset[group[block]];prepared=response_reduction(sector,cells[block],local)
                P,y,theta,pe,ye,te,n=prepared
                for j,jdx in enumerate(src):
                    red=(P,y[:,:,jdx],theta[:,jdx,:][:,:,jdx],pe,ye,te,n)
                    Y=np.column_stack((Yold,np.sqrt(7)*read[:,j]))
                    d,nu,_=margin_distances(red,Y,rad,scales)
                    improves=d>margins[ix,j];margins[ix[improves],j]=d[improves];ranks[ix[improves],j]=nu[improves]
    np.savez_compressed(folder/'update.npz',ids=ids,margins=margins,negative_index=ranks)
    # The independent generator is opened for scoring only after the decision arrays are written.
    centre=np.load(MAPS/'offline'/f"{rec['oid']}.npz")['centre']
    true_id=int(centre[1]*side+centre[0])
    save_json(folder/'update.json',dict(old=len(ids),after={name:int(np.sum(~(margins[:,j]>0))) for j,name in enumerate(NAMES)},
              keeps_generating_centre={name:bool(true_id in ids[~(margins[:,j]>0)]) for j,name in enumerate(NAMES)},
              new_error_radius=rec['new_radius'],**elapsed(start)))

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['design','acquire','read']);a.add_argument('side',type=int);a.add_argument('ids',nargs='*')
    cmd=a.parse_args()
    with threadpool_limits(limits=1):
        stock=json.loads((MAPS/'public/acquisition.json').read_text())['observations']
        for key,info in stock.items():
            if info['side']!=cmd.side or info['scene']!='weak' or cmd.ids and key.split('__')[1] not in cmd.ids:continue
            folder=HERE/key
            if cmd.stage=='design':propose(cmd.side,key,info['shape'])
            elif cmd.stage=='acquire':acquire(folder)
            else:update(folder)
            print(cmd.stage,key,flush=True)
