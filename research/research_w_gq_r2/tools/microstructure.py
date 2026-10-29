"""Resolved finite microstructure worlds. All targets fixed before measurements."""
from pathlib import Path
import json,gc,time,argparse,numpy as np
from r2.accounting import BudgetBook,plain
from r2.physics import Eval,S1_TX
from r2.worlds import micro_wall,target_defs,make_target,white
ROOT=Path(__file__).resolve().parents[1]
def key(p,n,tex,f,j):return f'p{p}_n{n}_t{tex}_f{f:g}_q{j}'
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--validate',action='store_true');a=ap.parse_args()
    b=BudgetBook();out=ROOT/'results/micro';out.mkdir(exist_ok=True)
    train,held=target_defs();targets=[*train,*held];rows=[]
    if not a.validate:
        for n in [8,12,16]:
            wa=micro_wall(0,1/n,0)
            for f in [1.,.85,1.15]:
                ev=Eval(b,'E1',f,S1_TX,wa.x,wa.y,h=.15 if f>1 else .17)
                tasks=(range(len(targets)) if n==16 else range(4)) if f==1 else [1]
                for p in range(2):
                    for tex in range(2):
                        wall=micro_wall(p,1/n,tex)
                        np.savez_compressed(out/f'WALL_p{p}_n{n}_t{tex}.npz',x=wall.x,y=wall.y,eps=wall.epsinf,ell=wall.ell)
                        for j in tasks:
                            path=out/(key(p,n,tex,f,j)+'.npz')
                            if path.exists():continue
                            before=b.snapshot();started=time.process_time()
                            if p==0 and n==16 and tex==0 and f==1 and j==1:
                                old=np.load(ROOT/'results/preflight/BENCHMARK.npz')
                                data=old['data'];fields=dict(u=old['u'],mesh_p=old['mesh_p'],mesh_t=old['mesh_t']);meta=dict(reused_paid_benchmark=True)
                            else:
                                s=ev.evaluate(wall,make_target(targets[j]),label=key(p,n,tex,f,j))
                                data=s.data.copy();fields=dict(u=s.u,mesh_p=ev.p.mesh.p,mesh_t=ev.p.mesh.t) if j==1 and f==1 else {}
                                meta=dict(backward=s.residual,ndofs=ev.p.ndofs,timings=s.timings)
                            np.savez_compressed(path,data=data,**fields,target=np.array(json.dumps(targets[j])),metadata=np.array(json.dumps(meta,default=plain)))
                            row=dict(parent=p,n=n,eta=1/n,texture=tex,f=f,target=j,path=str(path.relative_to(ROOT)),CPU=time.process_time()-started,cost_before=before,cost_after=b.snapshot())
                            with (out/'ACTIONS.jsonl').open('a') as log:log.write(json.dumps(row,default=plain)+'\n')
                            print(key(p,n,tex,f,j),flush=True)
                del ev;gc.collect()
    for n in [8,12,16]:
        for p in range(2):
            wa,wb=[micro_wall(p,1/n,t) for t in range(2)]
            diff=wa.epsinf-wb.epsinf;mask=wa.support|wb.support;material=float(np.sqrt(np.mean(diff[mask]**2)))
            for f in [1.,.85,1.15]:
                tasks=(range(len(targets)) if n==16 else range(4)) if f==1 else [1]
                for j in tasks:
                    paths=[out/(key(p,n,t,f,j)+'.npz') for t in range(2)]
                    if not all(q.exists() for q in paths):continue
                    ds=[np.load(q)['data'] for q in paths]
                    rows.append(dict(parent=p,eta=1/n,n=n,f=f,target=j,role='train' if j<2 else 'heldout',material_RMS_occupied=material,material_RMS_public=float(np.sqrt(np.mean(diff**2))),d_A0=float(np.linalg.norm(ds[0][:2]-ds[1][:2])/.01),d_A1=float(np.linalg.norm(ds[0]-ds[1])/.01)))
    if a.validate:
        vout=ROOT/'results/micro_validation';vout.mkdir(exist_ok=True);vr=[]
        for p in range(2):
            for f in [1.,1.15]:
                rr=[r for r in rows if r['parent']==p and r['n']==16 and r['f']==f and (r['role']=='heldout' if f==1 else True)]
                j=max(rr,key=lambda z:z['d_A1'])['target']
                w=micro_wall(p,1/16,0);ev=Eval(b,'E1',f,S1_TX,w.x,w.y,wall_refine=2,h=.145,targetstep=.1,max_dofs=180000)
                ref=[]
                for tex in range(2):
                    s=ev.evaluate(micro_wall(p,1/16,tex),make_target(targets[j]),label=f'independent micro p{p} f{f} tex{tex}')
                    coarse=np.load(out/(key(p,16,tex,f,j)+'.npz'))['data'];ref.append(s.data.copy())
                    row=dict(parent=p,f=f,target=j,texture=tex,endpoint_error_white=float(np.linalg.norm(s.data-coarse)/.01),backward=s.residual,ndofs=ev.p.ndofs)
                    np.savez_compressed(vout/f'p{p}_f{f}_t{tex}.npz',coarse=coarse,reference=s.data,u=s.u,mesh_p=ev.p.mesh.p,mesh_t=ev.p.mesh.t)
                    vr.append(row)
                vr[-1]['pair_distance_reference']=float(np.linalg.norm(ref[0]-ref[1])/.01)
                del ev,s;gc.collect()
        (vout/'RESULTS.json').write_text(json.dumps(vr,default=plain,indent=2))
    result=dict(rows=rows,scope='finite task set, selected frequency stresses; no continuous uniform family bound',newPDE=True)
    (out/'RESULTS.json').write_text(json.dumps(result,default=plain,indent=2));print(json.dumps(result,default=plain))
if __name__=='__main__':main()
