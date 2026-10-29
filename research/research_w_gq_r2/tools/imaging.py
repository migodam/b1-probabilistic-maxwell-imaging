"""Fixed finite unknown-macro inference and separate truth scoring, nested data."""
from pathlib import Path
import argparse,json,time,gc,csv,numpy as np
from scipy.special import logsumexp
from dataclasses import replace
from r2.accounting import BudgetBook,plain
from r2.physics import Eval,Target,S1_TX
from r2.worlds import public_wall,micro_wall,target_defs,make_target
ROOT=Path(__file__).resolve().parents[1]
def microfile(p,tex,f,j):return ROOT/f'results/micro/p{p}_n16_t{tex}_f{f:g}_q{j+2}.npz'
def saved_micro(p,tex,f,j):
    q=microfile(p,tex,f,j)
    if q.exists():return np.load(q)['data']
    return np.load(ROOT/f'results/imaging/micro_p{p}_t{tex}_f{f:g}_q{j}.npz')['data']
def bank(phase):
    b=BudgetBook();out=ROOT/'results/imaging';out.mkdir(exist_ok=True);targets=target_defs()[1]
    fs=[1.] if phase=='E1' else [.85,1.15]
    for f in fs:
        w=public_wall(0);ev=Eval(b,phase,f,S1_TX,w.x,w.y,h=.15 if f>1 else .17)
        for g in range(6):
            w=public_wall(g)
            for j,d in enumerate(targets+[dict(pose=[0,2.1],radius=.12,contrast=0.,loss=0.,shape='disk',orientation=0.)]):
                path=out/f'macro_g{g}_f{f:g}_q{j}.npz'
                if path.exists():continue
                s=ev.evaluate(w,make_target(d),label=f'public macro bank g{g} f{f} target{j}')
                np.savez_compressed(path,data=s.data,eps=w.epsinf,ell=w.ell,x=w.x,y=w.y,target=np.array(json.dumps(d)),backward=np.array(json.dumps(s.residual)))
        del ev;gc.collect()
        if phase=='E2':
            w=micro_wall(0,1/16,0);ev=Eval(b,phase,f,S1_TX,w.x,w.y,h=.15 if f>1 else .17)
            for p in range(2):
                for tex in range(2):
                    for j,d in enumerate(targets):
                        path=out/f'micro_p{p}_t{tex}_f{f:g}_q{j}.npz'
                        if path.exists():continue
                        s=ev.evaluate(micro_wall(p,1/16,tex),make_target(d),label=f'nested true/oracle p{p} tex{tex} f{f} target{j}')
                        np.savez_compressed(path,data=s.data,backward=np.array(json.dumps(s.residual)))
            del ev;gc.collect()
    (out/f'BANK_{phase}_SEALED.json').write_text(json.dumps(dict(phase=phase,completed=True,targets=targets,prior='six public macro walls x12target uniformcounting, approximate physical macro class',knownwall='oracle diagnostic separate micro array never online wall prior')))
def generate():
    out=ROOT/'results/imaging';truth=[];ys=[];rng=np.random.default_rng(28102911)
    targets=target_defs()[1]
    for j in range(12):
        p=j%2;tex=(j//2)%2
        mu=np.asarray([saved_micro(p,tex,f,j) for f in [.85,1.,1.15]])
        for n in range(2):
            y=mu+.01*(rng.normal(size=mu.shape)+1j*rng.normal(size=mu.shape))
            ys.append(y);truth.append(dict(parent=j,noise=n,wall_parent=p,texture=tex,target_id=j,target=targets[j]))
    np.save(out/'OBSERVATIONS.npy',np.asarray(ys));(out/'TRUTH_OFFLINE.json').write_text(json.dumps(truth,indent=2))
    (out/'OBSERVATION_META.json').write_text(json.dumps(dict(nested=True,packing='f=.85/1/1.15,source,RX;complex real/imag sigma.01 independently',public_prior='6x12 independent uniformcounting',seed=28102911)))
def infer():
    out=ROOT/'results/imaging';ys=np.load(out/'OBSERVATIONS.npy');targets=target_defs()[1]
    fs=[.85,1.,1.15]
    data=np.asarray([[[np.load(out/f'macro_g{g}_f{f:g}_q{j}.npz')['data'] for f in fs] for j in range(12)] for g in range(6)])
    empty=np.asarray([[np.load(out/f'macro_g{g}_f{f:g}_q12.npz')['data'] for f in fs] for g in range(6)])
    post=[];est=[];labels=[]
    for case,y in enumerate(ys):
        for acq in ['A0','A1','A2']:
            def sl(a):
                return a[...,1,:2,:] if acq=='A0' else a[...,1,:,:] if acq=='A1' else a
            ll=-.5*np.sum(abs(sl(data)-sl(y))**2,axis=tuple(range(2,sl(data).ndim)))/.01**2
            joint=np.exp(ll-logsumexp(ll));marginal=joint.sum(axis=0)
            wg=int(np.argmax(joint.sum(axis=1)))
            wl=-.5*np.sum(abs(sl(empty)-sl(y))**2,axis=tuple(range(1,sl(empty).ndim)))/.01**2
            g2=int(np.argmax(wl));gmap,jmap=np.unravel_index(np.argmax(joint),joint.shape)
            methods={'PRIOR_ONLY':(np.ones(12)/12,-1),'FINITE_MACRO_TARGET_MARGINAL':(marginal,-1),'JOINT_FINITE_MAP':(np.eye(12)[jmap],gmap),'MEASUREMENT_SELECTED_MACRO':(joint[wg]/joint[wg].sum(),wg),'WALL_FIRST_FIXED':(np.exp(ll[g2]-logsumexp(ll[g2])),g2)}
            for name,(pp,g) in methods.items():
                post.append(pp);est.append(dict(case=case,contract=acq,method=name,wall_id=g,map_id=int(np.argmax(pp)),posterior_mean_position=np.asarray([t['pose'] for t in targets]).T@pp));labels.append(name)
    np.savez_compressed(out/'ONLINE_OUTPUTS.npz',posterior=np.asarray(post),rows=np.array(json.dumps(est,default=plain)))
    (out/'ONLINE_SEALED.json').write_text(json.dumps(dict(time_unix=time.time(),no_truth_input=True,source='infer() only OBSERVATIONS/public prior/macro paid responses',PDE_speed_evidence=False,truth_isolation='separate CLI stage function inputs, not OS sandbox')))
def iou(a,b):
    xx,yy=np.meshgrid(np.linspace(-1,1,401),np.linspace(1.5,2.7,301));points=np.column_stack((xx.ravel(),yy.ravel()))
    def mask(t):
        q=points-np.asarray(t['pose']);ang=t.get('orientation',0.);rot=np.array([[np.cos(ang),np.sin(ang)],[-np.sin(ang),np.cos(ang)]]);q=q@rot.T
        return (q[:,0]**2+(q[:,1]/(2/3 if t['shape']=='ellipse' else 1.))**2)<=t['radius']**2
    ma,mb=mask(a),mask(b);return float(np.sum(ma&mb)/max(1,np.sum(ma|mb)))
def score():
    out=ROOT/'results/imaging'
    if not (out/'ONLINE_SEALED.json').exists():raise RuntimeError('online notsealed')
    truth=json.loads((out/'TRUTH_OFFLINE.json').read_text());targets=target_defs()[1]
    o=np.load(out/'ONLINE_OUTPUTS.npz');records=json.loads(o['rows'].item());rows=[]
    ys=np.load(out/'OBSERVATIONS.npy');oracle=[]
    for case,tt in enumerate(truth):
        p,tex=tt['wall_parent'],tt['texture'];y=ys[case]
        table=np.asarray([[saved_micro(p,tex,f,j) for f in [.85,1.,1.15]] for j in range(12)])
        for acq in ['A0','A1','A2']:
            pred=table[:,1,:2] if acq=='A0' else table[:,1] if acq=='A1' else table
            yy=y[1,:2] if acq=='A0' else y[1] if acq=='A1' else y
            ll=-.5*np.sum(abs(pred-yy)**2,axis=tuple(range(1,pred.ndim)))/.01**2;pp=np.exp(ll-logsumexp(ll));oracle.append(pp)
            records.append(dict(case=case,contract=acq,method='KNOWN_WALL_ORACLE',wall_id=-2,map_id=int(np.argmax(pp)),posterior_mean_position=np.asarray([t['pose'] for t in targets]).T@pp))
    posts=np.r_[o['posterior'],np.asarray(oracle)];positions=np.asarray([t['pose'] for t in targets])
    for rec,pp in zip(records,posts):
        tt=truth[rec['case']];gt=targets[tt['target_id']];estimate=targets[rec['map_id']]
        error=np.linalg.norm(np.asarray(estimate['pose'])-gt['pose']);me=np.linalg.norm(np.asarray(rec['posterior_mean_position'])-gt['pose'])
        material=abs((1+estimate['contrast']+1j*estimate['loss'])-(1+gt['contrast']+1j*gt['loss']))/abs(1+gt['contrast']+1j*gt['loss'])
        eps_err=abs(estimate['contrast']-gt['contrast']);loss_err=abs(estimate['loss']-gt['loss'])
        order=np.argsort(pp)[::-1];hpd=order[np.r_[0,np.cumsum(pp[order])[:-1]]<.95]
        rows.append(dict(**{k:v for k,v in rec.items() if k!='posterior_mean_position'},parent=tt['parent'],noise=tt['noise'],wall_parent=tt['wall_parent'],texture=tt['texture'],position_error=float(error),mean_position_error=float(me),radius_error=abs(estimate['radius']-gt['radius']),eps_real_error=float(eps_err),loss_error=float(loss_err),material_relative_error=float(material),shape_correct=int(estimate['shape']==gt['shape']),IoU=iou(estimate,gt),true_probability=float(pp[tt['target_id']]),HPD95_contains_true=int(tt['target_id'] in hpd),HPD95_size=len(hpd),posterior_position_risk=float(np.sum(pp*np.sum((positions-np.asarray(rec['posterior_mean_position']))**2,axis=1)))))
    with (out/'TARGET_RESULTS.csv').open('w') as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
    summary=[]
    for acq in ['A0','A1','A2']:
        for name in sorted(set(r['method'] for r in rows)):
            rr=[r for r in rows if r['contract']==acq and r['method']==name]
            summary.append(dict(contract=acq,method=name,n=len(rr),position_RMSE=float(np.sqrt(np.mean([r['position_error']**2 for r in rr]))),position_p90=float(np.quantile([r['position_error'] for r in rr],.9)),worst=float(max(r['position_error'] for r in rr)),posterior_mean_RMSE=float(np.sqrt(np.mean([r['mean_position_error']**2 for r in rr]))),mean_material_relative=float(np.mean([r['material_relative_error'] for r in rr])),mean_IoU=float(np.mean([r['IoU'] for r in rr])),coverage95=float(np.mean([r['HPD95_contains_true'] for r in rr]))))
    np.savez_compressed(out/'ALL_POSTERIORS_OFFLINE.npz',posterior=posts,rows=np.array(json.dumps(records,default=plain)))
    (out/'RESULTS.json').write_text(json.dumps(dict(summary=summary,scope='finite target grid, continuous microstructure truth outside public macro likelihood model; not physically exact reference posterior',n_independent_target_parents=12,n_independent_wall_parents=2),indent=2));print(json.dumps(summary))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['bank','generate','infer','score','all']);ap.add_argument('--phase',default='E1');a=ap.parse_args()
    if a.stage=='bank':bank(a.phase)
    elif a.stage=='all':
        generate();infer();score()
    else:globals()[a.stage]()
if __name__=='__main__':main()
