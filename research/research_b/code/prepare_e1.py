"""Single owner of Research B paid vector-Mie calls; resumable cached outputs."""
from __future__ import annotations
import json,sys
from pathlib import Path
import numpy as np
from scipy.stats import truncnorm
REPO=Path(__file__).resolve().parents[3];sys.path.insert(0,str(REPO))
from src.b1.physics import Forward
from src.b1_joint.material import decode,prior_sample,from_geometry
from src.b1_joint.response import Response,fit
from cost import Cost,ROOT,dump

SOURCE=REPO.parent/'b1-joint-p2/research/joint_p2/results'

def forward(cost):
    cfg=json.loads((REPO/'configs/protocol.json').read_text());cfg['nmax']=14
    backend=Forward(cfg,cost,ROOT/'results/forward_cache')
    def predict(x,k,kind='E1_FULL',nmax=14,implementation='matching'):
        cost.check_time()
        try:d=decode(x)
        except Exception as exc:
            def failed_setup():raise exc
            cost.evaluate('SETUP_FAILURE_'+kind,dict(latent=list(map(float,x)),k=float(k),kernel_not_started=True),failed_setup)
            raise
        return backend.run(-1,x[2],x[3],kind=kind,nmax=nmax,implementation=implementation,
            radii=d['radii'],eps_layers=d['eps_layers'],k=float(k),cache=True)
    return predict

def prepare(cost):
    cfg=json.loads((ROOT/'configs/E1_FROZEN.json').read_text());paid=forward(cost)
    prior=np.load(SOURCE/'COMMON_RESPONSE_BANK.npz');axes=[prior[f'axis_{j}'] for j in range(4)]
    mesh=np.stack([v.ravel() for v in np.meshgrid(*axes,indexing='ij')],axis=1)
    latent=np.array([from_geometry(*v) for v in mesh]);fields={k:prior[f'fields_{k:g}'] for k in (1.6,2.4)}
    for k in cfg['new_frequencies']:
        dest=ROOT/f'results/E1_BANK_{k:g}.npz'
        if dest.exists():
            fields[k]=np.load(dest)['fields'];continue
        rows=[]
        for i,x in enumerate(latent):
            rows.append(paid(x,k,kind='E1_CONTINUOUS_RESPONSE_BANK')[0])
            if i%200==199:print('bank',k,i+1,'charged',cost.count(),flush=True)
        fields[k]=np.array(rows).reshape(5,5,5,5,32)
        np.savez_compressed(dest,fields=fields[k],latent=latent,k=np.array(k))
        print('bank complete',k,'charged',cost.count(),flush=True)
    saved={'frequencies':np.array(cfg['frequencies']),**{f'axis_{j}':v for j,v in enumerate(axes)}}
    for k in cfg['frequencies']:
        saved[f'fields_{k:g}']=fields[k];saved[f'coeff_{k:g}']=fit(axes,fields[k])
    np.savez_compressed(ROOT/'results/E1_RESPONSE_BANK.npz',**saved)
    # Fixed reference object anchors, never candidate/truth-dependent noise.
    anchor=from_geometry(0,0);sigma=[];anchor_fields=[]
    old_sigma=np.load(SOURCE/'FRESH_FINAL_OBSERVATIONS.npz')['sigma']
    for k in cfg['frequencies']:
        if k in (1.6,2.4):
            j=0 if k==1.6 else 1
            field=np.load(SOURCE/'forward_cache'/f'forward_{j:05d}.npz')['data']
            sd=old_sigma[0 if k==1.6 else 32]
        else:
            field=paid(anchor,k,kind='E1_FIXED_NOISE_ANCHOR')[0];sd=.12*np.sqrt(np.mean(field**2))
        anchor_fields.append(field);sigma.append(sd)
    np.savez_compressed(ROOT/'results/E1_NOISE_ANCHORS.npz',frequencies=cfg['frequencies'],F=np.array(anchor_fields),sigma=np.array(sigma))
    # Independent holdout probes; their truths are never exposed to policies.
    response=Response(ROOT/'results/E1_RESPONSE_BANK.npz');rng=np.random.default_rng(cfg['independent_response_audit']['seed'])
    auditx=prior_sample(rng,64);audit=[];audit_fields=[]
    for k in cfg['new_frequencies']:
        full=[];pred=response.predict(auditx,(k,))
        sd=sigma[cfg['frequencies'].index(k)]
        for i,x in enumerate(auditx):
            f,a,b,info=paid(x,k,kind='E1_INDEPENDENT_FIELD_AUDIT');full.append(f)
            row=dict(k=k,point=i,relative=float(np.linalg.norm(f-pred[i])/np.linalg.norm(f)),whitened=float(np.linalg.norm((f-pred[i])/sd)),
                boundary_residual=info['boundary_backward_residual'],Qabs=info['Qabs'])
            if i<8:
                high=paid(x,k,kind='E1_N14_N18_CONVERGENCE',nmax=18)[0]
                rec=paid(x,k,kind='E1_RECURSION_AUDIT',nmax=18,implementation='recursive')[0]
                row.update(n14_n18_white=float(np.linalg.norm((f-high)/sd)),recursion_white=float(np.linalg.norm((rec-high)/sd)))
            audit.append(row)
        audit_fields.append(full)
        print('audit complete',k,'maxwhite',max(a['whitened'] for a in audit if a['k']==k),'charged',cost.count(),flush=True)
    dump(ROOT/'results/E1_FORWARD_AUDIT.json',dict(rows=audit,uniform_bound=False,scope='independent finite prior points; new-frequency qualification, not continuum certification'))
    np.savez_compressed(ROOT/'results/E1_FORWARD_AUDIT.npz',latent=auditx,frequencies=cfg['new_frequencies'],F_full=np.array(audit_fields))
    # Fresh simulated population created only after all policy definitions frozen.
    out=ROOT/'results/E1_FRESH_SIMULATOR.npz'
    if not out.exists():
        rng=np.random.default_rng(cfg['fresh_parents']['seed']);truth=prior_sample(rng,32);full=[]
        for i,x in enumerate(truth):
            full.append([paid(x,k,kind='E1_FRESH_TRUTH_RESPONSE')[0] for k in cfg['frequencies']])
        full=np.array(full);gain=truncnorm.ppf(rng.random((32,2)),-2,2,loc=1,scale=.15)
        # Potential outcomes share material/gain across each episode, distinct thermal noises.
        # Index0 supplies the original1.6 observation; indices1..3 are NEW blocks.
        noise=rng.normal(size=(32,2,4,6,32))
        outcomes=gain[:,:,None,None,None]*full[:,None,None]+noise*np.array(sigma)[None,None,None,:,None]
        calibration_sd=old_sigma[0]/np.linalg.norm(anchor_fields[cfg['frequencies'].index(1.6)])
        calibration=gain[:,:,None]+calibration_sd*rng.normal(size=(32,2,3))
        from src.b1_joint.material import void_volume
        context=void_volume(truth)[:,None,None]+.01*rng.normal(size=(32,2,3))
        np.savez_compressed(out,truth=truth,F_truth=full,gain=gain,outcomes=outcomes,calibration=calibration,
            calibration_sd=np.array(calibration_sd),context=context,sigma=np.array(sigma),frequencies=cfg['frequencies'],
            parent=np.arange(32),exposure=np.array('created after E1 policy freeze, hidden from acquisition functions'))
    dump(ROOT/'results/E1_PREPARATION.json',dict(status='COMPLETE',new_physics=cost.count(),bank_frequencies=cfg['frequencies'],fresh_parents=32,noise_episodes_per_parent=2,
        hardware_measurements=0,all_new_observations='independent numerical Maxwell simulation plus thermal-noise draws',historical_target_unchanged=True))
    print('E1 physical preparation complete',cost.summary(),flush=True)

if __name__=='__main__':
    cost=Cost()
    with cost.job('e1-physical-preparation-v1','E1'):prepare(cost)
