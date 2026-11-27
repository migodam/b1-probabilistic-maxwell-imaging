"""Premeasurement design, isolated synthetic acquisition and independent task bounds."""
from pathlib import Path
import argparse,json,sys
from time import perf_counter,process_time
import numpy as np
from threadpoolctl import threadpool_limits
from probe import HERE,Circuit,at_load,R
from passive_replay import stored_response,VerifiedWave,lower,fro,rd,up

LAYOUT_A=(np.arange(8),np.arange(8,16))
LAYOUT_B=(np.arange(0,16,2),np.arange(1,16,2))
ACTIONS=('repeat_A','swap_B','read_matched_arc','switch_existing_load','switch_new_probe','A_460')
def block(a,layout):tx,rx=layout;return a[np.ix_(rx,tx)]/np.sqrt(8)

def measured(s,action,phase):
    if action=='switch_new_probe':
        return block(at_load(s,np.exp(1j*phase)),LAYOUT_A)
    o=at_load(s,1,range(24))
    if action=='read_matched_arc':
        return block(o,(np.arange(8),np.arange(16,24)))
    if action=='switch_existing_load':
        # Original hold terminal 16 switched to a short, no new geometry.
        o=o-np.outer(o[:,16],o[16,:])/(1+o[16,16])
    return block(o,LAYOUT_B if action=='swap_B' else LAYOUT_A)

def name(side):return HERE/f'acquisition_{side}'
def models(side):
    t={32:2.,96:1.}[side]
    q=np.load(HERE/f'phase_{side}_two_shells.npz')
    return t,q['ix'],q['group'],q[f'plus_{t}'],q[f'minus_{t}'],float(q['r'])
def dep(side,t,band):
    return np.load(HERE/f'verified_responses_{side}_{t}_phase_{side}_two_shells_{band}.npz')

def design(side):
    path=name(side);path.mkdir(exist_ok=True);t,ix,group,xp,xm,r=models(side);p={};data={};clock=perf_counter()
    for band in (400,460):
        states=dep(side,t,band); plus,minus=states['s1'],states['s0']
        for a in ACTIONS:
            if (a=='A_460') != (band==460):continue
            phases=np.arange(256)*2*np.pi/256 if a=='switch_new_probe' else [0.]
            def pred(phi):
                return measured(plus,a,phi),measured(minus,a,phi)
            scored=[np.linalg.norm(pred(phi)[0]-pred(phi)[1]) for phi in phases]
            j=np.argmax(scored);arg=float(phases[j]);val=pred(arg)
            p[a]=dict(band=band,phase_radians=arg,pair_distance_over_2r=float(scored[j]/(2*r)),shots=8,complex_reads=64,
                      sum_available_W=1.,receiver_chains=8)
            data[a+'_high'],data[a+'_low']=val
    selected=max(p,key=lambda a:p[a]['pair_distance_over_2r'])
    old=(data['repeat_A_high']+data['repeat_A_low'])/2
    np.savez_compressed(path/'old_and_predictions.npz',old=old,**data)
    proposal=dict(side=side,core=t,r=r,selected=selected,actions=p,wall_s=perf_counter()-clock,
                 scope='Finite public competing specifications. Scores from fitted public exemplars; uniform separation on neighborhoods requires the separate two-endpoint directional replay.',
                 readout='total calibrated power waves. Empty subtraction is an identical known translation of each same-action pair.')
    if (path/'readings.npz').exists():raise FileExistsError('Use a fresh study copy for a new proposal')
    (path/'proposal.json').write_text(json.dumps(proposal,indent=2)+'\n')
    print(side,selected,p[selected],flush=True)

def acquire(side):
    path=name(side);prop=json.loads((path/'proposal.json').read_text())
    prev=np.load(path/'old_and_predictions.npz');t,ix,g,xp,xm,r=models(side);xys={}
    # Deterministic held specimen perturbations within publicly declared manufacturing tolerance.
    xy=np.unravel_index(ix,(side,side)); field=.0001*np.sin(.43*xy[0]+.32*xy[1]+.2)
    for i,(label,x) in enumerate((('specimen_1',xp),('specimen_2',xm))):
        xys[label]=np.clip(x+field*(1 if i==0 else -1),0,4.5)
    arrays={};records={};tic=perf_counter()
    for band in (400,460):
        w=Circuit(side,band,ix)
        dummy=dict(spectrum=w.symbol,grid_shape=(side,side),S=w.c,E0=w.b,low_rank_left=np.zeros((side*side,1)),low_rank_right=np.zeros((1,side*side)))
        ver=VerifiedWave(dummy)
        for j,(label,x) in enumerate(xys.items()):
            s,err,evidence=stored_response(w,ver,x)
            for i,(a,c) in enumerate(prop['actions'].items()):
                if c['band']!=band:continue
                y=measured(s,a,c['phase_radians'])
                rng=np.random.default_rng(side*103+j*29+i)
                noise=rng.normal(size=y.shape)+1j*rng.normal(size=y.shape)
                noise*=.60*r/np.linalg.norm(noise)
                arrays[label+'_'+a]=y+noise
                if a=='repeat_A':
                    # Error under new-probe open circuit; interval for the exact old discrepancy.
                    b=fro(s[:24,24]);f=fro(s[24,:24]);h=lower(1-abs(s[24,24])-err)
                    lft=up(err*(1+(b+f+err)/h+b*f/h**2)+rd(64)*(1+fro(s)))
                    records[label]=dict(old_feasible_residual_upper_over_r=up((fro(y-prev['old'])+lft)/r),
                                        material_box=[float(min(x)),float(max(x))],true_mean=float(np.mean(x[g==0])))
                if a==prop['selected']:
                    records[label]['point_independent_response_error_25port']=err
            arrays[label+'_material']=x
    assert all(v['old_feasible_residual_upper_over_r']<1 for v in records.values())
    np.savez_compressed(path/'readings.npz',**arrays)
    (path/'sampler_receipt.json').write_text(json.dumps(dict(physical_samples=records,wall_s=perf_counter()-tic,
        truth_access='specimens generated here only after sealed proposal; material arrays scored offline'),indent=2)+'\n')

def update(side):
    path=name(side);prop=json.loads((path/'proposal.json').read_text());r=prop['r'];obs=np.load(path/'readings.npz');pred=np.load(path/'old_and_predictions.npz')
    out={}
    for a in prop['actions']:
        out[a]={}
        for lab in ('specimen_1','specimen_2'):
            y=obs[lab+'_'+a]
            out[a][lab]={i:float(np.linalg.norm(y-pred[a+'_'+i])/r) for i in ('low','high')}
    (path/'outcomes.json').write_text(json.dumps(dict(arms=out,selected=prop['selected'],budget_stop='1 additional 8-pulse block',
        inference='Point discrepancies. Full class conclusion for selected frequency uses the independent spectral-direction lower bounds and the presence of exactly one of the two specification classes.'),indent=2)+'\n')
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('step',choices=['design','acquire','update']);a.add_argument('side',type=int)
    cmd=a.parse_args()
    with threadpool_limits(1):globals()[cmd.step](cmd.side)
