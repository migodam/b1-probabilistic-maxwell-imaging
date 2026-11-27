"""Precompute a one-direction class decision and validate an independent low-power readout."""
from pathlib import Path
from time import perf_counter,process_time
import argparse,json,numpy as np
from threadpoolctl import threadpool_limits
from measurement import models,dep,measured
from probe import HERE, Circuit
from passive_replay import lower,up,fro,rd


def real(z):return np.r_[z.real.ravel(),z.imag.ravel()]
def root(side):return HERE/f'low_energy_{side}'
def freeze(side):
    start=perf_counter(),process_time();p=root(side);p.mkdir(exist_ok=True)
    if (p/'readings.npz').exists():raise FileExistsError('Acquisition exists. Use a clean experiment directory.')
    t,ix,group,xp,xm,r=models(side);beta=.01
    info=json.loads((HERE/f'frequency_classes_{side}.json').read_text())['frequency']['460']
    transfer=json.loads((HERE/'frequency_classes_open.json').read_text())['records']
    item=next(rr for rr in transfer if rr['side']==side)['bands']['460']
    mapped=next(b for b in item['budgets'] if b['beta']==beta)
    row=next(b for b in info['budget'] if b['beta']==beta);data=np.load(HERE/f'frequency_classes_{side}.npz')
    lam=data[row['lambda_key']];assert fro(lam)<=1.
    g0=json.loads((HERE/f'kyp_verified_responses_{side}_{t}_phase_{side}_two_shells_400.json').read_text())['upper_over_2r']
    g1=mapped['open_lower_over_two_r']
    # Predeclared dyadic attenuations in eighths/sixteenths for an equal-duration instrument.
    energy=next(e/16 for e in range(1,17) if lower((e/16)**.5*g1)>1.04)
    def split(index,name):
        q=row['endpoint_terms'][index]
        center=real(data['p_460_'+name]+data['p_460_empty_total'])
        error=up(info[name]['response_error_upper']+item['coefficient_certificate']['empty_total_Fro_error_upper']
                 +mapped['endpoint_terms'][index]['response_error_fro_upper']
                 +rd(300)*(fro(center)+fro(data['p_460_empty_total'])+1))
        loss=up(q['safe_support']+q['derivative_term_error']+q['nonlinear']+error)
        mid=float(lam@center)
        return lower(mid-loss) if index==0 else up_signed(mid+loss)
    a,b=split(0,'plus'),split(1,'minus')
    high_limit=lower(energy**.5*a-r-rd(16)*(abs(a)+r+1))
    low_limit=up_signed(energy**.5*b+r+rd(16)*(abs(b)+r+1))
    if high_limit<=low_limit:raise ArithmeticError('Disjoint data classification failed')
    threshold=(high_limit+low_limit)/2
    np.savez_compressed(p/'selector.npz',direction=lam,threshold=threshold)
    out=dict(side=side,beta=beta,upper400=g0,lower460=g1,necessary_400=lower(1/g0**2),sufficiency_threshold_460=up(1/g1**2),
             ratio_lower=lower((g1/g0)**2),operating_energy=energy,actual_energy_ratio_lower=lower(1/(g0*g0*energy)),
             radius=r,high_read_score_lower=high_limit,low_read_score_upper=low_limit,decision_threshold=threshold,core=t,
             interval_high=[max(0,t-beta),min(4.5,t+beta)],interval_low=[0,beta],shots=8,reads=64,
             scope='Known membership in union of calibrated specifications; common 400 MHz initial data, fixed absolute deterministic receiver ball. No specimen label or added readouts in selection.',
             wall_s=perf_counter()-start[0],cpu_s=process_time()-start[1])
    (p/'proposal.json').write_text(json.dumps(out,indent=2)+'\n')
    print(side,energy,g0,g1,threshold)
def up_signed(x):return float(np.nextafter(x+abs(x)*rd(3),np.inf))
def acquire(side):
    start=perf_counter(),process_time();p=root(side);d=json.loads((p/'proposal.json').read_text())
    t,ix,g,xp,xm,r=models(side);w=Circuit(side,460,ix);xy=np.unravel_index(ix,(side,side))
    vals={};offline={}
    for j,(tag,x) in enumerate([('batch_leaf',xm),('batch_arch',xp)]):
        # Held sample distinct from the full-power pair; support is used only by the simulator.
        change=.006*np.cos(.29*xy[0]+.51*xy[1]+1.3+j)
        actual=np.clip(x+change,0,4.5)
        y0=measured(w.solve(actual),'A_460',0)*np.sqrt(d['operating_energy'])
        rng=np.random.default_rng(20261128+side+9*j);z=rng.normal(size=y0.shape)+1j*rng.normal(size=y0.shape)
        y=y0+(r*.85/np.linalg.norm(z))*z
        vals[tag]=y;offline[tag+'_physical']=actual
    np.savez_compressed(p/'readings.npz',**vals);np.savez_compressed(p/'offline_specimens.npz',**offline)
    (p/'acquisition.json').write_text(json.dumps(dict(frequency_MHz=460,relative_per_sample_noise='0.85 of fixed absolute ball',
        calls=w.count,wall_s=perf_counter()-start[0],cpu_s=process_time()-start[1]),indent=2)+'\n')
def decide(side):
    start=perf_counter(),process_time();p=root(side);d=json.loads((p/'proposal.json').read_text());ys=np.load(p/'readings.npz')
    q=np.load(p/'selector.npz');lam=q['direction']; assert fro(lam)<=1
    res={}
    for tag in ys:
        y=real(ys[tag]);v=float(lam@y);err=up(rd(len(y))*fro(y)*fro(lam))
        if lower(v-err)>d['decision_threshold']:answer='high'
        elif up_signed(v+err)<d['decision_threshold']:answer='low'
        else:answer='unresolved'
        res[tag]=dict(label=answer,ROI_interval=d.get('interval_'+answer,[0,4.5]),score=v,score_abs_error=err)
    (p/'decision.json').write_text(json.dumps(dict(results=res,wall_s=perf_counter()-start[0],cpu_s=process_time()-start[1],
        stop='score with rounding falls on one side of precomputed gap; calibrated two-class contract'),indent=2)+'\n')
def score(side):
    # Postdecision validation includes the actual old data under the chosen perturbation.
    p=root(side);dec=json.loads((p/'decision.json').read_text())['results'];d=json.loads((p/'proposal.json').read_text())
    actual=np.load(p/'offline_specimens.npz');old=np.load(HERE/f'acquisition_{side}'/'old_and_predictions.npz')['old']
    t,ix,g,xp,xm,r=models(side); w=Circuit(side,400,ix); rows={}
    from measurement import VerifiedWave,stored_response
    ver=VerifiedWave(dict(spectrum=w.symbol,grid_shape=(side,side),S=w.c,E0=w.b,low_rank_left=np.zeros((side*side,1)),low_rank_right=np.zeros((1,side*side))))
    for tag,ans in dec.items():
        x=actual[tag+'_physical'];s,e,_=stored_response(w,ver,x);b=fro(s[:24,24]);c=fro(s[24,:24]);gap=lower(1-abs(s[24,24])-e)
        lerr=up(e*(1+(b+c+e)/gap+b*c/gap**2)+rd(64)*(1+fro(s)))
        old_res=up((fro(measured(s,'repeat_A',0)-old)+lerr)/r)
        roi=float(np.mean(x[g==0]));mean_lo,mean_hi=ans['ROI_interval']
        rows[tag]=dict(measured_ROI_native=roi,in_reported_interval=bool(mean_lo<=roi<=mean_hi),old_compatible_upper=old_res,
                       maximum_distance_to_one_public_specification=float(min(max(abs(x-xp)),max(abs(x-xm)))))
        assert rows[tag]['in_reported_interval'] and old_res<=1
    (p/'offline_scores.json').write_text(json.dumps(rows,indent=2)+'\n')
if __name__=='__main__':
    arg=argparse.ArgumentParser();arg.add_argument('op',choices=['freeze','acquire','decide','score']);arg.add_argument('grid',type=int)
    a=arg.parse_args()
    with threadpool_limits(1):globals()[a.op](a.grid)
