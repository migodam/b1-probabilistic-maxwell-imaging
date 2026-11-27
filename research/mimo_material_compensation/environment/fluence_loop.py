"""Seal coherent waveforms, acquire new specimens, classify, and independently check both data contracts."""
from pathlib import Path
from datetime import datetime,timezone
from fractions import Fraction
import argparse
import secrets
import json
import numpy as np
from threadpoolctl import threadpool_limits
from coherent_modes import HERE, TARGET, DEN, material, read_pair, elapsed, now, fro, rd, up
from band_allocation import outward
from passive_replay import Circuit, VerifiedWave, stored_response, load_cover
from probe import at_load
from flint import ctx

IDS=("unit_415","unit_862")

def dump(path,data):
    with path.open('x') as fp:
        json.dump(data,fp,indent=2,ensure_ascii=False)
        fp.write('\n')

def save(path,values):
    with path.open('xb') as fp:np.savez_compressed(fp,**values)

def stamp():
    return datetime.now(timezone.utc).isoformat()

def construction(side):
    if side==96:
        mats,ix,core,r,_=read_pair()
        old=HERE/'low_energy_96_amplitude2'/'selector.npz'
        k='old_reading_total'
    else:
        _,_,r,_,ix,core,*mats=material(32)
        old=HERE/'acquisition_32'/'old_and_predictions.npz'
        k='old'
    with np.load(old,allow_pickle=False) as f: y=f[k]
    return mats,ix,core,r,y

def choice(side):
    clock=now()
    X=json.loads((TARGET/f'classes_{side}.json').read_text())
    r=Fraction(X['radius'])
    # Prefer the fewest physically distinct pulses at the guaranteed gap.
    rows=[b for b in X['rows'] if b['fullclass_lower_over_2r']>1.04]
    row=min(rows,key=lambda b:(b['shots'],-b['fullclass_lower_over_2r']))
    amp=next(Fraction(j,16) for j in range(1,17) if Fraction(j,16)*Fraction(row['fullclass_lower_over_2r']) > Fraction(104,100))
    high=amp*Fraction(row['unit_high_clean_score_lower'])-r
    low=amp*Fraction(row['unit_low_clean_score_upper'])+r
    if low>=high:raise ArithmeticError('Threshold intervals overlap')
    _,_,_,_,old=construction(side)
    key=row['key']
    with np.load(TARGET/f'actions_{side}.npz',allow_pickle=False) as a:
        factors={field:a[key+'_'+field] for field in ('tx_ticks','rx_ticks','amp_ticks')}
    dir=TARGET/f'n{side}'
    dir.mkdir(exist_ok=True)
    save(dir/'controls.npz',dict(**factors,old_observed=old))
    dump(dir/'proposal.json',dict(side=side,selected=key,shots=row['shots'],complex_reads=row['shots'],detector_radius=float(r),
        frequency_MHz=460,amplitude=float(amp),available_energy_upper=float(amp**2),
        low_noisy_upper=outward(low,True),high_noisy_lower=outward(high,False),beta=.01,
        expected_uniform_relative_distance_lower=outward(amp*Fraction(row['fullclass_lower_over_2r']),False),
        selector='fewest pulses then largest rigorous class margin; first amplitude in sixteenths above 1.04 normalized gap',
        receiver='one coherently combined complex reading per pulse. Eight physically connected transmit feeds and eight calibrated receive feeds, ideal lossless combiner, open probe, fixed joint additive ball.',
        generator='Two independently seeded whole-specification perturbations; 48 new exterior sites; shared real nondispersive material at both bands; complex noise at 0.80r.',
        sealed_utc=stamp(),**elapsed(clock)))
    print('sealed',side,key,float(amp**2),flush=True)

def factors(dir):
    with np.load(dir/'controls.npz',allow_pickle=False) as a:
        ph=lambda t:(t[0].astype(float)+1j*t[1])/DEN
        return ph(a['tx_ticks']),ph(a['rx_ticks']),a['amp_ticks']/DEN,a['old_observed']

def readout(ports,tx,rx,ap,amp):
    return np.sum(rx.conj()*(ports[8:16,:8]@tx),axis=0)*ap*amp

def acquire(side):
    clock=now()
    dir=TARGET/f'n{side}'
    prop=json.loads((dir/'proposal.json').read_text())
    tx,rx,ap,_=factors(dir)
    mats,ix,core,r,_=construction(side)
    # Isolated simulation. Seeds are created after the design is sealed and retained for replay.
    seed=secrets.randbits(120)
    rng=np.random.default_rng(seed)
    order=rng.permutation(2)
    ids=np.unique(np.r_[ix,rng.choice(np.setdiff1d(np.arange(side**2),ix),48,replace=False)])
    circuit=Circuit(side,460,ids)
    measured={}
    scores={}
    truth={}
    for j,(name,idx) in enumerate(zip(IDS,order)):
        baseline=mats[idx]
        actual=baseline.copy()
        delta=rng.uniform(-.0075,.0075,len(ids))
        actual[ids]=np.clip(baseline[ids]+delta,0,4.5)
        if np.max(abs(actual-baseline))>.008:raise ValueError('Specimen outside tolerance')
        s=circuit.solve(actual[ids])
        ports=at_load(s,1.)
        clean=readout(ports,tx,rx,ap,prop['amplitude'])
        g=rng.standard_normal(len(ap))+1j*rng.standard_normal(len(ap))
        measured[name]=clean+.80*r*g/np.linalg.norm(g)
        truth[name]=actual
        scores[name]=dict(label='high' if idx==0 else 'low',region_mean=float(actual[core].mean()))
    save(dir/'new_readings.npz',measured)
    save(dir/'offline_materials.npz',dict(**truth,active_union=ids))
    dump(dir/'offline_key.json',dict(seed=seed,answers=scores))
    dump(dir/'acquisition.json',dict(old_observation='in controls.npz',new_readings=list(measured),sealed_utc=stamp(),
        solver=circuit.count,noise_norm_fraction=.8,**elapsed(clock)))

def decide(side):
    clock=now()
    dir=TARGET/f'n{side}'
    prop=json.loads((dir/'proposal.json').read_text())
    with np.load(dir/'controls.npz',allow_pickle=False) as f: a=f['amp_ticks'].tolist()
    result={}
    with np.load(dir/'new_readings.npz',allow_pickle=False) as f:
        for name in f.files:
            value=sum((Fraction(int(x),DEN)*Fraction(float(z.real)) for x,z in zip(a,f[name])),Fraction())
            if value<=Fraction(prop['low_noisy_upper']): label='low'
            elif value>=Fraction(prop['high_noisy_lower']): label='high'
            else: label='inconsistent with both declared specifications'
            result[name]=dict(label=label,score_lower=outward(value,False),score_upper=outward(value,True),
                              conditional_mean_interval=([0,.010000000000000002] if label=='low' else
                              [1.9899999999999998,2.0100000000000002] if label=='high' else None))
    dump(dir/'decision.json',dict(sealed_utc=stamp(),accessed=['controls.npz::amp_ticks','proposal.json','new_readings.npz'],readings=result,**elapsed(clock)))

def replay(side, write_receipt=True):
    clock=now()
    dir=TARGET/f'n{side}'
    pr=json.loads((dir/'proposal.json').read_text())
    verdict=json.loads((dir/'decision.json').read_text())['readings']
    keys=json.loads((dir/'offline_key.json').read_text())['answers']
    tx,rx,ap,old=factors(dir)
    marg={}
    with np.load(dir/'offline_materials.npz',allow_pickle=False) as f, np.load(dir/'new_readings.npz',allow_pickle=False) as y:
        ids=f['active_union']
        for band in (400,460):
            w=Circuit(side,band,ids)
            # The 128-bit kernel audit also checks the native n-by-n convolution symbol.
            v=VerifiedWave(dict(spectrum=w.symbol,grid_shape=(side,side),S=w.c,E0=w.b,
                 low_rank_left=np.zeros((side**2,1)),low_rank_right=np.zeros((1,side**2))))
            for key in IDS:
                s,e,rec=stored_response(w,v,f[key][ids])
                transfer,src_error,_=load_cover(s,e,1.,0.,0.)
                if band==400:
                    prediction=transfer[8:16,:8]/np.sqrt(8);measurement=old
                    guard=up(src_error+rd(16)*fro(transfer))
                else:
                    prediction=readout(transfer,tx,rx,ap,pr['amplitude']);measurement=y[key]
                    guard=up(pr['amplitude']*src_error+rd(64)*fro(transfer)*(len(ap)+1))
                upper=up((fro(prediction-measurement)+guard+rd(4)*(fro(prediction)+fro(measurement)))/pr['detector_radius'])
                marg.setdefault(key,{})[str(band)]=dict(residual_upper_over_r=upper,**rec)
                if upper>=1:raise ArithmeticError(f'Outside error ball {key} {band}')
    if any(verdict[key]['label']!=keys[key]['label'] for key in IDS):raise ArithmeticError('Specimen misclassified')
    if write_receipt:
        dump(dir/'independent_check.json',dict(side=side,matched_labels=True,offline_scores=keys,
                         native_full_wave_errors=marg,arb_bits=ctx.prec,**elapsed(clock)))
    print('independent',side,{key:(verdict[key]['label'],[marg[key][b]['residual_upper_over_r'] for b in ('400','460')]) for key in IDS},flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('step',choices=('choice','acquire','decide','replay','check'))
    p.add_argument('sides',nargs='+',type=int,choices=(32,96))
    args=p.parse_args()
    with threadpool_limits(1):
        for side in args.sides:
            if args.step=='check':
                replay(side,write_receipt=False)
            else:
                globals()[args.step](side)
