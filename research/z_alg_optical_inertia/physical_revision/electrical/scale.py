"""Physical sizes, incident/receiver angular content and the fixed geometry's reactive spectrum."""
import os
for key in ('OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import csv, json
from pathlib import Path
from time import perf_counter, process_time
import numpy as np
from scipy.fft import fft2, ifft2
from scipy.linalg import eigh, svd, svdvals
from scipy.sparse.linalg import LinearOperator,eigsh

HERE=Path(__file__).resolve().parent
STUDY=HERE.parents[1]; ROOT=STUDY.parents[1]
HIST=ROOT/'research/z_alg_secant_tasks/highdim'
OUT=ROOT/'.artifacts/z_alg_physical_revision'
CAPS=(.25,.5,.75,1.,1.25,1.5,2.,5.)
def spectrum_values(size,apply,limit=48):
    # Only eigenvalues >1/5 matter for the declared contrasts; expand as needed.
    n=size;count=min(limit,n-2)
    while True:
        op=LinearOperator((n,n),matvec=apply,dtype=np.complex128)
        vals=eigsh(op,k=count,which='LA',v0=np.random.default_rng(48912).standard_normal(n),
                   tol=3.e-9,return_eigenvectors=False)
        vals=np.sort(vals)[::-1]
        if vals[-1]<.20 or count>=n-2: return vals
        count=min(n-2,count+28)
def apply_fft(symbol):
    sz=symbol.shape[0];side=sz//2
    def apply(x):
        return ifft2(fft2(x.reshape(side,side),s=(sz,sz))*symbol)[:side,:side].ravel()
    return apply
def effective(values):
    return [int(np.count_nonzero(values>q*values[0])) for q in (.001,.01,.05,.1)]
def analyse(label,side,width,freq,c,normlambda,sm,bm,green,multi,receivers,source_type,truth,doc,cap_list=CAPS):
    k=2*np.pi/normlambda; ka=k*width/np.sqrt(2)
    s,sv,vh=svd(sm,full_matrices=False,check_finite=False);del s
    inc=svdvals(bm,check_finite=False)
    actual_rank=np.count_nonzero(sv>sv[0]*np.finfo(float).eps*max(sm.shape))
    # Rank is numerical at machine threshold. Absolute public ring is full row rank.
    mu=sv[0]**2/100
    vg=vh[:actual_rank]; mul=(1+sv[:actual_rank]**2/mu)**(-.5)-1
    def narrow(x):return x+vg.conj().T@(mul*(vg@x))
    def dark(x):return x-vg.conj().T@(vg@x)
    if callable(green): react,rad=green()
    else:
        h=(green+green.conj().T)/2; rr=(green-green.conj().T)/(2j)
        react=lambda x:h@x;rad=lambda x:rr@x
    n=sm.shape[1];tick=lambda:np.array([perf_counter(),process_time()])
    start=tick()
    eigs=spectrum_values(n,react, min(32,n-2))
    weighted=spectrum_values(n,lambda x:narrow(react(narrow(x))),min(32,n-2))
    hidden=spectrum_values(n,lambda x:dark(react(dark(x))),min(32,n-2))
    # R0 and R spectral radii use the upper end of the Hermitian spectrum.
    seed=np.random.default_rng(4052).standard_normal(n)
    op=LinearOperator((n,n),matvec=lambda x:dark(rad(dark(x))),dtype=complex)
    r0=eigsh(op,k=1,which='LA',v0=seed,tol=1.e-10,return_eigenvectors=False)[0]
    op=LinearOperator((n,n),matvec=rad,dtype=complex)
    rmax=eigsh(op,k=1,which='LA',v0=seed,tol=1.e-10,return_eigenvectors=False)[0]
    rows=[]
    for u in cap_list:
        rows.append(dict(dataset=label,cap=u,body=int(sum(eigs>1/u)),state=int(sum(weighted>1/u)),
                   hidden=int(sum(hidden>1/u)),source_columns=sm.shape[0] if bm is None else bm.shape[1],
                   body_nearest_gap=float(min(abs(1-u*eigs))),weighted_nearest_gap=float(min(abs(1-u*weighted)))))
    item=dict(dataset=label,grid=side,side_in_m=width if freq is not None else '',side_in_wavelengths=width/normlambda,
              wavelength_m=normlambda if freq is not None else '',frequency_GHz=freq if freq is not None else '',
              ka=ka,two_ka=2*ka,grid_area_wavelengths2=(width/normlambda)**2,s=bm.shape[1],p=sm.shape[0],
              physical_ports=receivers,groups=multi,illumination=source_type,contrast=truth,reference=doc,
              effective_source_rank_1pct=effective(inc)[1],effective_receiver_rank_1pct=effective(sv)[1],
              source_ranks_0p1_1_5_10pct='/'.join(map(str,effective(inc))),receiver_ranks_0p1_1_5_10pct='/'.join(map(str,effective(sv))),
              R0_top_numeric=float(r0),R_top_numeric=float(rmax),R0_fraction=float(r0/rmax),mu_norm=mu,
              wall_seconds=(tick()-start)[0],cpu_seconds=(tick()-start)[1])
    print(label,'ka',round(ka,3),'nu body/K/ker at 1',[(r['body'],r['state'],r['hidden']) for r in rows if r['cap']==1],
          'R0',round(r0/rmax,5),flush=True)
    return item,rows
def save(rows,name):
    with (HERE/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
def main():
    inputs=[]; allrows=[]
    def add(*args):
        item,rows=analyse(*args);inputs.append(item);allrows.extend(rows);save(inputs,'geometry.csv');save(allrows,'state_counts.csv')
    g=np.load(HIST/'public_vie_inputs.npz')
    for case,label,W,lam,siz,frequency,desc,truth in [
        ('r10_twod_0001','Annular inclusion',.96,1.,24,None,'12 exterior line sources','0.1 to 12'),
        ('austria_eps4_400MHz','Austria',2.,.75,32,.4,'20 unit plane waves','0 to 3')]:
        S=g[case+'__S']; B=g[case+'__E0']; G=g[str(g[case+'__G_key'].item())]
        add(label,siz,W,frequency,None,lam,S,B,G,1,len(S),desc,truth,'highdim/results/'+case+'.json')
    native=np.load(STUDY/'synthetic/public/fine_geometry.npz')
    sy=native['spectrum']
    for key,label in [('PW360_10x20','Six shapes, full ring'),('PW180_10x20','Six shapes, half ring')]:
        add(label,96,2.,.4,None,.75,native[key+'__S'],native[key+'__B'],lambda:(apply_fft(sy.real),apply_fft(sy.imag)),
            1,20,'10 unit plane waves','0 to shape-specific peak','synthetic/public/fine_geometry.npz')
    man=ROOT/'.artifacts/z_alg_validation/manitoba/prepared_private.npz'
    cal=list(csv.DictReader((STUDY/'validation/manitoba/calibration_margins.csv').open()))
    with np.load(man) as z:
        for idx,freq in enumerate(np.arange(3.,6.01,.5)):
            spec=z[f'f{idx}_spectrum'];S=z[f'f{idx}_S'][1::2];B=z[f'f{idx}_B'][:,::2]
            add(f'Manitoba {freq:g} GHz',40,.14,freq,299792458.,299792458./(freq*1e9),S,B,
                lambda:(apply_fft(spec.real),apply_fft(spec.imag)),2,'24 antennae; 12 x 12 per colour','calibrated line sources',
                'nominal 0 to 1.3 (polyethylene), 0 to 2.03+0.03i (nylon)','validation/manitoba/physical.py')
    # Different masks and input fits make a representative s,p Gram arbitrary for Fresnel.
    # Record sampling and geometry only. Actual common-row groups are documented with each field.
    for freq in (4,6,8):
      for target,ss in [('Foam interior / exterior',8),('Twin dielectric',18)]:
        fname='FoamDielIntTM' if ss==8 else 'FoamTwinDielTM'
        with np.load(STUDY/'results/fresnel'/f'{fname}_{freq}GHz_ports.npz') as f:
            assert np.all(f['mask'].sum(axis=0)==241) and f['mask'].shape==(360,ss)
        lam=299792458./(freq*1e9);ka=2*np.pi*.20/(np.sqrt(2)*lam)
        d=dict.fromkeys(inputs[0],'');d.update(dataset=f'Fresnel {target} {freq} GHz',grid=60,side_in_m=.2,
                side_in_wavelengths=.2/lam,wavelength_m=lam,frequency_GHz=freq,ka=ka,two_ka=2*ka,
                grid_area_wavelengths2=(.2/lam)**2,s=ss,p=241,physical_ports='360 positions, 241 permitted per shot',
                groups='source-dependent mask',illumination='empty-field fitted cylindrical waves',contrast='nominal 0, 0.45, 2',
                reference='fresnel/physical.py')
        inputs.append(d)
    save(inputs,'geometry.csv')
    # Summarise the data-dependent, outward-checked older residuals at unchanged error scales.
    report=[]
    for freq in np.arange(3.,6.01,.5):
      for target in sorted(set(r['target'] for r in cal)):
        chosen=[r for r in cal if r['target']==target and float(r['frequency_GHz'])==freq and r['arm']=='positive_optical']
        above=[float(r['cap']) for r in chosen if float(r['margin_after_both_cross_checks'])>0.]
        r=next(r for r in chosen if float(r['cap'])==1.)
        report.append(dict(frequency_GHz=freq,object=target,highest_above_both_empirical_scales=max(above,default=0.),
                  u1_relative_distance=float(r['verified_lower_relative']),u1_signed_margin=float(r['margin_after_both_cross_checks']),
                  first_calibration=float(r['two_metal_references_relative']),second_calibration=float(r['empty_plus_held_metal_relative'])))
    save(report,'chamber_by_frequency.csv')
if __name__=='__main__':main()
