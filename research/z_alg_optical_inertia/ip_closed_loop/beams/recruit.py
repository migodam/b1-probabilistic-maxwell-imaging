"""Build the public candidate classes before any physical readings are acquired."""
import argparse
import numpy as np
from threadpoolctl import threadpool_limits

from beam_common import (HERE, OLD, N_TEST, BACKGROUND_WIDTH, FOREGROUND_WIDTH, clocks,
                         elapsed, atomic_write, require_absent, save, geometry,
                         energy_gram, shapes, profiles, forward, pair_jacobian, herm, eigh)


def recruit_one(geom, gram, index, kind, selector):
    t=clocks()
    cells,fg,bg,chi,xy=shapes(geom.side,index,kind)
    g,s,b=geom.restricted(cells)
    y,field,receive,j,born,_=forward(g,s,b,chi,auxiliary=True)
    funcs=profiles(xy,fg)
    jac=pair_jacobian(field,receive,funcs)
    old=jac[:,OLD,:]
    # Find a finite-dimensional material perturbation with strong future response
    # and weak initial response. All fields here are hypothetical class models.
    new=jac[:,1:,:].reshape(-1,funcs.shape[1])
    protected=herm(np.real(old.conj().T@old+0.0002*new.conj().T@new)) + .000002*(funcs.T@funcs)
    vals,vectors=eigh(herm(np.real(new.conj().T@new)),protected,check_finite=False)
    old_radius=.04*np.linalg.norm(y[:,OLD])
    rows=[]
    opt=None
    for rank in range(1,5):
        v=funcs@vectors[:,-rank]
        v /= max(abs(v))
        for change in (.4,.7,1.0,1.3):
            other=chi+change*v
            y2=forward(g,s,b,other)
            delta=y2-y
            discrepancy=np.linalg.norm(delta[:,OLD])
            cov=delta.conj().T@delta
            cov[OLD]=0
            cov[:,OLD]=0
            maxim=float(np.sqrt(max(0,eigh(herm(cov[1:,1:]),gram[1:,1:],eigvals_only=True)[-1])))
            angular=np.linalg.norm(delta[:,1:],axis=0)/np.sqrt(np.diag(gram)[1:].real)
            rec=dict(direction_rank=rank,peak_change=change,old_relative=discrepancy/old_radius,
                     strongest_coherent_relative=maxim/old_radius,
                     best_angle_relative=float(max(angular)/old_radius),
                     median_angle_relative=float(np.median(angular)/old_radius), 
                     maximum_angle_relative=float(max(angular)/old_radius))
            rows.append(rec)
            if rec['old_relative'] > .87 or np.min(other)<=0:
                continue
            score=maxim/old_radius
            if selector=='mixed':
                score-=.4*float(np.median(angular)/old_radius)
            if opt is None or score>opt[0]:
                opt=score,other,rec
    if opt is None:
        return dict(index=index,kind=kind,selected=False,reason='no old-compatible candidate',trials=rows, **elapsed(t)), None
    _,other,selection=opt
    widths=np.where(bg,BACKGROUND_WIDTH,FOREGROUND_WIDTH)
    centers=np.stack((chi,other))
    crit=dict(peak=[float(max(x)) for x in centers], support=int(fg.sum()), background=int(bg.sum()),
              old_radius=old_radius, candidate_new_noise_radius=.65*old_radius,
              born_relative=float(np.linalg.norm(y-born)/np.linalg.norm(y)),
              fullwave_nominal_old_response_norm=float(np.linalg.norm(y[:,OLD])))
    return (dict(index=index,kind=kind,selected=True, **selection, **crit, trials=rows, **elapsed(t)),
            dict(cells=cells,fg=fg,bg=bg,centers=centers,widths=widths))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('side',type=int,choices=[32,96])
    parser.add_argument('--native-count',type=int,default=2)
    parser.add_argument('--tag',default='main')
    args=parser.parse_args()
    folder=HERE/f'{args.tag}_{args.side}'
    folder.mkdir(exist_ok=True)
    path=folder/'registration.json'
    require_absent(path,folder/'public_classes.npz')
    t=clocks()
    geom=geometry(args.side)
    gram=energy_gram(geom)
    arrays=dict(power=gram,rx_indices=np.arange(16),tx_indices=np.arange(16))
    records=[]
    cases=[(i,('solid','ring','double','austria')[i%4]) for i in range(N_TEST)]
    if args.side==96:
        cases=[(200+i,'ring' if i%2==0 else 'solid') for i in range(args.native_count)]
    for index,kind in cases:
        rec,out=recruit_one(geom,gram,index,kind,'mixed')
        records.append(rec)
        if out:
            arrays.update({f'{index}__{k}':v for k,v in out.items()})
        print(index,kind,{k:round(rec[k],4) for k in
              ('old_relative','strongest_coherent_relative','best_angle_relative','born_relative') if k in rec},flush=True)
    save(folder/'public_classes.npz',**arrays)
    atomic_write(path,dict(side=args.side,stage='REGISTERED_BEFORE_DATA',
                           origin='Independent deterministic poses and synthetic shape families; class pairs recruited by modeled old ambiguity and hypothetical larger new response. No blinded material sample claim.',
                           foreground_contrast_uncertainty_half_width=FOREGROUND_WIDTH,
                           background_contrast_half_width=BACKGROUND_WIDTH,
                           background_support='Two weak unknown lobes with independent per-cell values; exact vacuum outside allowed supports.',
                           receivers='half ring (16 out of 32, 0 through 168.75 deg)',
                           illuminators='half dictionary (16 out of 32, 0 through 168.75 deg)',
                           frequency_hz=400000000., side_m=2., old_indices=[OLD],
                           energy='Fixed incident fluence on full square relative to original beam 0, one new shot', cases=records,
                           **elapsed(t)))


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
