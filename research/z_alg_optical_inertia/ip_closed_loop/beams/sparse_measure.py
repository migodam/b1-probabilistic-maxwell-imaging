"""Freeze, then read back a mechanistic sparse-port replay on the exposed cases."""
import argparse
from datetime import datetime, timezone
import json
from shutil import copyfile

import numpy as np
from flint import acb, arb, ctx
from threadpoolctl import threadpool_limits

from beam_common import HERE, atomic_write, clocks, digest, elapsed, geometry, normalize, require_absent, save
from interval_certify import ExactBox, ProofGrid, ball, column, endpoints, joined, lo, norm_upper, tr, trust_direction
from native_paired import green_kernel_roundoff
from native_residual import ResidualBoundBox
from sparse_ports import restore
from sparse_secant import bound as secant_bound


NAMES=('positive_best_one_port','positive_best_two_ports','positive_best_three_ports')


def freeze(side, dest, refinement):
    start=clocks()
    dest.mkdir(exist_ok=True)
    require_absent(dest/'proposal.json',dest/'proposal_arrays.npz',dest/'future_readings.npz')
    base=HERE.parent/'pair'/f'power_{side}'
    for n in ('registration.json','public_classes.npz','old_readings.npz'):
        copyfile(base/n,dest/n)
    confirm=json.loads((base/'proposal.json').read_text())
    suffix='' if refinement=='_direct' else refinement
    stem=(HERE/f'sparse_ports_power_{side}' if not refinement else
          HERE/f'sparse_secant_power_{side}{suffix}')
    sweep=json.loads(stem.with_suffix('.json').read_text())['rows']
    unscaled=np.load(stem.with_suffix('.npz'))
    bank=np.load(dest/'public_classes.npz')
    gram=bank['power']
    design={}
    arrays={}
    for i,main in confirm['rows'].items():
        energy=main['fluence']
        if refinement:
            key=i+'__'
            raw=unscaled[key+'best_two_actions']
            arrays[key+'actions']=normalize(raw,gram,energy)[:,None]
            arrays[key+'receive_direction']=unscaled[key+'best_two_receive_direction']
            names=['real_box_best_two_ports']
            support=[sweep[i]['best']['indices']]
            pre=sweep[i]['best']
        else:
            v=[normalize(unscaled[f'{i}__positive__k{k}__unscaled'],gram,energy) for k in (1,2,3)]
            arrays[i+'__actions']=np.column_stack(v)
            names=list(NAMES)
            support=[sweep[i]['cases']['positive']['supports'][str(k)]['indices'] for k in (1,2,3)]
            pre=None
        design[i]=dict(policies=names,supports=support,real_box_prospective=pre,
                       fluence=energy,old_radius=main['old_radius'],new_radius=main['new_radius'])
    save(dest/'proposal_arrays.npz',**arrays)
    atomic_write(dest/'proposal.json',dict(rows=design,
                    stage='FROZEN_BEFORE_THIS_REMEASUREMENT',
                    cohort_exposure='Previously published in fixed and powered experiments; same old data and frozen material class. This is a mechanistic paired replay.',
                    design_utc=datetime.now(timezone.utc).isoformat(),
                    old_file_digest=digest(np.frombuffer((dest/'old_readings.npz').read_bytes(),dtype=np.uint8)),
                    design_binary_digest=digest(np.frombuffer((dest/'proposal_arrays.npz').read_bytes(),dtype=np.uint8)),
                    **elapsed(start)))


def readout(side,dest,refinement):
    start=clocks()
    prop=json.loads((dest/'proposal.json').read_text())
    actions=np.load(dest/'proposal_arrays.npz')
    old=np.load(dest/'old_readings.npz')
    new=np.load(dest/'future_readings.npz')
    truth=np.load(dest/'future_offline.npz')
    bank=np.load(dest/'public_classes.npz')
    qbank=np.load(HERE/f'paired_replay_power_{side}.npz')
    bounds=json.loads((HERE/f'sparse_ports_power_{side}.json').read_text())['rows']
    ctx.prec=128
    ctx.threads=2
    geo=geometry(side)
    grid=ProofGrid(geo)
    b=ball(geo.B)
    gram=tr(b)*b
    gram=acb(1/gram[0,0].real)*gram
    eps=None if geo.spectrum is None else green_kernel_roundoff(grid)
    rows={}
    data={}
    for i,r in prop['rows'].items():
        tick=clocks()
        prefix=i+'__'
        cells=bank[prefix+'cells']
        center=bank[prefix+'centers']
        width=bank[prefix+'widths']
        fg=bank[prefix+'fg']
        numeric,exact=grid.restricted(cells)
        models=[ExactBox(exact,numeric,c,width) if eps is None
                else ResidualBoundBox(exact,numeric,c,width,eps*len(cells)) for c in center]
        read=new[prefix+'Y']
        xold=old[prefix+'Y']
        exact_y=ball(read)
        v=actions[prefix+'actions']
        if not np.array_equal(v,new[prefix+'V']):
            raise ArithmeticError('Fresh readings differ from prospective sparse input.')
        fine=ball(v)
        m=fine.nrows()
        gs=tr(fine)*gram*fine
        normalized=[acb((arb(r['fluence'])/gs[j,j].real).sqrt())*column(fine,j) for j in range(fine.ncols())]
        old_v=np.eye(m)[:,[0]]
        physical_v=np.column_stack((old_v,v))
        exact_v=joined(ball(old_v),*normalized)
        fields=np.column_stack((xold,read))
        radii=[r['old_radius']]+[r['new_radius']]*v.shape[1]
        witnesses=[]
        for cls,model in enumerate(models):
            c=center[cls]
            var=np.where(fg,-.46,.56)*width
            locus=np.stack((c,c+var,c-var))
            rec=[]
            for w in locus:
                measure,_,_,_=ExactBox.physical_candidate(model,numeric,w,physical_v,fields,radii,exact_vectors=exact_v)
                rec.append(measure)
            witnesses.append(np.array(rec))
            data[prefix+f'witness_materials_{cls}']=locus
            data[prefix+f'witness_residual_over_radius_upper_{cls}']=np.array(rec)
        policies={}
        for j,name in enumerate(r['policies']):
            stat=[]
            for cls in (0,1):
                arr=qbank[prefix+f'positive__{cls}__q_midpoint']
                err=qbank[prefix+f'positive__{cls}__q_entry_radii']
                independent=restore(arr,err)
                bound=trust_direction(independent,m,normalized[j],column(exact_y,j),r['new_radius'])
                wit=int(np.count_nonzero((witnesses[cls][:,0]<=1.)&(witnesses[cls][:,j+1]<=1.)))
                exclusion=bound is not None and bound[0]>0
                secant_margin=None
                if refinement:
                    q=ball(actions[prefix+'receive_direction'])
                    lower,upper,*_=secant_bound(models[cls],normalized[j],q)
                    projected=(tr(q)*column(exact_y,j))[0,0].real
                    noise=arb(r['new_radius'])*norm_upper(q)
                    decision=max((projected-upper).lower(),(lower-projected).lower())-noise
                    secant_margin=lo(decision.lower())
                    exclusion|=bool(decision>0)
                if exclusion and wit:
                    raise ArithmeticError('A material witness contradicts quadratic exclusion.')
                stat.append(dict(excludes=exclusion,witness_count=wit,
                                 state='EXCLUDED' if exclusion else 'WITNESSED' if wit else 'UNRESOLVED',
                                 bound=None if bound is None else bound[0],
                                 real_box_measured_margin_lower=secant_margin))
            frozen=r['real_box_prospective'] if refinement else bounds[i]['cases']['positive']['supports'][str(j+1)]['result']
            policies[name]=dict(states=stat,ports=r['supports'][j],frozen_distance_certificate=frozen)
        real=int(truth[prefix+'label'])
        truth_residual,*_=ExactBox.physical_candidate(models[real],numeric,truth[prefix+'chi'],
                                                       physical_v,fields,radii,exact_vectors=exact_v)
        if max(truth_residual)>=1 or any(v['states'][real]['excludes'] for v in policies.values()):
            raise ArithmeticError('Physical generator audit failed.')
        rows[i]=dict(policies=policies,exact_common_fluence=r['fluence'],
                     original_material_label=real,physical_generator_all_residual_over_radius_upper=truth_residual,
                     numeric_coefficients_energy_intervals=[endpoints(gs[j,j].real) for j in range(gs.nrows())],
                     **elapsed(tick))
        print(i,[(n,[s['state'] for s in v['states']]) for n,v in policies.items()],flush=True)
    atomic_write(dest/'results.json',dict(rows=rows,side=side,
                 stage='ORIGINAL_OPERATOR_REPLAY_WITH_ACTUAL_SPARSE_READINGS',
                 **elapsed(start)))
    save(dest/'verified_materials.npz',**data)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('stage',choices=('freeze','readout'))
    parser.add_argument('side',type=int,choices=(32,96))
    parser.add_argument('--refined',choices=('direct_receive','optimized_receive'))
    args=parser.parse_args()
    folder=HERE/f'sparse_measured_{args.side}{'_refined' if args.refined else ''}'
    suffix=('_receiver' if args.refined=='optimized_receive' else '' if args.refined is None else '_direct')
    (freeze if args.stage=='freeze' else readout)(args.side,folder,suffix)


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
