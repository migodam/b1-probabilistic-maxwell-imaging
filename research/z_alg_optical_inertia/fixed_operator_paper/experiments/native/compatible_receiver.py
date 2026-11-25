"""A strictly compatible native material with a non-positive receiving block."""
import json
import time

from flint import arb, ctx
import numpy as np
from threadpoolctl import threadpool_limits

from source_index import CONTROL_ID, HERE, STUDY, ProofGrid, ball, diag, diagonal_times
from run import CONTROL, native_geometry, read
from verify import hi, lo, midpoint

CASE = 'solid_half_1.5'
OLD_COUNT = 2


def clocks():
    return time.perf_counter(), time.process_time()


def elapsed(start):
    return dict(wall_seconds=time.perf_counter()-start[0], cpu_seconds=time.process_time()-start[1])


def norm(a):
    return sum((abs(z)**2 for z in a.entries()), arb(0)).nonnegative_part().sqrt()


def main():
    start=clocks()
    contract=json.loads((HERE/'contract.json').read_text())
    previous=json.loads((HERE/'source_index.json').read_text())
    _,width,weight=next(case for case in contract['cases'] if case[0]==CASE)
    state_proof=next(row for row in previous['rows']
                     if row['box']==CASE and row['object']=='disk' and row['s']==OLD_COUNT)
    # Read disk normalization and geometry only. No previously observed Y enters.
    r=float(read(CONTROL/f'{CONTROL_ID}.json')['radius'])
    with np.load(STUDY/'ip_final_v5/strong_control/strong_ring_arrays.npz',allow_pickle=False) as data:
        ids, material = data['ids'], data['c_solid']
    lower=np.maximum(0.,material-width)
    upper=material+width
    if not np.all((lower<=material)&(material<=upper)):
        raise ArithmeticError('The synthetic generating material is outside the licensed box.')
    geo=native_geometry()
    order=sorted(range(geo.B.shape[1]),key=lambda i:int(f'{i:05b}'[::-1],2))
    pick=order[:OLD_COUNT]
    n,p=len(ids),geo.S.shape[0]
    if n!=previous['n'] or p!=state_proof['receiver_directions']:
        raise ArithmeticError('Native proof dimensions do not match the forward instance.')
    setup=clocks()
    _,(g,s,_)=ProofGrid(geo).restricted(ids)
    rhs_sources=ball(geo.B[np.ix_(ids,pick)])
    matrix=diag([1]*n)-diagonal_times([arb(float(c)) for c in material],g)
    rhs=diagonal_times([arb(float(c)) for c in material],rhs_sources)
    build_time=elapsed(setup)

    solve_start=clocks()
    current=matrix.solve(rhs)
    exact_ports=s*current
    # These actual binary64 complex values are the new old-data matrix.
    y=midpoint(exact_ports)
    error=norm(exact_ports-ball(y))
    margin=arb(r)-error
    if not margin>0:
        raise ArithmeticError('Native material did not strictly realize the old observation.')
    state_defect=norm(matrix*current-rhs)
    data_norm=norm(ball(y))
    finish_time=elapsed(solve_start)

    # This state proof is independent of the archived receiver blocks and their
    # observed fields. Only source-only counts are carried over from the archive.
    physical_negative=state_proof['physical_negative']
    source_pos,source_neg=state_proof['source_inertia']
    nu=state_proof['nu']
    if nu!=0 or source_pos+source_neg!=OLD_COUNT:
        raise ArithmeticError('Borrowed state/source enclosures lack the needed complete signs.')
    # At this witness the material sector is <= 0, and K > 0, so
    # Z_o^* Q Z_o <= E^* E. Therefore P_o = Z_o^* Q Z_o-r^2 I is strictly negative.
    old_gap=arb(r)**2-abs(error)**2
    if not old_gap>0:
        raise ArithmeticError('The strict old negative gap could not be established.')
    # One elimination of the augmented form gives K_mat, source-only block,
    # and I_p; the other gives K, P_o, and R_y. Both preceding pivots are
    # nonsingular and complete, so the receiver count includes zero nullity.
    left=(n-physical_negative+source_pos+p,physical_negative+source_neg)
    positive=left[0]-n
    negative=left[1]-OLD_COUNT
    if positive+negative!=p or negative!=state_proof['deficit_source_only']:
        raise ArithmeticError('Conditional receiver dimension calculation mismatch.')

    out=dict(status='STRICT_COMPATIBILITY_AND_RECEIVER_INERTIA_BY_CONGRUENCE',
             geometry='native 400 MHz, 96×96 over 2 m square, vacuum outside stated disk',
             operator_contract=contract['operator_contract'],interval_precision_bits=ctx.prec,
             arb_threads=ctx.threads,blas_threads=1,box=CASE,half_width=width,sector_weight=weight,
             material_source='ip_final_v5/strong_control/strong_ring_arrays.npz:c_solid',
             material_min_max_on_disk=[float(min(material)),float(max(material))],
             lower_endpoint_was_clipped=bool(np.any(material-width<0)),
             physical_active_cells=n,full_grid_cells=geo.side**2,receiver_count=p,
             source_indices=pick,physical_forward_source_columns=OLD_COUNT,scalar_complex_readings=int(y.size),
             old_data='binary64 rounding of the enclosed center-material forward response',
             rounded_old_receivers_real=y.real.tolist(),rounded_old_receivers_imag=y.imag.tolist(),
             historical_disk_Y_not_used=True,disk_reference_radius_id=CONTROL_ID,
             unchanged_joint_radius=r,observed_norm_interval=[max(0.,lo(data_norm)),hi(data_norm)],
             forward_error_upper=hi(error),strict_radius_margin_lower=lo(margin),
             relative_error_upper=hi(error/arb(r)),state_residual_upper=hi(state_defect),
             state_solve='128-bit Arb complex matrix verified solve',
             carried_source_only_proof=dict(file='source_index.json',state_negative=nu,
                 sector_state_negative=physical_negative,source_only_positive=source_pos,
                 source_only_negative=source_neg,source_deficit=state_proof['deficit_source_only']),
             old_positive=0,old_negative=OLD_COUNT,old_zero=0,old_saturated=True,
             strict_negative_old_pivot_gap_lower=lo(old_gap),
             augmented_positive=left[0],augmented_negative=left[1],augmented_zero=0,
             receiver_positive=positive,receiver_negative=negative,receiver_zero=0,
             receiver_count_method='two exact congruences using the same source-only receipts and the new strict compatibility bound',
             setup=build_time,forward_and_verification=finish_time,total=elapsed(start))
    target=HERE/'compatible_receiver.json'
    tmp=target.with_suffix('.tmp')
    tmp.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    tmp.replace(target)
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    ctx.prec=128
    ctx.threads=1
    with threadpool_limits(limits=1):
        main()
