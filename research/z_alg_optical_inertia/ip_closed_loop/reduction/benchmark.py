"""Equal-certificate moving-support benchmark on frozen 32² and 96² TM observations."""
import argparse
import json
import os
from pathlib import Path
import sys
from time import perf_counter, process_time

for var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[var] = '1'

from flint import arb, ctx
import numpy as np
from scipy.linalg import eigvalsh, solve
from scipy.sparse.linalg import LinearOperator, svds
from threadpoolctl import threadpool_info, threadpool_limits

from interval_backend import (BallCache, ball, blocks, diagonal, exact_inverse,
                              lower as ball_lower, norm as ball_norm, pivot_ball,
                              row_scaled, trace_ball, verify_signs)
from moving_trace import adj, bk_counts, herm, physical_update
from moving_trace import TraceCache

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
sys.path.insert(0, str(STUDY/'ip_guided_v4/boxes'))
from run import coordinates, SMALL_ID, load, native_geometry
from verify import ProofGrid


def clock():
    return perf_counter(), process_time()


def elapsed(t):
    w,c = clock()
    return {'wall_seconds': w-t[0], 'cpu_seconds': c-t[1]}


def add(*parts):
    return {k: sum(p[k] for p in parts) for k in parts[0]}


def write(path, record):
    tmp = path.with_suffix('.part')
    tmp.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False)+'\n')
    tmp.replace(path)


def load_model(side, level):
    if side == 32:
        obs = load(STUDY/'ip_guided_v4/boxes/results', [SMALL_ID])[SMALL_ID]
        geom, y, radius = obs.geometry, obs.Y, obs.radius
        ids = np.arange(geom.side**2)
        update_mask = np.linalg.norm(coordinates(side), axis=1) <= .35
        base = np.where(update_mask, level, 1/32)
        name = obs.id
    else:
        geom = native_geometry()
        with np.load(STUDY/'ip_final_v5/strong_control/strong_ring_arrays.npz', allow_pickle=False) as data:
            ids, y, radius = data['ids'], data['Y'], float(data['radius'][0])
        x, yy = np.divmod(ids, side)
        update_mask = (x%2 == 0) & (yy%2 == 0)
        base = np.full(len(ids), level)
        name = 'strong_ring350__pw360_32x32__n02'
    numeric = geom.restricted(ids)
    return geom, ids, numeric, y, radius, base, np.flatnonzero(update_mask), name


def original(g, s, b, y, lo, up, w):
    """Independent original-matrix assembly, in current coordinates."""
    mid = w*(lo+up)/2
    product = w*lo*up
    mg = mid[:, None]*g
    pos = np.flatnonzero(product != 0)
    gpos, bpos = g[pos], b[pos]
    kj = adj(gpos)@(product[pos,None]*gpos)
    hj = adj(gpos)@(product[pos,None]*bpos)
    cj = adj(bpos)@(product[pos,None]*bpos)
    return (herm(adj(s)@s + np.diag(w)-mg-adj(mg)+kj),
            adj(s)@y + mid[:,None]*b-hj, herm(adj(y)@y+cj))


def original_queries(base, trace, affine, candidates, t):
    """Direct original-dimensional comparator caching the same fixed terms."""
    k,h,c = base
    base_border = np.block([[k,-h],[-adj(h),c]])
    border_shift = np.zeros(len(base_border))
    border_shift[len(k):] = t
    base_border[np.diag_indices(len(base_border))] -= border_shift
    trace_and_e = np.column_stack((trace, affine))
    for cand in candidates:
        local = np.ravel(np.column_stack((2*cand['local'],2*cand['local']+1)))
        tr = trace_and_e[local]
        w = cand['weight']
        a = -w*(cand['lo']+cand['up']-cand['base'])/2
        b = w*cand['lo']*cand['up']
        delta = np.zeros((len(local),len(local)))
        k2 = np.arange(0,len(local),2)
        delta[k2,k2+1] = delta[k2+1,k2] = a
        delta[k2+1,k2+1] = b
        full = base_border + adj(tr)@(delta@tr)
        yield bk_counts(full[:len(k),:len(k)]), bk_counts(full)


def trace_numeric(g, b, sites):
    t = np.zeros((2*len(sites),len(g)), dtype=complex)
    t[::2,sites] = np.eye(len(sites))
    t[1::2] = g[sites]
    e = np.zeros((len(t),b.shape[1]), dtype=complex)
    e[1::2] = b[sites]
    return t,e


def prepare_candidates(side, ids, sites, base, weights, count, seed):
    """Distinct moving triples. Two rank-one changes and one rank-two change."""
    rng = np.random.default_rng(seed)
    xy = np.column_stack(np.divmod(ids[sites], side))
    spacing = 1 if side == 32 else 2
    neighbors = np.max(abs(xy[:,None,:]-xy[None,:,:]),axis=2)
    neighborhoods = [np.flatnonzero((d > 0) & (d <= 3*spacing)) for d in neighbors]
    boxes, seen = [], set()
    while len(boxes) < count:
        center = int(rng.integers(len(sites)))
        selected = rng.choice(neighborhoods[center], 2, replace=False)
        triplet = (center, int(selected[0]), int(selected[1]))
        sorted_triplet = tuple(sorted(triplet))
        if sorted_triplet in seen:
            continue
        seen.add(sorted_triplet)
        site = sites[list(triplet)]
        cap = base[site]
        fraction = rng.integers(2,15,size=2)/32
        l = np.r_[cap[:2]*fraction, 0.]
        u = cap - l
        # A third upper endpoint grows. It breaks nesting with the common parent
        # and with every other unique support in this benchmark.
        u[2] += int(rng.integers(1,17))/16
        zeros = np.zeros(len(site))
        w = weights[site]
        update = physical_update(triplet,l,u,zeros,cap,w)
        boxes.append(dict(index=len(boxes), cells=ids[site], local=np.array(triplet),
                          state_ids=site, lo=l, up=u, base=cap, weight=w,
                          update=update, base_lo=zeros))
    return boxes


def endpoints(box, base_lo, base_hi):
    l,u = base_lo.copy(),base_hi.copy()
    l[box['state_ids']],u[box['state_ids']] = box['lo'],box['up']
    return l,u


def full_interval(g,s,b,y,left,right,weight,t, numeric=None):
    p,h,c = pivot_ball(g,s,b,y,left,right,weight)
    border = blocks([[p,-h],[-h.conjugate().transpose(),c-diagonal([t]*c.nrows())]])
    sp = verify_signs(p, numeric[0] if numeric is not None else None)
    bp = verify_signs(border, numeric[1] if numeric is not None else None)
    return sp,bp


def count_at(signs):
    return signs['positive_at_least'], signs['negative_at_least'], signs['unresolved']


def certified(box, cache, interval, t, n, radius, j):
    u = box['update']
    inv = exact_inverse(u,box['lo'],box['up'],box['base_lo'],box['base'],box['weight'])
    state = interval.proof(u,inv,None,cache.reduced(u,None))
    border = interval.proof(u,inv,t,cache.reduced(u,float(t)))
    value = ball_lower((arb(j)*t/arb(radius)**2).sqrt())
    strict = border['positive_at_least'] >= n+j and arb(j)*t > arb(radius)**2
    return state,border,strict,(value if strict else None)


def neumann_probe(g, s, b, y, left, right, radius, exact_g):
    """One-sided failure test for the unweighted centered Neumann norm."""
    del s,b,y,radius
    from scipy.linalg import lu_factor, lu_solve
    mid,r = (left+right)*.5,(right-left)*.5
    a = np.eye(len(g))-mid[:,None]*g
    rg = r[:,None]*g
    factor = lu_factor(a,check_finite=False)
    op = LinearOperator(g.shape, matvec=lambda z: rg@lu_solve(factor,z,check_finite=False),
                        rmatvec=lambda z: lu_solve(factor,adj(rg)@z,trans=2,check_finite=False),dtype=complex)
    _,_,vh = svds(op,k=1,which='LM',tol=1.e-8,rng=169)
    p = ball(lu_solve(factor,adj(vh),check_finite=False))
    gp = exact_g*p
    m = [(arb(float(l))+arb(float(u)))/2 for l,u in zip(left,right,strict=True)]
    rr = [(arb(float(u))-arb(float(l)))/2 for l,u in zip(left,right,strict=True)]
    v = ball_norm(row_scaled(gp,rr)).lower()/ball_norm(p-row_scaled(gp,m)).upper()
    return {'norm_lower':ball_lower(v),'unweighted_norm_contraction_disproved':bool(v>1),
            'scope':'||diag(halfwidth) G (I-diag(midpoint)G)^(-1)||_2; lower norm probe, not a spectral-radius or weighted-contraction claim'}


def freeze_config():
    return json.loads((HERE/'freeze.json').read_text())


def run(side, args):
    spec = freeze_config()
    task = spec['scales'][str(side)]
    suffix = '' if args.count == 1000 else f'_{args.count}'
    result_path = args.output_dir/f'result_{side}{suffix}.json'
    st = clock()
    geom, ids, (g,s,b), y, radius, base, sites, name = load_model(side,task['candidate_base_upper'])
    geometry_time = elapsed(st)
    st = clock()
    weights = np.full(len(g),task['weight'])
    base_lo = np.zeros(len(g))
    j = task['tail_j']
    threshold = float(task['threshold_over_radius_squared']*radius**2)
    boxes = prepare_candidates(side,ids,sites,base,weights,args.count,task['seed'])
    generation = elapsed(st)
    st = clock()
    k,h,c = original(g,s,b,y,base_lo,base,weights)
    tr,aff = trace_numeric(g,b,sites)
    common_orig = elapsed(st)
    st = clock()
    numeric = TraceCache.from_dense(k,h,c,tr,aff)
    float_setup = elapsed(st)
    st = clock()
    float_counts, fixed_times = [], {}
    for i,box in enumerate(boxes,1):
        u = box['update']
        float_counts.append((numeric.counts(u,None),numeric.counts(u,threshold)))
        if i in (100,args.count):
            fixed_times[str(i)] = elapsed(st)
    st = clock()
    fixed_base_counts = bk_counts(np.block([[k,-h],[-adj(h),c-threshold*np.eye(len(c))]]))
    fixed_base = {'state':numeric.pivot_counts,'border':fixed_base_counts,
                  'would_exclude_with_same_numeric_rule':fixed_base_counts[0]>=len(g)+j}
    print(f'{side} numeric base {fixed_base} results {sum(bc[0]>=len(g)+j for sc,bc in float_counts)}/{len(boxes)}',flush=True)
    if args.proposals_only:
        return
    grid = ProofGrid(geom)
    _,(eg,es,eb) = grid.restricted(ids)
    ey = ball(y)
    interval_input_time = elapsed(st)
    st = clock()
    ek,eh,ec = pivot_ball(eg,es,eb,ey,base_lo,base,weights)
    et,ee = trace_ball(eg,eb,sites)
    strict_cache = BallCache.from_base(ek,eh,ec,et,ee,numeric_pivot=k)
    strict_setup = elapsed(st)
    st = clock()
    base_proof = verify_signs(strict_cache.source-diagonal([arb(threshold)]*y.shape[1]),
                              numeric.source-threshold*np.eye(y.shape[1]))
    base_proof['positive_at_least'] += strict_cache.pivot_proof['positive_at_least']
    base_proof['negative_at_least'] += strict_cache.pivot_proof['negative_at_least']
    base_proof['strict_exclusion_at_this_threshold'] = (
        base_proof['positive_at_least'] >= len(g)+j and arb(j)*arb(threshold) > arb(radius)**2
    )
    base_verification_time = elapsed(st)
    st = clock()
    interval_counts, strict_times = [], {}
    for i,box in enumerate(boxes,1):
        interval_counts.append(certified(box,numeric,strict_cache,arb(threshold),len(g),radius,j))
        if i in (100,args.count):
            strict_times[str(i)] = elapsed(st)
    print(f'{side} verified {sum(v[2] for v in interval_counts)}/{len(boxes)}, starting direct replay',flush=True)
    if any((count_at(sc)!=npair[0] or count_at(bc)!=npair[1]) and not (sc['unresolved'] or bc['unresolved'])
           for (sc,bc,strict,value),npair in zip(interval_counts,float_counts,strict=True)):
        raise AssertionError('one-sided exact and numerical moving trace disagree')
    # Store a prefix before the dense same-rule reference. It can be recovered from the log if interrupted.
    direct_times = {}
    st = clock()
    direct = []
    for i,pair in enumerate(original_queries((k,h,c),tr,aff,boxes,threshold),1):
        if pair != float_counts[i-1] and not (any(pair[q][2] or float_counts[i-1][q][2] for q in (0,1))):
            raise AssertionError(f'original and reduced numerical disagreement at {i-1}')
        direct.append(pair)
        if i in (100,args.count):
            direct_times[str(i)] = elapsed(st)
            print(f'{side} completed {i} original numerical queries in {direct_times[str(i)]}',flush=True)
    st = clock()
    replay_ids = sorted(set(i for i in task['independent_samples'] if i < args.count))
    replay = []
    for i in replay_ids:
        node_clock = clock()
        box = boxes[i]
        l,u = endpoints(box,base_lo,base)
        nk,nh,nc = original(g,s,b,y,l,u,weights)
        nb = np.block([[nk,-nh],[-adj(nh),nc-threshold*np.eye(len(nc))]])
        a,p = full_interval(eg,es,eb,ey,l,u,weights,arb(threshold),(nk,nb))
        fixed_a,fixed_b,*_ = interval_counts[i]
        if not (a['unresolved'] or p['unresolved'] or fixed_a['unresolved'] or fixed_b['unresolved']):
            if count_at(a) != count_at(fixed_a) or count_at(p) != count_at(fixed_b):
                raise AssertionError(f'original and reduced interval disagreement at {i}')
        replay.append(dict(index=i,state=a,border=p,independent_full_counts_agree=(
            count_at(a)==count_at(fixed_a) and count_at(p)==count_at(fixed_b)
            and not (a['unresolved'] or p['unresolved'])),cost=elapsed(node_clock)))
        print(f'{side} original interval replay box {i}',flush=True)
    replay_time = elapsed(st)
    st = clock()
    neumann_id = next((i for i in replay_ids if interval_counts[i][2]), 0)
    l,u = endpoints(boxes[neumann_id],base_lo,base)
    probe = neumann_probe(g,s,b,y,l,u,radius,eg)
    probe['index'] = neumann_id
    neumann_time = elapsed(st)
    per_box = []
    for box,npair,fpair,(sc,bc,strict,value) in zip(boxes,float_counts,direct,interval_counts,strict=True):
        state_sizes = [int(len(g)), int(len(g)+y.shape[1]), len(box['update'].trace_rows), len(box['update'].trace_rows)+y.shape[1]]
        per_box.append(dict(index=box['index'],global_cells=box['cells'].tolist(),
                            bounds=np.column_stack((box['lo'],box['up'])).tolist(),
                            update_dimension=box['update'].full_block_dimension,
                            update_rank=len(box['update'].trace_rows), dimensions=state_sizes,
                            numeric_direct=fpair,numeric_moving=npair,verified_state=sc,verified_border=bc,
                            strict_exclusion=strict,residual_over_radius_lower=value))
    summaries = {}
    for key,tf in fixed_times.items():
        count = int(key)
        ta,td = strict_times[key],direct_times[key]
        sampled_costs = [r['cost'] for r in replay if r['index'] < count]
        selected = add(*sampled_costs) if sampled_costs else dict(wall_seconds=0,cpu_seconds=0)
        summary = dict(unique_endpoint_vectors=count,unique_changed_supports=count,
                       strict_exclusions=sum(row['strict_exclusion'] for row in per_box[:count]),
                       verified_full_bordered_counts=sum(row['verified_border']['unresolved']==0 for row in per_box[:count]),
                       numeric_full_agreement=sum(row['numeric_direct']==row['numeric_moving'] for row in per_box[:count]),
                       numerical_original=td,numerical_moving=tf,verified_moving=ta,
                       wall_speed_ratio_same_numeric_rule=td['wall_seconds']/tf['wall_seconds'],
                       exact_moving_with_own_cold_setup=add(float_setup,interval_input_time,strict_setup,base_verification_time,ta),
                       total_own_cold_moving_with_common_io=add(geometry_time,generation,common_orig,float_setup,interval_input_time,strict_setup,base_verification_time,ta),
                       selected_independent_replay=selected,
                       total_verified_moving_with_common_io_and_selected_replay=add(
                           geometry_time,generation,common_orig,float_setup,interval_input_time,strict_setup,
                           base_verification_time,ta,selected),
                       total_cold_original_with_common_io=add(geometry_time,generation,common_orig,td))
        summaries[key] = summary
    result = dict(scope='same fixed physical sector and one-threshold state/border sign rule; exact inequalities on exact discrete model',
                  observation=name, radius=radius, side=side, original_grid_cells=side**2,
                  current_dimension=len(g), source_channels=y.shape[1], receiver_channels=y.shape[0],
                  trace_bank_cells=len(sites), known_exterior_vacuum=(side==96), task=task,
                  base=fixed_base, base_interval_pivot=strict_cache.pivot_proof,
                  base_interval_border=base_proof, base_verification_time=base_verification_time,
                  tail_j=j,threshold=threshold,
                  threshold_test_evaluations_per_node=1, interval_precision_bits=ctx.prec,
                  threads=threadpool_info(), shared_operator_and_observation_load=geometry_time,
                  design_batch_once=generation, shared_operator_algebra_setup=common_orig, moving_numeric_setup=float_setup,
                  moving_interval_input_and_green=interval_input_time,moving_interval_setup=strict_setup,
                  prefixes=summaries, interval_original_select_nodes=replay,
                  interval_original_select_nodes_time=replay_time, neumann_probe=probe,neumann_probe_time=neumann_time,
                  boxes=per_box)
    write(result_path,result)
    print(f'{side} result {result_path.name}',flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--side',type=int,choices=(32,96),required=True)
    parser.add_argument('--count',type=int,choices=(100,1000),default=1000)
    parser.add_argument('--output-dir',type=Path,default=HERE)
    parser.add_argument('--proposals-only',action='store_true')
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    ctx.prec,ctx.threads = 128,1
    with threadpool_limits(limits=1):
        run(args.side,args)
