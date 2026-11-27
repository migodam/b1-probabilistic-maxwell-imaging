"""Validate a connected exact-data path in the supplied six-cell model.

Every interval calculation uses Arb. Floating point only proposes a global
free coordinate and the four-by-four preconditioners. Parameter and gauge
follow the same continuous polygon throughout, including at chart seams.
"""
from pathlib import Path
from time import perf_counter, process_time
import json
import numpy as np
from flint import acb, acb_mat, arb, ctx

HERE = Path(__file__).resolve().parent
N = 6
SETS = [[i for i in range(N) if m & (1 << i)] for m in range(1 << N)]
ctx.prec = 144


def q(a, b=1):
    return arb(a) / b


def upper(x):
    return float(np.nextafter(float(x.upper()), np.inf))


def make_equations():
    # det(I-X(G-B_j C)) - (1+y_j) det(I-XG), where C=2i B_3^T.
    b = [[q(19*v,1000) for v in row] for row in
         [[10,2,8],[8,10,1],[2,7,10],[-5,9,7],[6,-4,10],[-7,-5,3]]]
    x0 = [q(i,10) for i in (17,15,20,13,18,16)]
    g = [[acb([q(8,100), q(6,100), q(25,1000), q(1,100)][min(abs(i-j),N-abs(i-j))],
              sum(b[i][s]*b[j][s] for s in range(3)) + (q(24,100) if i==j else 0))
          for j in range(N)] for i in range(N)]
    coeff = []
    for port in (-1,0,1):
        mat = [[g[i][j] - (acb(0,2)*b[i][port]*b[j][2] if port>=0 else 0)
                for j in range(N)] for i in range(N)]
        coeff.append([(-1)**len(I)*acb_mat([[mat[i][j] for j in I] for i in I]).det()
                      if I else acb(1) for I in SETS])
    at_x = [sum(a[m]*prod(x0[i] for i in I) for m,I in enumerate(SETS)) for a in coeff]
    if at_x[0].contains(0):
        raise ArithmeticError('unverified denominator at target')
    mixed = [[v-(at_x[k]/at_x[0])*w for v,w in zip(coeff[k],coeff[0])] for k in (1,2)]
    return [[z.real for z in row] for row in mixed]+[[z.imag for z in row] for row in mixed], x0


def prod(items):
    out = arb(1)
    for a in items:
        out *= a
    return out


def centered(coeff, center):
    r = list(coeff)
    for i,c in enumerate(center):
        for mask in range(1<<N):
            if not mask & (1<<i):
                r[mask] += r[mask|(1<<i)]*c
    return r


def run():
    wc, pc = perf_counter(), process_time()
    coeff, x0 = make_equations()
    data = json.loads((HERE/'compensation_small.json').read_text())
    knots = {0:x0}
    for s in data['long_continuation']:
        t = round(200*s['t'])
        knots[t] = [x0[0]+q(t,200)]+[arb(float(a)) for a in s['x'][1:]]
    ticks = sorted(knots)
    assert ticks == list(range(min(ticks),max(ticks)+1))
    centers, derivatives = [], []
    for a,b in zip(ticks,ticks[1:]):
        mid = [(x+y)/2 for x,y in zip(knots[a],knots[b])]
        cs = [centered(row,mid) for row in coeff]
        centers.append(cs)
        derivatives.append(np.array([[float(row[1<<i].mid()) for i in range(N)] for row in cs]))
    choices = {k:[i for i in range(1,N) if i!=k] for k in range(1,N)}
    conditions = {k:max(np.linalg.cond(j[:,ix]) for j in derivatives) for k,ix in choices.items()}
    free = min(conditions,key=conditions.get)
    ix = choices[free]
    z = q(1,8192)
    records = []
    for a,b,cs,jmat in zip(ticks,ticks[1:],centers,derivatives):
        d = [(y-x)/2 for x,y in zip(knots[a],knots[b])]
        r = [abs(v)+(z if i in ix else 0) for i,v in enumerate(d)]
        dm = [prod(d[i] for i in I) for I in SETS]
        rm = [prod(r[i] for i in I) for I in SETS]
        pre = np.linalg.inv(jmat[:,ix])
        forcing, contraction = arb(0), arb(0)
        for j,weights in enumerate(pre):
            weights = [arb(float(v)) for v in weights]
            row = [sum(w*v for w,v in zip(weights,vals)) for vals in zip(*cs)]
            # Signed coefficients along the common parameter, before taking absolute values.
            line = [arb(0) for _ in range(N+1)]
            for m,I in enumerate(SETS):
                line[len(I)] += row[m]*dm[m]
            forcing = forcing.max(sum(map(abs,line)))
            lips = arb(0)
            for jj,i in enumerate(ix):
                bit=1<<i
                lips += abs(int(j==jj)-row[bit]) + sum(abs(row[m])*rm[m^bit]
                         for m in range(1<<N) if m!=bit and m&bit)
            contraction = contraction.max(lips)
        low = min(min(x,y)-(z if i in ix else 0) for i,(x,y) in enumerate(zip(knots[a],knots[b])))
        high = max(max(x,y)+(z if i in ix else 0) for i,(x,y) in enumerate(zip(knots[a],knots[b])))
        valid = bool(low>0 and high<q(9,2) and contraction<1 and forcing+contraction*z<z)
        records.append(dict(tick=a,forcing_upper=upper(forcing),contraction_upper=upper(contraction),
                            maps_ball=valid))
    for row in records:
        if not row['maps_ball']:
            break
    result=dict(model='supplied reciprocal dissipative six-cell control',
                real_data_coordinates=4,chart_fixed_exterior_pixel_zero_based=free,
                global_free_coordinate_condition_numbers=conditions,
                precision_bits=ctx.prec,tick_unit='1/200',common_corrector_radius='1/8192',
                all_intervals_verified=all(row['maps_ball'] for row in records),
                number_intervals=len(records),parameter_ticks=[ticks[0],ticks[-1]],
                physical_ROI_endpoints=[(340+ticks[0])/200,(340+ticks[-1])/200],
                verified_width_if_all=(ticks[-1]-ticks[0])/200,
                wall_seconds=perf_counter()-wc,cpu_seconds=process_time()-pc,segments=records)
    (HERE/'same_data_path.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print({k:v for k,v in result.items() if k!='segments'},flush=True)
    print('largest forcing',max(r['forcing_upper'] for r in records),'largest contraction',
          max(r['contraction_upper'] for r in records),'failed',sum(not r['maps_ball'] for r in records))
    return result


if __name__=='__main__':
    run()
