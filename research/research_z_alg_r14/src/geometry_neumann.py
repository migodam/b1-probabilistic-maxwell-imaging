"""R13 section 5 absolute response disks on the archived frequency-cover groups."""

from collections import Counter
import json
from itertools import combinations_with_replacement
from time import perf_counter, process_time

from flint import acb, acb_mat, arb, ctx
import numpy as np

from data14 import ACTIONS, FrequencyTM, ROOT, load_cases
from geometry_bounds import cmatrix, eye, frobenius, spectral_bounds, upper
from material_controls import GRID, HI, LO

ctx.prec = 128
OUT = ROOT/"results"
ACTIONS_CHECKED = ("one_minus_0.02","one_plus_0.02","split_0.02")
CASES = ("spectral-00","spectral-06","spectral-12","spectral-18","spectral-24","spectral-30",
         "fresh6-00","fresh6-14","confirm6-00","confirm6-05")


def opnorm(x):
    return arb(spectral_bounds(x,False)[1])


def absolute_disks(matrices, green_norm, lo, hi, gamma):
    green, incident, receiver = matrices
    d = len(lo)
    interval_step = (arb(HI)-arb(LO))/GRID
    point = [arb(LO)+(int(a)+int(b))*interval_step/2 for a,b in zip(lo,hi,strict=True)]
    radius = [(int(b)-int(a))*interval_step/2 for a,b in zip(lo,hi,strict=True)]
    material = acb_mat([[point[i] if i==j else 0 for j in range(d)] for i in range(d)])
    radial = acb_mat([[radius[i] if i==j else 0 for j in range(d)] for i in range(d)])
    try:
        q0 = (eye(d)-material*green).inv()
    except (ZeroDivisionError, ValueError):
        return None, None, "center_inverse_open"
    q0_norm = opnorm(q0)
    qr_norm = opnorm(q0*radial)
    q = min((qr_norm*green_norm).upper(),(q0_norm*opnorm(radial*green)).upper())
    if not q<1:
        return None, upper(q), "neumann_open"
    currents = q0*material*incident
    fields = incident+green*currents
    y0 = receiver*currents
    t = receiver*q0
    rv = frobenius(receiver)
    y, disk = [], []
    for col in range(incident.ncols()):
        e = acb_mat([[fields[i,col]] for i in range(d)])
        pushed = min((qr_norm*frobenius(e)).upper(),(q0_norm*frobenius(radial*e)).upper())
        tail = (rv*q/(1-q)*pushed).upper()
        for row in range(receiver.nrows()):
            linear_disk = sum(abs(t[row,i]*e[i,0])*radius[i] for i in range(d)).upper()
            c = y0[row,col]
            y.append(c)
            disk.append(((1+gamma)*(linear_disk+tail)+gamma*abs(c)).upper())
    return (y,disk),upper(q),"bounded"


def pair_gap_squared(a,b):
    if a is None or b is None:
        return arb(0)
    ca,ra = a
    cb,rb = b
    total = arb(0)
    for x,rx,y,ry in zip(ca,ra,cb,rb,strict=True):
        gap = (abs(x-y)-rx-ry).lower()
        if gap>0:
            total += gap**2
    return total.lower()


def study(case, source):
    ref = json.loads((OUT/source/(case.case_id+".json")).read_text())
    g = np.asarray(ref["groups"],dtype=np.int64)
    begin,cpu = perf_counter(),process_time()
    model = FrequencyTM(case.geometry)
    gamma = arb(min(case.gain_radius,case.hard_gain_radius or 1.))
    values = {}
    notes = {}
    for offset in (-.02,.02):
        matrices = tuple(cmatrix(x) for x in model.kernels(offset)[:3])
        norm = opnorm(matrices[0])
        result = [absolute_disks(matrices,norm,l,h,gamma) for l,h in g]
        values[offset] = [v[0] for v in result]
        ratios = [v[1] for v in result if v[1] is not None]
        notes[str(offset)] = dict(counts=dict(Counter(v[2] for v in result)),
                                  q_min=min(ratios,default=None),q_max=max(ratios,default=None))
    pairs = list(combinations_with_replacement(range(len(g)),2))
    sqdist = [sum(int(t)**2 for t in np.maximum(g[i,1]-g[j,0],g[j,1]-g[i,0]))
              for i,j in pairs]
    by_frequency = {d:[pair_gap_squared(v[i],v[j]) for i,j in pairs] for d,v in values.items()}
    result = {}
    for name in ACTIONS_CHECKED:
        action = ACTIONS[name]
        comparison = ref["actions"][name]
        eps2 = (2*arb(comparison["new_noise_radius"]))**2
        possible = []
        for k in range(len(pairs)):
            sep = sum(arb(e)*by_frequency[d][k] for d,e in zip(action.deltas,action.energies,strict=True))
            possible.append(not sep>eps2)
        eligible = [v for keep,v in zip(possible,sqdist,strict=True) if keep]
        upper_d = upper((arb(HI)-arb(LO))/GRID*(arb(max(eligible,default=0))/g.shape[-1]).sqrt())
        result[name] = dict(
            D_upper=upper_d, excluded_pairs=len(possible)-sum(possible),
            correlated_D_upper=comparison["D_upper"],
            correlated_excluded_pairs=comparison["excluded_pairs"],
            bernstein_absolute_D_upper=comparison["absolute_only_D_upper"],
            bernstein_absolute_excluded_pairs=comparison["absolute_excluded"],
            radius=comparison["new_noise_radius"])
    return dict(case_id=case.case_id,cohort=case.cohort,source=source,group_count=len(g),pair_count=len(pairs),
                old_D_upper=ref["zero_new_D_upper"],actions=result,response_ranges=notes,
                resolvent_inverses=2*len(g),wall_s=perf_counter()-begin,cpu_s=process_time()-cpu,
                saved_correlated_menu_cpu_s=ref["cpu_s"])


if __name__=="__main__":
    cases = {c.case_id:c for c in load_cases()}
    records=[]
    for source in ("frequency_covers","frequency_covers_r13"):
        for name in CASES:
            row=study(cases[name],source)
            records.append(row)
            print(source,name,row["group_count"],[(n,round(v["D_upper"],5),round(v["correlated_D_upper"],5),
                  v["excluded_pairs"],v["correlated_excluded_pairs"]) for n,v in row["actions"].items()],
                  row["response_ranges"],round(row["cpu_s"],4),flush=True)
    out=dict(definition="R13 linear term plus norm remainder per complex channel; one full complex gain disk",
             arithmetic="Arb 128-bit, stored binary64 operators, identical archived integer groups and noise radii",
             cases=records)
    (OUT/"geometry_neumann.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
