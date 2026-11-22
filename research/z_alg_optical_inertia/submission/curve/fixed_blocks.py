"""Cross-evaluate the two outward-certified 4×4 profiles without retuning."""

from itertools import chain
import json
from time import perf_counter, process_time

import numpy as np
from flint import acb, acb_mat, arb, ctx
from scipy.linalg import eigh
from scipy.sparse.linalg import LinearOperator, cg
from threadpoolctl import threadpool_limits

from curve_inputs import HERE, full_ring_inputs
from native import Geometry, Native, Weights, ball, frobenius_squared, gershgorin
from native import high, low, mid, mult_rows, star
from search import h, save_json, times

SINGLE = HERE.parent/'single'
CAP = .75


class SavedProfile(Native):
    """Use the same normalized binary64 weights as the inherited exact proof."""

    def __init__(self,g,profile,cap):
        self.mu=profile['mu']
        self.normalized_real=np.asarray(profile['sector']).ravel()
        self.normalized_optical=np.asarray(profile['optical']).ravel()
        super().__init__(g,Weights(cap,self.mu,tuple(self.mu*self.normalized_real),
                                  tuple(self.mu*self.normalized_optical)))
        self.wr=self.normalized_real[g.blocks]
        self.z=(cap*self.wr+1j*self.normalized_optical[g.blocks])/2

    def exact_scaled_coefficients(self):
        cap=arb(self.w.cap)
        real=[arb(float(x)) for x in self.normalized_real]
        z=[acb(cap*wr,arb(float(opt)))/2 for wr,opt in zip(real,self.normalized_optical,strict=True)]
        return [real[i] for i in self.g.blocks],[z[i] for i in self.g.blocks]


def vertical(a,b):
    return acb_mat(a.nrows()+b.nrows(),a.ncols(),a.entries()+b.entries())


def reusable_compact(p):
    start,cpu=perf_counter(),process_time()
    g=p.g
    rec=len(g.sensor)
    ns=g.s
    rhs=np.column_stack((g.sensor_star/p.mu,p.z[:,None]*g.incident))
    pre=LinearOperator((g.n,g.n),matvec=p.precondition,dtype=np.complex128)
    cols,iterations=[],[]
    for l in rhs.T:
        track=[]
        x,info=cg(p.operator,l,M=pre,rtol=2e-13,atol=0.,maxiter=600,
                  callback=lambda _:track.append(1))
        if info:
            raise ArithmeticError(f'Common RHS CG did not converge at column {len(cols)}: {info}.')
        cols.append(x)
        iterations.append(len(track))
    x=ball(np.column_stack(cols))
    ax=p.exact_action(x)
    _,z=p.exact_scaled_coefficients()
    s=star(g.sensor_ball)/arb(p.mu)
    zb=mult_rows(g.incident_ball,z)
    l=acb_mat(g.n,rec+ns,list(chain.from_iterable(
        [s[i,j] for j in range(rec)]+[zb[i,j] for j in range(ns)] for i in range(g.n))))
    r=ax-l
    m=star(x)*ax-star(x)*l-star(l)*x
    defect=star(r)*r
    return m,defect,dict(common_RHS_count=rec+ns, CG_iterations=iterations,
                         shared_Frobenius_defect_upper=high(frobenius_squared(r).sqrt()),**times(start,cpu))


def one_row(p,obs,m, defect, gap):
    start,cpu=perf_counter(),process_time()
    ns=p.g.s
    yb=ball(obs.y)
    eye=acb_mat(ns,ns)
    for i in range(ns):
        eye[i,i]=1
    c=vertical(yb,eye)
    h0=star(yb)*yb+arb(p.mu)*(star(c)*m*c)
    dr=star(c)*defect*c
    error=(arb(p.mu)*sum((dr[i,i].real for i in range(ns)),arb(0))/arb(gap)).upper()
    vals,v=eigh(h(mid(h0)),check_finite=False,driver='evd')
    vals,v=vals[::-1],v[:,::-1]
    vb=ball(v)
    gram=star(vb)*vb
    hproj=star(vb)*h0*vb
    lower=[]
    for k in range(1,ns+1):
        if vals[k-1]<=0:
            break
        ff=acb_mat(k,k,[hproj[i,j] for i in range(k) for j in range(k)])
        gg=acb_mat(k,k,[gram[i,j] for i in range(k) for j in range(k)])
        num=gershgorin(ff,upper=False)-error
        if num<=0:
            break
        lower.append((num/gershgorin(gg,upper=True)).lower())
    gmin=gershgorin(gram,upper=False)
    if gmin<=0:
        raise ArithmeticError('Source mode basis condition is unresolved.')
    for i,d in enumerate(vals):
        hproj[i,i]-=arb(float(d))
    upper_precision=frobenius_squared(hproj).sqrt().upper()
    eigen_uppers=[max(arb(0),(arb(float(v0))+upper_precision).upper())/gmin for v0 in vals]
    d=low(sum(lower,arb(0)).sqrt())
    d=max(0.,d)
    upper=high(sum(eigen_uppers,arb(0)).sqrt())
    coherent_upper=high(eigen_uppers[0].sqrt())
    coherent_lower=max(0.,low(lower[0].sqrt())) if lower else 0.
    return dict(id=obs.id,object_id=obs.object_id,noise_pct=obs.noise_pct,cap=CAP,nu=0,
                radius=obs.padded_radius, distance_lower=d,distance_upper=upper,
                coherent_distance_lower=coherent_lower,coherent_distance_upper=coherent_upper,
                margin_lower=low(arb(d)-arb(obs.padded_radius)),
                strictly_excludes=d>obs.padded_radius,
                coherent_excludes=coherent_lower>obs.padded_radius,
                positive_source_lower=[low(e) for e in lower],
                source_matrix_eigen_uppers=[high(e) for e in eigen_uppers],
                shared_source_solution_correction_upper=high(error),**times(start,cpu))


def main():
    ctx.prec,ctx.threads=112,2
    start,cpu=perf_counter(),process_time()
    spec,s,b,data=full_ring_inputs()
    profiles=json.loads((SINGLE/'coherent_general_profiles.json').read_text())
    verified=json.loads((SINGLE/'coherent_general_verified.json').read_text())
    all_records=[]
    g=Geometry(spec,s,b,block_side=4)
    out=dict(scope='Frozen 4×4 profiles cross-applied with no tuning to six materials and five paired noise levels',
             native_dimension=g.n,cap=CAP,proof_source='../single/coherent_general_verified.json',
             coefficient_source='../single/coherent_general_profiles.json',profiles=all_records,
             inherited_common_proof_timing={k:verified[k] for k in ('wall_seconds','cpu_seconds')})
    for original_id,proof in verified['cases'].items():
        index=proof['profile_index']
        p=SavedProfile(g,profiles[index],CAP)
        gap=proof['proved_K_scaled_gap']
        if proof['cap']!=CAP or gap<=0:
            raise ArithmeticError('An eligible common positive native proof is needed.')
        cost=g.cost.copy()
        t0,c0=perf_counter(),process_time()
        m,defect,timing=reusable_compact(p)
        rows=[one_row(p,obs,m,defect,gap) for obs in data]
        all_records.append(dict(profile_index=index,inherited_owner_id=original_id,cap=CAP,
                               gap_lower=gap,specific_inherited_proof_wall_seconds=proof.get('wall_seconds'),
                               common_response=timing,rows=rows,
                               operations=dict(g.cost-cost),**times(t0,c0)))
        out.update(times(start,cpu))
        save_json(HERE/'fixed_blocks.75.json',out)
        print('fixed 4×4 profile',index,'passed',[(r['object_id'],r['noise_pct']) for r in rows
                                                   if r['strictly_excludes']],flush=True)
    highcap=[]
    for index in sorted(q['profile_index'] for q in verified['cases'].values()):
        for u in (1.,1.25,1.5):
            t0,c0=perf_counter(),process_time()
            p=SavedProfile(g,profiles[index],u)
            result=p.modes(10)
            status='numerically_below_source_budget'
            nlow=0
            negative_upper=None
            if p.nu>=g.s:
                negative_upper=p.ensure_negative(g.s)
                nlow=g.s
                status='ten_negative_directions_proved'
            highcap.append(dict(profile_index=index,cap=u,native_negative_lower_count=nlow,
                                numeric_negative_count=p.nu,negative_congruence_upper=negative_upper,
                                status=status,operations=dict(g.cost),numeric_mode_result=result,**times(t0,c0)))
            save_json(HERE/'fixed_blocks_higher_caps.json',dict(
                scope='Inertia probes for the fixed 4×4 weights across the higher caps only',
                cases=highcap,**times(start,cpu)))
            print('fixed 4×4 inertia',index,u,status,p.nu,flush=True)


if __name__=='__main__':
    with threadpool_limits(1):
        main()
