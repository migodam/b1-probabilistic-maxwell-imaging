"""Validate the turning of the material branch on the line through the two recorded data.

Use the sixth material as coordinate. Enclose the other five and signed data
offset by a parametric Krawczyk inclusion. Evaluate both implicit derivatives
with the same bounded physical Green matrix. Neighbouring branches are glued
by verifying an independently enclosed common root at their shared face.
"""
from pathlib import Path
import json
import time
import numpy as np
from scipy.optimize import root
from flint import arb, arb_mat, acb, acb_mat, ctx
from point_fibers import exact, from_bounds, bounds
from certify_fold import (response, diagonal, identity, pack_rows, pulse_operators,
                          full_matrix, sigma_min_enclosure)
from host_model import reduced_forward, pack
from noisy_cover import path_box

OUT = Path(__file__).resolve().parents[1]/'results'


def reals(col):
    return [col[i,0] for i in range(col.nrows())]


def column(items):
    return arb_mat([[x] for x in items])


def realify(row):
    return [part for v in row for part in (v.real,v.imag)]


def center_matrix(A):
    return arb_mat([[v.mid() for v in row] for row in A.tolist()])


def normal_graph():
    ctx.prec=144
    data=np.load(OUT/'host_operators.npz')
    box=json.loads((OUT/'host_operators_interval.json').read_text())
    hh={k:from_bounds(box[k]) for k in ('H','B','R','y0')}
    H,B,R=(acb_mat(hh[k]) for k in ('H','B','R'))
    y0=acb_mat([hh['y0']])
    # Exact line through the two saved observations. Midpoint and unit of
    # parameter are mathematical values of their binary64 inputs.
    plus=[exact(x) for x in data['y_plus']]
    minus=[exact(x) for x in data['y_minus']]
    tau=arb(float(data['tau']))
    yc=column(realify([(a+b)/2 for a,b in zip(plus,minus)]))
    n=column(realify([(a-b)/(2*tau) for a,b in zip(plus,minus)]))
    nf=np.array([float(v.mid()) for v in reals(n)])
    yf=np.array([float(v.mid()) for v in reals(yc)])
    zc=data['zc'].copy()
    z6=arb(float(zc[-1]))
    Afloat=(data['H'],data['B'],data['R'],data['y0'])
    I6=arb_mat([[int(i==j) for j in range(6)] for i in range(6)])

    def full_block(J):
        return arb_mat([[J[i,j] for j in range(5)]+[-n[i,0]] for i in range(6)])

    def fields(z):
        # Shift the resolvent to the center before widening. Direct interval
        # elimination on I-D_z H loses cancellations in its six pivots.
        m=[q.mid() for q in z]
        delta=diagonal([a-b for a,b in zip(z,m)])
        f0,_,C,E0,V0=response(m,H,B,R,y0)
        K0=H*C
        W=(identity(6)-K0*delta).inv()
        E,K=W*E0,W*K0
        V=V0+V0*delta*K
        J=pack_rows(acb_mat([[V[0,j]*E[j,s] for j in range(6)] for s in range(3)]))
        return f0+V0*delta*E,J,K,E,V

    def centered_solve(A,b):
        Y=center_matrix(center_matrix(A).inv())
        seed=center_matrix(Y*center_matrix(b))
        return seed+(Y*A).solve(Y*(b-A*seed))

    def fit(t):
        def fun(xx):
            z=np.r_[xx[:5],zc[-1]+t]
            yy,J,_,_=reduced_forward(z,*Afloat,derivatives=True)
            return pack(yy)-yf-nf*xx[-1],np.column_stack((J[:,:5],-nf))
        u=np.r_[zc[:5]+data['null'][:5]/data['null'][-1]*t,0.]
        sol=root(fun,u,jac=True,tol=1e-11)
        xx=sol.x
        f,jac=fun(xx)
        xx-=np.linalg.solve(jac,f)
        z=np.r_[xx[:5],zc[-1]+t]
        _,J,_,_=reduced_forward(z,*Afloat,derivatives=True)
        slope=np.linalg.solve(np.column_stack((J[:,:5],-nf)),-J[:,-1])
        return [arb(float(x)) for x in xx], [arb(float(x)) for x in slope]

    def details(tmid,h,*,point=False):
        x,dx=fit(tmid)
        t0=arb(tmid)
        width=arb(0,h)
        u=column(x)
        tangent=column(dx)
        center=x[:5]+[z6+t0]
        y,J,_,_,_=fields(center)
        F=column(realify(y.tolist()[0]))-yc-n*x[-1]
        Z=[center[i]+dx[i]*width for i in range(5)]+[z6+t0+width]
        # Bound residual of affine predictor by a mean-value enclosure.
        _,Jp,_,_,_=fields(Z)
        residual=F+(full_block(Jp)*tangent+column([Jp[i,5] for i in range(6)]))*width
        Y=center_matrix(full_block(J).inv())
        corrected=Y*residual
        base_r=[max(5e-13 if i<5 else 5e-16, float(v.abs_upper())*1.5)
                for i,v in enumerate(reals(corrected))]
        for factor in (1,2,4,8,16):
            radii=[arb(0,v*factor) for v in base_r]
            Zc=[Z[i]+radii[i] for i in range(5)]+[Z[5]]
            _,JJ,K,E,rec=fields(Zc)
            S=full_block(JJ)
            M=I6-Y*S
            kr=-corrected+M*column(radii)
            if all(a.contains_interior(b) for a,b in zip(radii,reals(kr))):
                ratio=max(sum(M[i,j].abs_upper()*radii[j].rad()/radii[i].rad()
                              for j in range(6)) for i in range(6))
                if ratio<1:
                    break
        else:
            raise ArithmeticError(f'no parameter enclosure t={tmid} h={h}')
        if point:
            return [p+q for p,q in zip(x,reals(kr))]
        first=centered_solve(S,column([-JJ[i,5] for i in range(6)]))
        dz=reals(first)[:5]+[arb(1)]
        D=diagonal(dz)
        second=2*rec*D*K*D*E
        acc=centered_solve(S,column([-q for q in realify(second.tolist()[0])]))
        return dict(t=tmid,h=h,x=x,dx=dx,rad=radii,der=first[5,0],
                    d2=acc[5,0],ratio=ratio,mid_value=x[-1],
                    values=u[5,0]+tangent[5,0]*width+radii[5])

    def cover(a,b,rows):
        mid=(a+b)/2
        try:
            slab=details(mid,(b-a)/2)
            if not slab['d2']<0:
                raise ArithmeticError('curvature sign inconclusive')
        except (ZeroDivisionError,ArithmeticError,ValueError):
            if b-a<2**-19:
                raise
            cover(a,mid,rows)
            cover(mid,b,rows)
            return
        rows.append(slab)

    rows=[]
    extent=.625
    cover(-extent,extent,rows)
    nodes=[q['t']-q['h'] for q in rows]+[extent]
    point_roots=[]
    for i,t in enumerate(nodes):
        p=details(t,0,point=True)
        point_roots.append(p)
        for j in (i-1,i):
            if j<0 or j>=len(rows):
                continue
            q=rows[j]
            tube=[m+s*(t-q['t'])+r for m,s,r in zip(q['x'],q['dx'],q['rad'])]
            if not all(a.contains_interior(b) for a,b in zip(tube,p)):
                raise ArithmeticError(f'glue at {t} failed')
    def point_slope(p,t):
        _,J,_,_,_=fields(p[:5]+[z6+arb(t)])
        return centered_solve(full_block(J),column([-J[i,5] for i in range(6)]))[5,0]
    end_slopes=[point_slope(point_roots[i],nodes[i]) for i in (0,-1)]
    assert end_slopes[0]>0 and end_slopes[1]<0
    assert all(point_roots[i][-1]<-tau for i in [0,-1])
    # A point on the negative-curvature curve and its tangent bound the single
    # maximum: g(0) <= max g <= g(0)+g'(0)^2/(2 * uniform |g''| lower).
    p0=details(0.,0,point=True)
    z=p0[:5]+[z6]
    _,J,_,_,_=response(z,H,B,R,y0)
    dz0=full_block(J).solve(column([-J[i,5] for i in range(6)]))
    curvature_margin=min((-r['d2']).lower() for r in rows)
    top=p0[-1].union(p0[-1]+dz0[5,0].abs_upper()**2/(2*curvature_margin))
    assert -tau<top and top<tau
    turning_radius=2.**-31
    assert dz0[5,0].abs_upper()/curvature_margin<arb(turning_radius)
    turn=details(0.,turning_radius)
    tbox=arb(0,turning_radius)
    zz=[x+s*tbox+r for x,s,r in zip(turn['x'][:5],turn['dx'][:5],turn['rad'][:5])]+[z6+tbox]
    G,_,_=pulse_operators()
    alm=full_matrix(zz,G)
    Adata=np.eye(len(data['m0']))-np.r_[zc,data['m0'][6:]][:,None]*data['G']
    vec=np.linalg.svd(Adata)[2][-1].conjugate()
    floor,upper,_=sigma_min_enclosure(alm,'0.26124',vec)
    receiver=fields(zz)[4]
    line_fold={'material_box':[bounds(q) for q in zz],
               'full_state_sigma_min':[bounds(floor)[0],bounds(upper)[1]],
               'strict_lower_decimal':'0.26124',
               'receiver_node_lower':min(bounds(q.abs_lower())[0] for q in receiver.entries())}
    # Every line datum between the endpoints lies in the same certified
    # noise ball, so each of its physical roots lies in a returned outer cell.
    # At fixed z6, eliminate p by comparing five transverse data components.
    # Prove uniqueness between each candidate and the already constructed curve.
    record=json.loads((OUT/'noisy_cover.json').read_text())
    order=record['axis_order_zero_based']
    cells=np.load(OUT/'noisy_cover_cells.npz')
    paths=np.r_[cells['physical_direct_centered_inside'],cells['physical_direct_centered_undecided']]
    projection=arb_mat([[int(i==j)-int(j==0)*n[i,0]/n[0,0] for j in range(6)] for i in range(1,6)])
    checked=0
    failures=[]
    for path in paths:
        loc,den=path_box(int(path),order)
        a=[.5+29.5*i/k for i,k in zip(loc,den)]
        b=[.5+29.5*(i+1)/k for i,k in zip(loc,den)]
        # All material bounds, z6 and the slab endpoints are binary dyadics.
        # Subtraction is exact here by Sterbenz (positive, within factor two).
        lo6,hi6=a[5]-zc[5],b[5]-zc[5]
        if lo6 < -extent or hi6 > extent:
            failures.append(('outside section',int(path)))
            continue
        for q in rows:
            lo=max(lo6,q['t']-q['h'])
            hi=min(hi6,q['t']+q['h'])
            if lo>hi:
                continue
            t=arb(lo).union(arb(hi))
            delta=t-q['t']
            crest=[x+s*delta+r for x,s,r in zip(q['x'][:5],q['dx'][:5],q['rad'][:5])]
            Z=[v.union(arb(l)).union(arb(u)) for l,u,v in zip(a,b,crest)]+[z6+t]
            _,J,_,_,_=fields(Z)
            transverse=projection*arb_mat([[J[i,j] for j in range(5)] for i in range(6)])
            # P encloses a constant nonsingular midpoint preconditioner.
            # Enclosing ALL matrices in the interval also encloses the segment mean.
            P=center_matrix(transverse).inv()
            try:
                inv=(P*transverse).inv()
                if not all(v.is_finite() for v in inv.entries()):
                    raise ValueError('nonfinite')
            except (ZeroDivisionError,ValueError):
                failures.append((int(path),lo,hi))
            checked+=1
    global_line={'outer_method':'physical_direct_centered','outer_cells':len(paths),
                 'cell_slab_intersections':checked,'unproved_intersections':failures,
                 'entire_box_covered':not failures}
    assert all(v>arb(.5) and v<arb(30) for q in rows for v in
               [x+s*arb(0,q['h'])+r for x,s,r in zip(q['x'][:5],q['dx'][:5],q['rad'][:5])]+
               [z6+q['t']+arb(0,q['h'])])
    return dict(material_coordinate='z6-z6_reference', interval=[-extent,extent],
                data_line='midpoint(y+,y-) + p*(y+ - y-)/(2*tau), exact recorded endpoints',
                glued_slabs=len(rows),critical_offset_bounds=bounds(top),
                pprime_zero=bounds(dz0[5,0]),global_second_derivative_bounds=[
                min(bounds(q['d2'])[0] for q in rows),max(bounds(q['d2'])[1] for q in rows)],
                max_contraction=max(float(q['ratio'].upper()) for q in rows),
                left_and_right_endpoint_values=[bounds(point_roots[i][-1]) for i in [0,-1]],
                endpoint_slopes=[bounds(q) for q in end_slopes],
                curve=[dict(t=t,p=bounds(p[-1]),z=[bounds(v) for v in p[:5]+[z6+arb(t)]])
                       for t,p in zip(nodes,point_roots)],
                slabs=[dict(t=q['t'],halfwidth=q['h'],slope=bounds(q['der']),
                            curvature=bounds(q['d2']),value=bounds(q['values'])) for q in rows],
                global_line=global_line,line_fold=line_fold,
                domain='normal observation segment only; global iff entire_box_covered; the two endpoints are also independently counted by the quotient')

if __name__=='__main__':
    start=time.perf_counter()
    result=normal_graph()
    result['wall_seconds']=time.perf_counter()-start
    (OUT/'normal_section.json').write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items() if k not in ('curve','slabs','global_line')})
    g=result['global_line']
    print('Global intersections',g['cell_slab_intersections'],'unproved',len(g['unproved_intersections']))
