"""Independent one-sided readout from the stored full-wave operator and fields."""
import argparse
import json

import numpy as np
from flint import acb, acb_mat, arb, ctx
from scipy.optimize import minimize_scalar
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

from beam_common import (HERE, OLD, atomic_write, clocks, elapsed, forward, geometry, herm, save)
from verify import ProofGrid, ball, diagonal_times, midpoint, lo, hi, tr
from native_residual import ResidualBoundBox


def norm_upper(a):
    return sum((abs(z).upper()**2 for z in a.entries()),arb(0)).sqrt().upper()


def eye(n):
    return acb_mat(n,n,[int(i==j) for i in range(n) for j in range(n)])


def column(a,j):
    return acb_mat(a.nrows(),1,[a[i,j] for i in range(a.nrows())])


def joined(*mats):
    return acb_mat(mats[0].nrows(),sum(a.ncols() for a in mats),
                   [z for rows in zip(*(x.tolist() for x in mats),strict=True) for row in rows for z in row])


def endpoints(a):
    return [lo(a),hi(a)]


class ExactBox:
    def __init__(self, exact, numeric, center, width):
        start=clocks()
        g,s,b=exact
        self.g,self.s,self.b=g,s,b
        self.n= g.nrows()
        self.p=s.nrows()
        self.m=b.ncols()
        self.center=[arb(float(x)) for x in center]
        self.width=[arb(float(x)) for x in width]
        self.center_np=center
        self.width_np=width
        self.A=eye(self.n)-diagonal_times(self.center,g)
        # A successful Arb solve independently encloses the exact inverse and
        # verifies joint uniqueness of the center state on the stored operator.
        self.U=self.A.solve(eye(self.n))
        self.D=s*self.U
        self.T=g*self.U
        self.Tw=tr(diagonal_times(self.width,tr(self.T)))
        self.wT=diagonal_times(self.width,self.T)
        self.J=self.U*diagonal_times(self.center,b)
        self.Y=s*self.J
        self.E=b+g*self.J
        self.eta=norm_upper(self.Tw)
        self.eta_reverse=norm_upper(self.wT)
        self.Ubound=norm_upper(self.U)
        self.Sbound=norm_upper(s)
        if not self.eta<1 or not self.eta_reverse<1:
            raise ArithmeticError('Contraction needed for entire-box secant or witness failed.')
        self.time=elapsed(start)

    def support(self,a,w):
        response=tr(w)*self.Y*a
        direction_weights=tr(diagonal_times(self.width,tr(tr(w)*self.D)))
        reach=norm_upper(direction_weights)*norm_upper(self.E*a)/(1-self.eta)
        center=response[0,0].real
        return (center-reach).lower(),(center+reach).upper()

    def physical_candidate(self, numeric, material, vectors, data, radii, *, exact_vectors=None):
        # The finite candidate current is recomputed, then verified using the
        # exact operator, not a fitted residual in the response enclosure.
        if not all(abs(arb(float(v))-c)<r for v,c,r in zip(material,self.center,self.width,strict=True)):
            raise ValueError('Witness outside the material box.')
        G,S,B=numeric
        c=[arb(float(v)) for v in material]
        b=self.b*(ball(vectors) if exact_vectors is None else exact_vectors)
        y,_,_,j,_,_=forward(G,S,B@vectors,material,auxiliary=True)
        j=ball(j)
        defect=diagonal_times(c,b+self.g*j)-j
        pred=self.s*j
        target=ball(data)
        errors=[]
        pcs=[]
        clean=[]
        for k,r in enumerate(radii):
            eps=self.Sbound*self.Ubound*norm_upper(column(defect,k))/(1-self.eta_reverse)
            err=norm_upper(column(pred-target,k))+eps
            errors.append(hi(err/arb(float(r))))
            pcs.append(hi(eps))
            clean.append(hi(norm_upper(column(defect,k))))
        return errors,pcs,clean,midpoint(pred)

    def optical_base(self, mu):
        weight=acb(float(mu))
        we=diagonal_times(self.width,self.E)
        pa=eye(self.n)-tr(self.wT)*self.wT
        A0=weight*pa+tr(self.D)*self.D
        right_opt=-weight*tr(self.wT)*we
        response=joined(self.Y,-eye(self.p))
        L=joined(right_opt,acb_mat(self.n,self.p))+tr(self.D)*response
        C=acb_mat(self.m+self.p,self.m+self.p)
        little=-weight*tr(we)*we
        for i in range(self.m):
            for j in range(self.m):
                C[i,j]=little[i,j]
        C=C+tr(response)*response
        return A0,L,C

    def optical_q(self, form, *, base):
        """Rank-two exact optical change. Positive base follows from eta < 1."""
        mu=form['mu']
        loc=form['cell']
        ratio=form['ratio']
        A0,L,C=base
        ei=acb_mat(self.n,1,[int(i==loc) for i in range(self.n)])
        trow=acb_mat(1,self.n,[self.T[loc,i] for i in range(self.n)])
        srow=acb_mat(1,self.m+self.p,[self.E[loc,i] for i in range(self.m)]+[0]*self.p)
        delta=-acb(0,1)*acb(arb(float(mu))*arb(float(ratio))/2)
        U=joined(ei,tr(trow))
        H=acb_mat(2,2,[0,delta,delta.conjugate(),0])
        border=L+delta*(ei*srow)
        solved=A0.solve(joined(border,U))
        split=border.ncols()
        right=acb_mat(self.n,split,[solved[i,j] for i in range(self.n) for j in range(split)])
        u=acb_mat(self.n,2,[solved[i,split+j] for i in range(self.n) for j in range(2)])
        small=H.inv()+tr(U)*u
        det=(small[0,0]*small[1,1]-small[0,1]*small[1,0]).real
        diag=small[0,0].real
        # Haynsworth: n_-(A0+U H U*)=n_+(H^-1+U*A0^-1U)-n_+(H^-1).
        if det>0 and diag>0:
            nu=1
        elif det<0:
            nu=0
        else:
            return None,dict(state_negative=None,woodbury_determinant=endpoints(det))
        stationary=right-u*small.solve(tr(U)*right)
        q=C-tr(border)*stationary
        return q,dict(state_negative=nu,woodbury_determinant=endpoints(det),cell=loc,mu=mu,ratio=ratio)


def positive_lower(a):
    nominal=herm(midpoint(a))
    try:
        f=ball(np.linalg.cholesky(nominal))
    except np.linalg.LinAlgError:
        return None
    err=norm_upper(a-f*tr(f))
    floor=1/norm_upper(f.inv())**2-err
    return floor.lower() if floor>0 else None


def full_dictionary_impossibility(models,power,sigma,energy,old_witnesses):
    """Necessary fluence for worst-case discrimination, even with coherent sources."""
    contrast=models[1].Y-models[0].Y
    perturbation=sum((getattr(x,'y_error',arb(0)) for x in models),arb(0))
    gram_floor=positive_lower(power)
    if gram_floor is None:
        return dict(verified=False)
    midpoint_delta=midpoint(contrast)
    nominal=float(eigh(herm(midpoint_delta.conj().T@midpoint_delta),herm(midpoint(power)),
                       eigvals_only=True,check_finite=False)[-1])
    Q=tr(contrast)*contrast
    for d in (1.00001,1.0001,1.001,1.01,1.1,2.):
        q=arb(nominal)*arb(d)
        if positive_lower(acb(q)*power-Q) is None:
            continue
        maximum_gain=q.sqrt()+perturbation/gram_floor.sqrt()
        observed=energy.sqrt()*maximum_gain
        return dict(verified=True,
                    full_dictionary_includes_old_beam=True,
                    old_both_centers_fit=old_witnesses,
                    response_operator_norm_per_sqrt_fluence_upper=hi(maximum_gain),
                    maximum_center_difference_over_new_radius_upper=hi(observed/arb(sigma)),
                    worst_case_separation_impossible=bool(old_witnesses and observed<2*arb(sigma)),
                    necessary_fluence_lower=lo((2*arb(sigma)/maximum_gain)**2))
    return dict(verified=False)


def trust_direction(q,m,action,measured,r):
    receive=acb_mat(q.nrows()-m,q.ncols()-m,[q[i,j] for i in range(m,q.nrows()) for j in range(m,q.ncols())])
    both=acb_mat(m+measured.nrows(),1,action.entries()+measured.entries())
    val=(tr(both)*q*both)[0,0].real
    ell=acb_mat(receive.nrows(),1,(q*both).entries()[m:])
    r_np=herm(midpoint(receive))
    vals,vec=eigh(r_np)
    powers=np.abs(vec.conj().T@midpoint(ell).ravel())**2
    minimum=max(0.,-float(vals[0]))+1.e-10
    def objective(eta):
        return float(val.mid())-eta*r*r-np.sum(powers/(vals+eta))
    high=max(minimum*1.01,minimum+1.e-5,np.linalg.norm(midpoint(ell))/r*2+minimum)
    proposal=minimize_scalar(lambda ep:-objective(ep),bounds=(minimum,high),method='bounded',
                             options={'xatol':1.e-12}).x
    lower=None
    for eta in (proposal,minimum,high):
        if eta<0 or eta+vals[0]<=0:
            continue
        try:
            f=ball(np.linalg.cholesky(r_np+eta*np.eye(len(vals))))
        except np.linalg.LinAlgError:
            continue
        err=norm_upper(receive+acb(float(eta))*eye(len(vals))-f*tr(f))
        # For f*f*: bound sigma_min from the Arb enclosure of f^-1.
        df=1/norm_upper(f.inv())**2
        if not df>err:
            continue
        penalty=(tr(ell)*(receive+acb(float(eta))*eye(len(vals))).solve(ell))[0,0].real
        bound=val-arb(float(eta))*arb(float(r))**2-penalty
        candidate=(lo(bound),float(eta),lo(df-err))
        lower=candidate if lower is None or candidate[0]>lower[0] else lower
    return lower


def conditional_q(q, old, radius, m):
    left=q*old
    pivot=(tr(old)*left)[0,0].real-arb(float(radius))**2
    if pivot>0:
        return q-acb(1/pivot)*(left*tr(left)),pivot,1
    if pivot<0:
        return q-acb(1/pivot)*(left*tr(left)),pivot,0
    return None,pivot,None


def verify_case(exact, bank, pbank, record, proposal, xold, xnew, truth, certify_optical, green_error, power):
    idx=str(record['index'])
    material,widths,fg,bg=[bank[idx+'__'+j] for j in ('centers','widths','fg','bg')]
    vec=pbank[idx+'__actions']
    directions=pbank[idx+'__directions']
    old_v=np.eye(vec.shape[0],dtype=complex)[:,[OLD]]
    action_all=np.column_stack((old_v,vec))
    all_data=np.column_stack((xold,xnew))
    rad=proposal['prescribed_noise_balls']
    sigma=rad['new']
    radii=[rad['old']]+[sigma]*vec.shape[1]
    yold=ball(xold)
    old_pair=acb_mat(vec.shape[0]+len(xold),1,ball(old_v).entries()+yold.entries())
    target=arb(proposal['same_fluence'])
    V_raw=ball(vec)
    raw_fluence=tr(V_raw)*power*V_raw
    scales=[(target/raw_fluence[k,k].real).sqrt() for k in range(V_raw.ncols())]
    V=joined(*(acb(scale)*column(V_raw,k) for k,scale in enumerate(scales)))
    action_all_exact=joined(ball(old_v),V)
    class_records=[]
    storage={}
    secant_ends=[]
    optical_bounds=[]
    verified_models=[]
    midpoint_witnesses=[]
    for class_idx, center in enumerate(material):
        ncase=clocks()
        model=ExactBox(exact[1],exact[0],center,widths) if green_error is None else ResidualBoundBox(exact[1],exact[0],center,widths,green_error)
        verified_models.append(model)
        witness_list=[center]
        alt=np.where(fg,-.46,.56)
        witness_list.extend((center+widths*alt,center-widths*alt))
        wrs=[]
        predicted=[]
        for v in witness_list:
            r,eps,defect,clean=ExactBox.physical_candidate(model,exact[0],v,action_all,all_data,radii,exact_vectors=action_all_exact)
            wrs.append(dict(residual_over_radius_upper=r,verified_old_compatible=bool(r[0]<=1.),
                            receiver_forward_error_upper=eps,state_defect_upper=defect,
                            material_rmse_from_center=float(np.linalg.norm(v-center)/np.sqrt(len(v)))))
            predicted.append(clean)
        storage[f'{class_idx}__witness_materials']=np.stack(witness_list)
        storage[f'{class_idx}__witness_clean_midpoints']=np.stack(predicted)
        midpoint_witnesses.append((predicted[0],wrs[0]['receiver_forward_error_upper']))
        row_boxes=[]
        row_opt=[]
        if certify_optical:
            # Selection uses only the design bank of old-conditioned forms.
            frozen=proposal['designs']['source_shaped']['optical_support_candidate']
            cross=frozen[class_idx] if frozen is not None else None
            if cross is None:
                forms=[]
            else:
                forms=[cross]
        else:
            forms=[]
        for choice in forms:
            tick=clocks()
            q,cert=model.optical_q(choice,base=model.optical_base(choice['mu']))
            if q is None:
                row_opt.append(dict(**cert,seconds=elapsed(tick)))
                continue
            qcond,piv,k=conditional_q(q,old_pair,rad['old'],vec.shape[0])
            cert.update(old_pivot_interval=endpoints(piv),old_positive=k,saturated=k==cert.get('state_negative'))
            statements=[]
            if cert['saturated']:
                for j,name in enumerate(proposal['policies']):
                    read=trust_direction(qcond,vec.shape[0],column(V,j),ball(xnew[:,j,None]),sigma)
                    statements.append(dict(policy=name, bound=None if read is None else read[0],
                                           multiplier=None if read is None else read[1],
                                           verified_shift_gap=None if read is None else read[2],
                                           excludes=bool(read is not None and read[0]>0)))
                storage[f'{class_idx}__q_conditional_midpoint']=midpoint(qcond)
            cert.update(outcomes=statements,**elapsed(tick))
            row_opt.append(cert)
        optical_bounds.append(row_opt)
        for j in range(vec.shape[1]):
            w=ball(directions[:,j,None])
            a=column(V,j)
            lower,upper=model.support(a,w)
            wr=norm_upper(w)*arb(float(sigma))
            p=(tr(w)*ball(xnew[:,j,None]))[0,0].real
            dist=max((lower-p).lower(),(p-upper).lower())
            row_boxes.append((lower,upper,wr,dist))
        secant_ends.append(row_boxes)
        class_records.append(dict(eta_upper=hi(model.eta),reverse_eta_upper=hi(model.eta_reverse),
                                  verified_old_witnesses=sum(x['verified_old_compatible'] for x in wrs),
                                  class_witnesses=wrs,setup=model.time,optical=row_opt,**elapsed(ncase)))
    outcomes={}
    for j,name in enumerate(proposal['policies']):
        statuses=[]
        for cls,res in enumerate(class_records):
            lower,upper,wr,dist=secant_ends[cls][j]
            raw0=dist-wr
            excl=bool(raw0>0)
            app=[w for w in res['class_witnesses'] if w['residual_over_radius_upper'][0]<=1
                 and w['residual_over_radius_upper'][j+1]<=1]
            opt=any(any(v['policy']==name and v['excludes'] for v in opt.get('outcomes',[]))
                    for opt in optical_bounds[cls])
            if app and (excl or opt):
                raise ArithmeticError('Independent physical witness contradicts exclusion.')
            statuses.append(dict(state='EXCLUDED' if excl or opt else 'WITNESSED' if app else 'UNRESOLVED',
                                 interval_exclusion_margin_over_new_radius_lower=lo(raw0/arb(float(sigma))),
                                 source_optical_exclusion=opt,number_joint_physical_witnesses=len(app)))
        A,B=secant_ends[0][j],secant_ends[1][j]
        gap=(B[0]-A[1]-A[2]-B[2])/arb(float(sigma))
        clean_a,err_a=midpoint_witnesses[0]
        clean_b,err_b=midpoint_witnesses[1]
        clean_a=ball(clean_a[:,j+1,None])
        clean_b=ball(clean_b[:,j+1,None])
        observed_center=ball(midpoint((clean_a+clean_b)*acb(arb(1)/2)))
        center_dist=norm_upper(clean_a-clean_b)+arb(err_a[j+1])+arb(err_b[j+1])
        both_radius=[hi((norm_upper(x-observed_center)+arb(e[j+1]))/arb(sigma)) for x,e in ((clean_a,err_a),(clean_b,err_b))]
        end_gap=(B[0]-A[1])
        outcomes[name]=dict(classes=statuses,
                            class_center_distance_over_noise_upper=hi(center_dist/arb(sigma)),
                            common_adversarial_midpoint=[complex(z).real for z in observed_center.entries()],
                            common_adversarial_midpoint_imag=[complex(z).imag for z in observed_center.entries()],
                            common_midpoint_residual_over_radius_upper=both_radius,
                            midpoint_guarantee_impossible=bool(max(both_radius)<1 and all(r['verified_old_witnesses']>=1 for r in class_records)),
                            fixed_ray_necessary_energy_as_fraction_of_current_lower=lo((2*arb(sigma)/center_dist)**2),
                            fixed_ray_sufficient_energy_as_fraction_of_current_upper=hi(((A[2]+B[2])/end_gap)**2) if end_gap>0 else None,verified_prospective_full_box_gap_over_new_radius_lower=lo(gap),
                            actually_exactly_one_witnessed_and_other_excluded=sorted(c['state'] for c in statuses)==['EXCLUDED','WITNESSED'],
                            ambiguity_lower_witness_statement='joint observed' if all(c['state']=='WITNESSED' for c in statuses) else None)
    # Score the generator only after all decisions are made. Its independent
    # forward enclosure also verifies the injected noise and all baselines.
    label=int(truth[idx+'__label'])
    gen=verified_models[label]
    truth_r,truth_err,_,_=ExactBox.physical_candidate(gen,exact[0],truth[idx+'__chi'],action_all,all_data,radii,exact_vectors=action_all_exact)
    if max(truth_r)>=1:
        raise ArithmeticError('Generator outside acquired contract.')
    for decision in outcomes.values():
        if decision['classes'][label]['state']=='EXCLUDED':
            raise ArithmeticError('Generator excluded by strict audit.')
    bound=full_dictionary_impossibility(verified_models,power,rad['new'],target,
                                        all(c['verified_old_witnesses']>=1 for c in class_records))
    return dict(classes=class_records,policies=outcomes,truth_label_offline=label,
                exact_common_target_fluence=proposal['same_fluence'],
                exact_fluence_normalization_enclosures=[endpoints(x) for x in scales],
                all_coherent_actions=bound,
                truth_residual_over_radius_upper=truth_r,truth_forward_error_upper=truth_err,
                ),storage


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('side',type=int,choices=(32,96))
    parser.add_argument('--tag',default='main')
    parser.add_argument('--no-optical',action='store_true')
    parser.add_argument('--fast-native',action='store_true')
    parser.add_argument('--folder',type=type(HERE))
    parser.add_argument('--output',default=None)
    args=parser.parse_args()
    ctx.prec=128
    ctx.threads=3
    folder=args.folder or HERE/f'{args.tag}_{args.side}'
    reg=json.loads((folder/'registration.json').read_text())
    prop=json.loads((folder/'proposal.json').read_text())
    bank=np.load(folder/'public_classes.npz')
    pbank=np.load(folder/'proposal_arrays.npz')
    old=np.load(folder/'old_readings.npz')
    new=np.load(folder/'future_readings.npz')
    truth=np.load(folder/'future_offline.npz')
    tick=clocks()
    geom=geometry(args.side)
    pgrid=ProofGrid(geom)
    green_err_per_entry=None
    if args.fast_native:
        if geom.spectrum is None or not args.no_optical:
            raise ValueError('--fast-native requires 96 and --no-optical')
        kernel=np.fft.ifft2(geom.spectrum).ravel()
        green_err_per_entry=max(abs(v-acb(float(n.real),float(n.imag))).upper() for v,n in zip(pgrid.kernel,kernel,strict=True))
    power=tr(ball(geom.B))*ball(geom.B)
    power=acb(1/power[0,0].real)*power
    init=elapsed(tick)
    rows={}
    arrays={}
    name=args.output or ('results_secant.json' if args.no_optical else 'results.json')
    path=folder/name
    for record in reg['cases']:
        i=str(record['index'])
        if i not in prop['rows']:
            continue
        row_start=clocks()
        cells=bank[i+'__cells']
        matched=pgrid.restricted(cells)
        coords=ball(pbank[i+'__actions'])
        ed=tr(coords)*power*coords
        energy_intervals=[endpoints(ed[j,j].real) for j in range(ed.nrows())]
        if not np.array_equal(new[i+'__V'],pbank[i+'__actions']):
            raise ArithmeticError('Fresh data action differs from frozen proposal.')
        green_error=None if green_err_per_entry is None else green_err_per_entry*len(cells)
        row,store=verify_case(matched,bank,pbank,record,prop['rows'][i],old[i+'__Y'],new[i+'__Y'],truth,
                              certify_optical=not args.no_optical,green_error=green_error,power=power)
        row.pop('wall_seconds',None)
        row.pop('cpu_seconds',None)
        row.update(incident_energy_intervals=energy_intervals,**elapsed(row_start))
        rows[i]=row
        arrays.update({i+'__'+k:v for k,v in store.items()})
        out=dict(side=args.side,method='128-bit Arb exact stored operator, entire real boxes and original pointwise Maxwell witnesses',
                 strict_optical=not args.no_optical,verified_float_native_inverse=args.fast_native,
                 proof_setup=init,rows=rows,**elapsed(tick))
        atomic_write(path,out)
        print(i,[c['verified_old_witnesses'] for c in row['classes']],
              [(name, [v['state'] for v in outcome['classes']]) for name,outcome in row['policies'].items()],
              row['wall_seconds'],flush=True)
    save(folder/('verification_secant_arrays.npz' if args.no_optical else 'verification_arrays.npz'),**arrays)


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
