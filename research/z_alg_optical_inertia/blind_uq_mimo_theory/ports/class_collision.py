"""Constructive search for a shared 8x8 data ball at the widest material boxes."""
import argparse
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter, process_time
import subprocess

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.optimize import lsq_linear, minimize
from flint import acb, acb_mat, arb, ctx

from class_separation import endpoints, physical_arrays, row_format, spectral_from_gram, witness_matrices
from native_verify import ExplicitDots
from verify_witness import ball, frobenius, outward

HERE=Path(__file__).resolve().parent


class Model:
    def __init__(self,op):
        self.G,self.B,self.S=[op[k] for k in ("G","E0","S")]
        self.n=len(self.B)
        self.cost=dict(lu=0,forward_rhs=0,tangent_rhs=0,adjoint_rhs=0)

    def factor(self,m):
        A=np.eye(self.n)-m[:,None]*self.G
        lu=lu_factor(A,overwrite_a=True,check_finite=False)
        j=lu_solve(lu,m[:,None]*self.B,check_finite=False)
        self.cost["lu"]+=1
        self.cost["forward_rhs"]+=self.B.shape[1]
        return self.S@j,self.B+self.G@j,lu

    def jacobian(self,field,lu):
        L=lu_solve(lu,self.S.T,trans=1,check_finite=False).T
        self.cost["tangent_rhs"]+=self.S.shape[0]
        return (L[:,None,:]*field.T[None,:,:]).reshape(self.S.shape[0]*self.B.shape[1],-1)

    def pullback(self,field,lu,delta):
        v=lu_solve(lu,self.S.conj().T@delta,trans=2,check_finite=False)
        self.cost["adjoint_rhs"]+=delta.shape[1]
        return np.real(np.sum(v.conj()*field,axis=1))


class StopAtWitness(Exception):
    pass


def calibration_direction(h,eta,dots,direction):
    real=0.5*(h.imag+h.imag.T)
    ev,V=np.linalg.eigh(real)
    if np.max(np.abs(ev))<1e-18:
        return np.zeros_like(h)
    proposed=(V*np.sign(ev))@V.T
    proposed=0.5*(proposed+proposed.T)
    bound=float(np.nextafter(float(spectral_from_gram(proposed,dots,tight=True)),np.inf))
    change=1j*direction*(eta*(1-1e-12)/max(1.,bound))*proposed
    assert np.array_equal(change,change.T) and np.all(change.real==0)
    assert spectral_from_gram(change,dots,tight=True)<arb(".005")
    return change


def enclosed_port_model(point,original,D):
    base=ball(np.eye(D.shape[0]))
    a=ball(original["inverse_empty"])
    ea=original["error_inverse_empty"].upper()
    ru=point["full"]["eF"].upper()
    n=len(D)
    def allow_component(x,error):
        return acb(arb(float(x.real),error),arb(float(x.imag),error))
    a=acb_mat(n,n,[allow_component(z,ea) for z in original["inverse_empty"].ravel()])
    diff=acb_mat(n,n,[allow_component(z,ru) for z in point["full"]["F"].ravel()])
    b=a-diff/2
    Q=ball(D)
    dd=2*((base+a*Q).solve(a)-(base+b*Q).solve(b))
    row,col=original["train_index"],original["source_index"]
    scale=acb(arb(len(col)).sqrt())
    return acb_mat(len(row),len(col),[dd[i,j]/scale for i in row for j in col])


def two_ball_center(a,b,r):
    midpoint=(row_format(a)+row_format(b))/2
    y=ball(midpoint)
    errors=[frobenius(x-y)/r for x in (a,b)]
    return dict(shared_observation_outward=all(x<1 for x in errors),
                midpoint_distance_over_radius_interval=[outward(v) for v in errors]),midpoint


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--width",default=".05")
    args=p.parse_args()
    started=perf_counter(),process_time()
    ctx.prec,ctx.threads=128,1
    with np.load(HERE/"results/mechanism_witness.npz") as d:
        centers=np.array([d[k].ravel() for k in ("with_region","without_region")])
        radius=arb(12)*arb(float(d["reference_sigma"]))
    with np.load(HERE/"public/operators/n32_f400_a08.npz") as a:
        op=dict(a)
    m=Model(op)
    pert_bounds=list(zip(*(endpoints(c,float(args.width)-1e-12) for c in centers),strict=True))
    bounds=np.column_stack([np.concatenate(v) for v in pert_bounds])
    a,b=[m.factor(c) for c in centers]
    initial_diff=a[0]-b[0]
    scale=np.linalg.norm(initial_diff)**2
    linear=np.concatenate((m.jacobian(a[1],a[2]),-m.jacobian(b[1],b[2])),axis=1)
    basis=np.vstack((linear.real,linear.imag))
    d=initial_diff.ravel()
    fit=lsq_linear(basis,-np.r_[d.real,d.imag],bounds=bounds.T,method="trf",lsq_solver="lsmr",
                   max_iter=100,tol=1e-10)
    best, best_h=np.inf,None
    trace=[]
    target=.999*2*float(radius)
    def value(h):
        nonlocal best,best_h
        x,y=centers+h.reshape(centers.shape)
        A,B=[m.factor(c) for c in (x,y)]
        diff=A[0]-B[0]
        distance=np.linalg.norm(diff)
        trace.append(float(distance/(2*float(radius))))
        if distance<best:
            best=distance
            best_h=h.copy()
        if distance<target:
            raise StopAtWitness
        lhs=m.pullback(A[1],A[2],diff)
        rhs=m.pullback(B[1],B[2],diff)
        return .5*(distance**2)/scale,np.r_[lhs,-rhs]/scale
    try:
        outcome=minimize(value,fit.x,method="L-BFGS-B",jac=True,bounds=bounds,
                         options=dict(maxiter=100,ftol=2e-13,gtol=2e-8,maxls=20))
        stop=outcome.message
    except StopAtWitness:
        stop="candidate distance below .999 times 2r"
    materials=centers+best_h.reshape(centers.shape)
    maximum_shift=arb(0)
    for c,z in zip(centers.ravel(),materials.ravel(),strict=True):
        delta=abs(arb(float(c))-arb(float(z)))
        maximum_shift=max(maximum_shift,delta)
        if not (delta<=arb(args.width) and 0<=z<=4.5):
            raise ArithmeticError("Point outside frozen prescribed class")
    verification=perf_counter()
    with TemporaryDirectory() as td:
        lib=Path(td)/"verify.dylib"
        subprocess.run(["clang","-O3","-fno-fast-math","-ffp-contract=off","-dynamiclib",
                        "-o",str(lib),str(HERE/"certify_dot.c")],check=True)
        dots=ExplicitDots(lib,4)
        orig=physical_arrays(op,dots)
        a,b=[witness_matrices(orig,z,dots) for z in materials]
        no_error=np.zeros(orig["inverse_empty"].shape,complex)
        left,right=[enclosed_port_model(z,orig,no_error) for z in (a,b)]
        proof,mid=two_ball_center(left,right,radius)
        port_gap=row_format(left-right)
        U=port_gap/max(np.linalg.norm(port_gap),1e-300)
        C=np.zeros_like(no_error)
        C[np.ix_(orig["train_index"],orig["source_index"])]=U
        K_a=orig["inverse_empty"]
        sensitivity=[]
        for z in (a,b):
            P=z["full"]["F"]
            K_b=K_a-P/2
            sensitivity.append(P@C.conj().T@K_b+K_a@C.conj().T@P)
        D=calibration_direction(sensitivity[0]-sensitivity[1],.005,dots,-1.)
        Da=calibration_direction(sensitivity[0],.005,dots,-1.)
        Db=calibration_direction(sensitivity[1],.005,dots,1.)
        modes=dict(same_error=(D,D),independent_errors=(Da,Db))
        saved=dict(materials=materials,midpoint_nominal=mid,radius=float(radius),width=args.width)
        for mode,errors in modes.items():
            outputs=[enclosed_port_model(z,orig,D0) for z,D0 in zip((a,b),errors,strict=True)]
            calibrated,mid=two_ball_center(*outputs,radius)
            calibrated["normalized_error_norm_upper"]=[
                     outward(spectral_from_gram(D0,dots,tight=True))[1] for D0 in errors]
            proof[mode]=calibrated
            saved[mode+"_midpoint"]=mid
            saved[mode+"_D"]=np.asarray(errors)
        proof["costs"]=[a["raw"],b["raw"]]
        proof["seconds_wall"]=perf_counter()-verification
        np.savez_compressed(HERE/"results"/"mechanism_class_collision.npz",**saved)
    result=dict(width=args.width, maximum_shift_upper=outward(maximum_shift)[1], stop=str(stop),
                nominal_point_separation_over_two_r=best/(2*float(radius)),
                relaxed_affine_point_separation_over_two_r=np.linalg.norm(d+linear@fit.x)/(2*float(radius)),
                proof=proof, array_count=8, band_MHz=400, available_W=1, costs=m.cost, trace=trace,
                seconds_wall=perf_counter()-started[0],seconds_CPU=process_time()-started[1])
    (HERE/"results"/"mechanism_class_collision.json").write_text(json.dumps(result,indent=2))
    print({k:v for k,v in result.items() if k!="trace"},flush=True)


if __name__=="__main__":
    main()
