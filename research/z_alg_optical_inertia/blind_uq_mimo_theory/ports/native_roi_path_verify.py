"""Independent activity-block replay of the 96² zero-to-peak material path."""
import argparse
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
from time import perf_counter, process_time

import numpy as np
from flint import acb, acb_mat, arb, ctx
from scipy.linalg import lu_factor, lu_solve
from scipy.ndimage import label

from native_verify import ExplicitDots, gamma, tiny
from verify_witness import as_acb, ball, exact_kernel, frobenius, outward

HERE=Path(__file__).resolve().parent
CASE="n96_c01"
PROFILE=HERE.parent/"results/gn_profiles"/CASE/"roi_02.npz"
PATH_RECEIPT=PROFILE.with_name("roi_02_inner.json")
EXTENDED=HERE.parent/"results/extended"/f"{CASE}.npz"


def times(start):
    return dict(wall_seconds=perf_counter()-start[0],CPU_seconds=process_time()-start[1])


def midpoint(a):
    return np.array([complex(v) for v in a.entries()]).reshape(a.nrows(),a.ncols())


def square(z):
    return z.real*z.real+z.imag*z.imag


def squared(a):
    return sum((square(z) for z in a.entries()),arb(0))


def block_lags(left,right,shape,pad):
    yl,xl=np.unravel_index(left,shape)
    yr,xr=np.unravel_index(right,shape)
    return ((yl[:,None]-yr[None,:])%pad[0])*pad[1]+((xl[:,None]-xr[None,:])%pad[1])


class RestrictedModel:
    """A matrix restriction only at exactly zero current pixels of a listed point."""
    def __init__(self,data):
        self.data=data
        start=perf_counter(),process_time()
        self.kernel=exact_kernel(data["spectrum"])
        self.kernel_center=np.asarray([complex(k) for k in self.kernel])
        self.kernel_error=max(abs(k-as_acb(z)).upper()
                              for k,z in zip(self.kernel,self.kernel_center,strict=True))
        self.kernel_cost=times(start)
        self.shape=tuple(data["grid_shape"])
        self.pad=tuple(data["pad"])

    def verify_system(self,m,ids,dots):
        start=perf_counter(),process_time()
        n=len(ids)
        x=m[ids]
        max_x=arb(float(np.max(x)))
        # Build bare and loaded factors in blocks so the 9k system remains
        # below peak RAM of an explicit full pair of coordinate matrices.
        G=np.empty((n,n),complex)
        bare_sq,rank_sq=arb(0),arb(0)
        rank_error=arb(0)
        step=512
        for i in range(0,n,step):
            rows=ids[i:i+step]
            lags=block_lags(rows,ids,self.shape,self.pad)
            bare=self.kernel_center[lags]
            rank,r_bound=dots.product(self.data["low_rank_left"][rows],self.data["low_rank_right"][:,ids])
            bare_sq+=dots.norm(bare)**2
            rank_sq+=dots.norm(rank)**2
            rank_error+=r_bound**2
            G[i:i+step]=bare+rank
        del lags,bare,rank
        g_norm=dots.norm(G)
        g_error=n*self.kernel_error+rank_error.sqrt()+gamma(5)*(bare_sq.sqrt()+rank_sq.sqrt())+8*n*tiny
        A=-x[:,None]*G
        del G
        A.flat[::n+1]+=1
        a_error=max_x*g_error+gamma(6)*(arb(n).sqrt()+max_x*g_norm)+8*n*tiny
        assembly=times(start)
        start=perf_counter(),process_time()
        lu=lu_factor(np.asfortranarray(A),overwrite_a=True,check_finite=False)
        factor_cost=times(start)
        solve_cost=dict(wall_seconds=0.,CPU_seconds=0.)
        product_cost=dict(wall_seconds=0.,CPU_seconds=0.)
        check_cols=(512 if n>7000 else n)
        inverse_sq,defect_sq=arb(0),arb(0)
        for i in range(0,n,check_cols):
            indices=np.arange(i,min(i+check_cols,n))
            k=len(indices)
            eye_cols=np.zeros((n,k),complex)
            eye_cols[indices,np.arange(k)]=1.
            start=perf_counter(),process_time()
            inv=lu_solve(lu,eye_cols,check_finite=False)
            inverse_sq+=dots.norm(inv)**2
            for key,value in times(start).items():
                solve_cost[key]+=value
            start=perf_counter(),process_time()
            ax,product_err=dots.product(A,inv)
            del inv
            block_defect=(dots.norm(eye_cols-ax)+gamma(5)*(arb(k).sqrt()+dots.norm(ax))
                          +8*n*tiny+product_err)
            defect_sq+=block_defect**2
            for key,value in times(start).items():
                product_cost[key]+=value
            del ax,eye_cols
        inv_norm=inverse_sq.sqrt()
        defect=defect_sq.sqrt()+a_error*inv_norm
        lu_cost={key:factor_cost[key]+solve_cost[key] for key in factor_cost}
        if not defect<1:
            raise ArithmeticError(f"Right inverse defect did not prove invertibility: {defect}")
        inverse=inv_norm/(1-defect)
        return A, lu, a_error, inverse, dict(active_cells=n,dense_LU=1,paid_inverse_columns=n,
            operator_matrix_error_upper=outward(a_error)[1],inverse_defect_upper=outward(defect)[1],inverse_columns_per_chunk=check_cols,
            inverse_2_upper_from_Fro=outward(inverse)[1],assembly=assembly,LU_and_inverse=lu_cost,
            verified_binary64_product=product_cost)

    def exact_block(self,rows,cols):
        lags=block_lags(rows,cols,self.shape,self.pad)
        bare=acb_mat(len(rows),len(cols),[self.kernel[int(k)] for k in lags.ravel()])
        return bare+ball(self.data["low_rank_left"][rows])*ball(self.data["low_rank_right"][:,cols])


def checked_solve(A,lu,eA,inv,rhs,erhs,dots):
    start=perf_counter(),process_time()
    j=lu_solve(lu,rhs,check_finite=False)
    reconstructed,product_error=dots.product(A,j)
    defect=(dots.norm(reconstructed-rhs)+gamma(5)*(dots.norm(reconstructed)+dots.norm(rhs))
            +eA*dots.norm(j)+erhs+product_error+8*rhs.size*tiny)
    return j,inv*defect,dict(physical_rhs=rhs.shape[1],residual_Fro_upper=outward(defect)[1],**times(start))


def scaled_rhs(m,b,dots):
    return m[:,None]*b,gamma(5)*arb(float(np.max(m)))*dots.norm(b)+8*b.size*tiny


def product_enclosure(A,eA,B,eB,dots):
    rounded,r=dots.product(A,B)
    error=eA*(dots.norm(B)+eB)+dots.norm(A)*eB+r
    return rounded,error


def enclose_matrix(values,error):
    e=error.upper()
    return acb_mat(values.shape[0],values.shape[1],
                    [acb(arb(float(v.real),e),arb(float(v.imag),e)) for v in values.ravel()])


def exact_small_center(a):
    values=midpoint(a)
    return values, frobenius(a-ball(values)).upper()


def bernstein_three(a,b,c,max_depth=4):
    """Refine the range of a quadratic expressed by three Bernstein coefficients."""
    pending=[(a,b,c)]
    count=1
    for _ in range(max_depth):
        if all(a>0 and b>0 and c>0 for a,b,c in pending):
            break
        if all(a<0 and b<0 and c<0 for a,b,c in pending):
            break
        more=[]
        for a,b,c in pending:
            d=(a+b)/2
            e=(b+c)/2
            f=(d+e)/2
            more.extend([(a,d,f),(f,e,c)])
        pending=more
        count*=2
    coefficients=[v for p in pending for v in p]
    return min(v.lower() for v in coefficients),max(v.upper() for v in coefficients),count


def replay_path(model,obs,dots):
    start=perf_counter(),process_time()
    p=json.loads(PATH_RECEIPT.read_text())["selected"]["pixel"]
    with np.load(PROFILE) as a:
        x=a["material"].ravel()
        roi=a["roi"].ravel()
    if x.shape!=(9216,) or roi.sum()!=624 or x[roi].any() or not roi[p] or np.any((x<0)|(x>4.5)):
        raise ValueError("The frozen profile must have exactly zero contrast on the full ROI")
    ids=np.flatnonzero(x)
    pos=np.array([p])
    A,lu,eA,inv,charges=model.verify_system(x,ids,dots)
    data=model.data
    B=data["E0"][ids]
    S=data["S"][:,ids]
    rr,err_rr=scaled_rhs(x[ids],B,dots)
    g,eg=exact_small_center(model.exact_block(ids,pos))
    v,err_v=scaled_rhs(x[ids],g,dots)
    J,eJ,base=checked_solve(A,lu,eA,inv,rr,err_rr,dots)
    W,eW,tangent=checked_solve(A,lu,eA,inv,v,err_v+arb(float(max(x)))*eg,dots)
    del A,lu
    Y0v,eY0=product_enclosure(S,arb(0),J,eJ,dots)
    sensor_v,e_sensor=product_enclosure(S,arb(0),W,eW,dots)
    gp,egp=exact_small_center(model.exact_block(pos,ids))
    field_v,e_field=product_enclosure(gp,egp,J,eJ,dots)
    h_v,e_h=product_enclosure(gp,egp,W,eW,dots)
    h=(model.exact_block(pos,pos)+enclose_matrix(h_v,e_h))[0,0]
    field=ball(data["E0"][pos])+enclose_matrix(field_v,e_field)
    sensor=ball(data["S"][:,pos])+enclose_matrix(sensor_v,e_sensor)
    variation=sensor*field
    baseline=enclose_matrix(Y0v,eY0)-ball(obs["Y"])
    r=arb(float(obs["radius"]))
    length=arb("4.5")
    delta=variation-acb(h)*baseline
    a=squared(baseline)-r*r
    b=(baseline.conjugate().transpose()*delta)
    inner=sum((b[i,i].real for i in range(b.ncols())),arb(0))
    center=a+length*(inner+r*r*h.real)
    end=squared(baseline+acb(length)*delta)-r*r*square(acb(1)-length*h)
    den=bernstein_three(arb(1),arb(1)-length*h.real,square(acb(1)-length*h))
    res=bernstein_three(a,center,end)
    if not den[0]>0 or not res[1]<0:
        raise ArithmeticError(f"No uniform one-cell path proof, denominator {den}, quadratic {res}")
    checkpoints=[]
    for t in (arb(0),length/2,length):
        D=baseline+acb(t/(1-t*h))*variation
        checkpoints.append(dict(contrast=float(t),residual_over_r=outward(frobenius(D)/r)))
    return dict(case_id=CASE,profile=str(PROFILE.relative_to(HERE.parent)),roi_cells=int(roi.sum()),
                changed_index_flat_zero_based=int(p),changed_index_row_col=list(map(int,np.unravel_index(p,model.shape))),
                state_dimension_global=9216,strict_all_segment_feasible=True,
                attainable_peak_range=[0.,4.5],noise_radius_sqrtW=outward(r),
                loading_scalar_H_real_interval=outward(h.real),loading_scalar_H_imag_interval=outward(h.imag),
                path_denominator_squared_min_lower=outward(den[0])[0],
                quadratic_over_radius_squared_max_upper=outward(res[1]/(r*r))[1],
                denominator_segments=den[2],quadratic_segments=res[2],
                checkpoints=checkpoints,represented_data="stored binary64 FFT spectrum, exact low-rank port loading and stored E0,S,Y",
                scope="One fixed whole-material path; exact zero currents eliminated for point replay only. Training data ball.",
                costs=dict(system=charges,base=base,increment_column=tangent,**times(start)))


def replay_extended(model,obs,dots):
    start=perf_counter(),process_time()
    with np.load(EXTENDED) as p:
        x=p["material"].ravel()
        patch=p["patch"]
        roi=p["roi"].ravel()
    if np.any((x<0)|(x>4.5)) or not np.all(roi[patch]) or np.min(x[patch])<2 or len(patch)!=156:
        raise ValueError("Frozen material is outside its box or violates the patch lower bound")
    with np.load(PROFILE) as saved:
        base=saved["material"]
        if not np.array_equal(saved["roi"].ravel(),roi):
            raise ValueError("The extended witness must replay the same reported ROI")
    mask=np.zeros(9216,dtype=bool)
    mask[patch]=True
    if np.any(base[patch]!=0) or label(mask.reshape(96,96))[1]!=1:
        raise ValueError("The certified comparison point must be zero on the connected patch")
    ids=np.flatnonzero(x)
    A,lu,eA,inv,charges=model.verify_system(x,ids,dots)
    rhs,erhs=scaled_rhs(x[ids],model.data["E0"][ids],dots)
    j,ej,state=checked_solve(A,lu,eA,inv,rhs,erhs,dots)
    del A,lu
    out={}
    for name,k,k_radius in (("train","Y","radius"),("hold","Y_hold","radius_hold")):
        obs_resp=obs[k]
        mid,err=product_enclosure(model.data["S" if name=="train" else "S_hold"][:,ids],arb(0),j,ej,dots)
        r=arb(float(obs[k_radius]))
        f=frobenius(enclose_matrix(mid,err)-ball(obs_resp))
        out[name]=dict(residual_over_r=outward(f/r),in_ball=bool(f<r),port_error_Fro_upper=outward(err)[1])
    if not out["train"]["in_ball"]:
        raise ArithmeticError("Material endpoint is not in the original training ball")
    area_cm2=156*(200/96)**2
    return dict(case_id=CASE,material=str(EXTENDED.relative_to(HERE.parent)),roi_cells=int(roi.sum()),
                changed_patch_cells=len(patch),patch_connected_four_neighbors=True,patch_area_cm2=area_cm2,
                min_patch_contrast=float(min(x[patch])),
                mean_patch_contrast=float(x[patch].mean()),mean_ROI_contrast=float(x[roi].mean()),
                delta_relative_to_zero_path_base="Exact zero on the whole patch for the profiled base material",
                costs=dict(system=charges,forward=state,**times(start)),**out)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--extended",action="store_true")
    p.add_argument("--threads",type=int,default=8)
    args=p.parse_args()
    ctx.prec,ctx.threads=128,1
    started=perf_counter(),process_time()
    obs=dict(np.load(HERE/"public"/"cases"/f"{CASE}.npz"))
    model=RestrictedModel(dict(np.load(HERE/"public/operators/n96_f400_a08.npz")))
    with TemporaryDirectory() as tmp:
        dy=Path(tmp)/"explicit_dots.dylib"
        subprocess.run(["clang","-O3","-fno-fast-math","-ffp-contract=off","-dynamiclib",
                        "-o",str(dy),str(HERE/"certify_dot.c")],check=True)
        dots=ExplicitDots(dy,args.threads)
        if args.extended:
            receipt=replay_extended(model,obs,dots)
            output="native_c01_extended_independent_arb.json"
        else:
            receipt=replay_path(model,obs,dots)
            output="native_c01_roi_path_independent_arb.json"
    report=dict(arb_precision_bits=ctx.prec,threads=args.threads,kernel_cost=model.kernel_cost,
                bare_kernel_error_upper=outward(model.kernel_error)[1],result=receipt,**times(started))
    (HERE/"results"/output).write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2),flush=True)


if __name__=="__main__":
    main()
