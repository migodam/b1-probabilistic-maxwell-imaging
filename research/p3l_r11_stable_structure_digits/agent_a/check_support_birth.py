"""One essential independent inactive-cell differential check on R10 physics."""
import os
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
import json
import numpy as np
from scipy.linalg import lu_factor,lu_solve
from threadpoolctl import threadpool_limits
from benchmark import HERE,Ledger,load_host,inputs_r10,SIGMA,write
from solver import Solver


def run():
    if (HERE/"SUPPORT_BIRTH_CHECK.json").exists():
        raise RuntimeError("preserve existing physical check")
    book=Ledger("essential_birth_check")
    try:
        st=load_host();identity,a,ds=next(inputs_r10(st["ROI"]));delta=complex(ds[0]);S=np.flatnonzero(a>0);I=np.flatnonzero(a==0)
        # A deterministic illuminated inactive cell, chosen without data or GT.
        j=int(I[np.argmax(np.linalg.norm(st["ER"][I],axis=1))]);da=np.zeros(256);da[j]=1.
        solver=Solver(st["C"],st["O"],st["ER"],st["Fh"],SIGMA)
        book.reserve(len(S),40,1,"analytic complete occupancy/property derivative")
        out=solver.evaluate(a,delta,need_jac=True)
        book.reserve(256,8,1,"independent full-current state for implicit tangent")
        fac=lu_factor(np.eye(256)-delta*a[:,None]*st["C"],check_finite=False)
        h=lu_solve(fac,delta*a[:,None]*st["ER"],check_finite=False)
        e=st["ER"]+st["C"]@h
        book.reserve(256,8,0,"independent support-birth tangent")
        dh=lu_solve(fac,delta*da[:,None]*e,check_finite=False)
        birth=(st["O"]@dh).T.reshape(-1)
        book.reserve(256,8,0,"independent holomorphic property tangent")
        dhd=lu_solve(fac,a[:,None]*e,check_finite=False)
        prop=(st["O"]@dhd).T.reshape(-1)
        step=1e-4;anew=a+step*da
        book.reserve(len(S)+1,8,1,"actual positive occupancy birth")
        new=solver.evaluate(anew,delta)
        actual=((new["F"]-out["F"])/step).T.reshape(-1)
        result=dict(status="COMPLETE",shape=identity,cell=j,a_before=float(a[j]),a_after=float(anew[j]),
            support_before=len(S),support_after=len(S)+1,
            exact_support_rule="a>0, no truncation",delta=[delta.real,delta.imag],step=step,
            inactive_column_norm=float(np.linalg.norm(out["J_occ"][:,j])),
            independent_implicit_birth_relative_error=float(np.linalg.norm(out["J_occ"][:,j]-birth)/np.linalg.norm(birth)),
            independent_implicit_property_relative_error=float(np.linalg.norm(out["J_delta"]-prop)/np.linalg.norm(prop)),
            positive_birth_linearization_relative_remainder=float(np.linalg.norm(actual-birth)/np.linalg.norm(birth)),
            zeroed_inactive_derivative_would_have_relative_error=1.,
            interpretation="one real support birth and independent full-state implicit derivatives; finite-step remainder is not an exact equality test",GT_access=False)
        write(HERE/"SUPPORT_BIRTH_CHECK.json",result);book.event("SUPPORT_BIRTH_RESULT",**result);book.close("COMPLETE")
        print(json.dumps(result,indent=2))
    except BaseException as e:
        book.close("FAILED",repr(e));raise


if __name__=="__main__":
    with threadpool_limits(limits=1):run()
