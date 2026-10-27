"""One explicitly paid 8+32 RHS replay to retain complete coefficient arrays."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from shared.physics import OUT,Book,Target,write,save

def main():
    b=Book('A','compensation_saved_vectors',role='PHYSICS_MECHANISM',max_wall=120)
    status='COMPLETE';error=None
    try:
        t=Target(b);q=np.r_[np.zeros(64),.4,.45];q[0]=1.
        p=t.evaluate(q,'paid_compensation_vector_replay');J=t.jacobian(p)
        t.persist_cache(p,OUT/'track_a/raw/COMPENSATION_BASE_CACHE.npz')
        c=-2.74625*(1.2*J[:,64]+.0975/.15*J[:,65]);JV=np.column_stack((J[:,64],J[:,65]/.15))
        rows=[]
        # Every selected datum retains its real and imaginary component. This
        # confirms the stated complex-Hilbert/conformal qualification.
        for nsource in [1,2,4,8]:
            keep=np.arange(nsource*64);A=JV[keep];coef=c[keep]
            gram=A.T@A;scale=np.trace(gram)/2
            q2,_=np.linalg.qr(A);q1=A[:,:1]/np.linalg.norm(A[:,:1])
            rows.append(dict(sources=nsource,complex_receivers_each=32,
                conformal_Gram_relative_error=float(np.linalg.norm(gram-scale*np.eye(2))/scale),
                two_property_normal_ratio=float(np.linalg.norm(coef-q2@(q2.T@coef))/np.linalg.norm(coef)),
                one_property_normal_ratio=float(np.linalg.norm(coef-q1@(q1.T@coef))/np.linalg.norm(coef))))
        save(OUT/'track_a/raw/COMPENSATION_VECTORS.npz',c_white=c,J_material_white=JV,Dg=t.prior.derivative(q),q=q)
        write(OUT/'track_a/COMPENSATION_VECTOR_AUDIT.json',dict(rows=rows,new_RHS=40,
            theoretical_one_property_ratio=.0975/abs(complex(1.2,.0975)),
            qualification='discrete float64, original P1 and complete complex-channel subsets; standard holomorphic response consequence, not global field conformality'))
    except BaseException as exc:status='FAILED';error=repr(exc);raise
    finally:b.close(status,error)

if __name__=='__main__':main()
