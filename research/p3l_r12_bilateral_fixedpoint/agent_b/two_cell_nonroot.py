"""Additional asymmetric two-cell physical residual-curvature control."""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import null_space
from two_cell_controls import Book, module, inverse_spd, SIGMA


def main():
    out = Path(__file__).resolve().parent/'two_cell_nonroot'
    out.mkdir(exist_ok=False)
    book = Book()
    model = module.Model(book=book,
        centers=np.array([[-.12, -.02], [.15, .03]]),
        sources=np.array([[-1., -.6], [1., .7]]),
        receivers=np.array([[-.8, .9], [.9, -.8], [.3, 1.2]]), radius=.04)
    x = np.array([.65, 1.1])
    base = model.forward(x, label='ASYMMETRIC_BASE', jac=True)
    y, J = base['y'].copy(), base['J'].copy()*np.sqrt(2)/SIGMA
    h = 1e-5
    nodes = []
    DJ = []
    for i in range(2):
        ei = np.eye(2)[i]
        qp = model.forward(x+h*ei, label=f'ASYMMETRIC_PLUS_{i}', jac=True)
        yp, jp = qp['y'].copy(), qp['J'].copy()*np.sqrt(2)/SIGMA
        qm = model.forward(x-h*ei, label=f'ASYMMETRIC_MINUS_{i}', jac=True)
        ym, jm = qm['y'].copy(), qm['J'].copy()*np.sqrt(2)/SIGMA
        nodes.append((yp,jp,ym,jm))
        DJ.append((jp-jm)/(2*h))
    N = null_space(J.T)
    curvature = np.column_stack([d[:,j] for d in DJ for j in range(2)])
    u, strength, _ = np.linalg.svd(N.T@curvature, full_matrices=False)
    r = N@u[:,0]
    observation = y-model.unpack(r*SIGMA/np.sqrt(2))
    A = J.T@J+np.eye(2)
    g = J.T@r
    gnstep = inverse_spd(A, g)
    exact, omitted, fd = [], [], []
    for i in range(2):
        ei = np.eye(2)[i]
        DA = DJ[i].T@J+J.T@DJ[i]
        exact.append(ei-inverse_spd(A,J.T@J@ei+DJ[i].T@r)+inverse_spd(A,DA@gnstep))
        omitted.append(ei-inverse_spd(A,J.T@J@ei)+inverse_spd(A,DA@gnstep))
        yp,jp,ym,jm = nodes[i]
        rp = model.pack(yp-observation)*np.sqrt(2)/SIGMA
        rm = model.pack(ym-observation)*np.sqrt(2)/SIGMA
        tp = x+h*ei-inverse_spd(jp.T@jp+np.eye(2),jp.T@rp)
        tm = x-h*ei-inverse_spd(jm.T@jm+np.eye(2),jm.T@rm)
        fd.append((tp-tm)/(2*h))
    exact, omitted, fd = map(np.column_stack,(exact,omitted,fd))
    result = dict(status='PASS', actual_discrete_Maxwell=True, mu=1.,
        residual_norm=float(np.linalg.norm(r)), gradient_norm=float(np.linalg.norm(g)),
        residual_curvature_singular_values=strength.tolist(),
        exact_DT_FD_error=float(np.linalg.norm(fd-exact)),
        omitted_residual_term_FD_error=float(np.linalg.norm(fd-omitted)),
        residual_term_norm=float(np.linalg.norm(exact-omitted)),
        DT_eigenvalues=np.linalg.eigvals(exact).real.tolist(),
        cost=dict(RHS=book.rhs,LU=book.lu,tickets=book.tickets),
        interpretation='stationary positive-residual point; DJ transpose r term is measurable; no exact-data root or fiber claim')
    if result['exact_DT_FD_error'] > 1e-6 or result['residual_term_norm'] < 1e-5:
        result['status']='CONTROL_CONDITION_NOT_MET'
    np.savez_compressed(out/'ARRAYS.npz',J=J,DJ=np.array(DJ),residual=r,
                        exact_DT=exact,omitted_DT=omitted,FD_DT=fd,observation=observation)
    (out/'RESULTS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
