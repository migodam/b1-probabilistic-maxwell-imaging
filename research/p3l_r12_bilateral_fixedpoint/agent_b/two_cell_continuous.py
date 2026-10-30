"""Actual three-real-DOF / one-complex-port two-cell continuous control.

This is a separate low-dimensional passive-material contract, NOT original P1.
The explicit fixed-response curve and uniform elementary bounds establish a
local finite-dimensional family. Paid node solves verify its implementation.
"""
from pathlib import Path
import json
import time
import numpy as np
from scipy.linalg import lu_factor, lu_solve, cho_factor, cho_solve
from two_cell_controls import Book, module, SIGMA


def main():
    out = Path(__file__).resolve().parent/'two_cell_continuous'
    out.mkdir(exist_ok=False)
    book = Book()
    model = module.Model(book=book, centers=np.array([[-.12,-.02],[.15,.03]]),
        sources=np.array([[-1.,-.6]]), receivers=np.array([[.9,-.8]]), radius=.04)
    fixed_loss1 = .08
    def physical(u, jac=False, label='CURVE'):
        m = np.array([u[0]+1j*fixed_loss1, u[1]+1j*u[2]])
        with book.op('COMPLEX_TWO_CELL_FORWARD', rhs=1, LU=1, label=label):
            A = np.eye(2)-m[:,None]*model.G
            fac = lu_factor(A, check_finite=False)
            current = lu_solve(fac, m[:,None]*model.E0, check_finite=False)
            F = (model.S@current).item()
        q = dict(u=np.asarray(u).copy(), m=m, F=F, A=A)
        if jac:
            with book.op('COMPLEX_TWO_CELL_ADJOINT', rhs=1, adjoint_RHS=1, label=label):
                transfer = lu_solve(fac, model.S.T, trans=1, check_finite=False).T
                field = model.E0+model.G@current
                Jm = (transfer*field.T).ravel()
                Dm = np.array([[1.,0.,0.],[0.,1.,1j]])
                Jc = Jm@Dm
                q['J'] = np.vstack((Jc.real,Jc.imag))*np.sqrt(2)/SIGMA
        return q
    u0 = np.array([.65,1.1,.2])
    base = physical(u0, True, 'CURVE_BASE')
    target = base['F']
    G, e, s = model.G, model.E0.ravel(), model.S.ravel()
    a, b = s[0]*e[0], s[1]*e[1]
    c = -s[0]*G[1,1]*e[0]+s[0]*G[0,1]*e[1]+s[1]*G[1,0]*e[0]-s[1]*G[0,0]*e[1]
    detG = G[0,0]*G[1,1]-G[0,1]*G[1,0]
    m10,m20 = base['m']
    n0 = target-(target*G[0,0]+a)*m10
    n1 = -(target*G[0,0]+a)
    d1 = c-target*detG
    d0 = target*G[1,1]+b+d1*m10
    derivative_numerator = n1-m20*d1
    sigma_A = np.linalg.svd(base['A'], compute_uv=False)[-1]
    norm_G = np.linalg.norm(G,2)
    # Analytic a-priori interval rule, determined before any nonbase curve solve.
    delta_bound_target = min(.25*min(m20.real,m20.imag), .25*sigma_A/norm_G)
    tau = min(.05, .8*delta_bound_target,
              .8*delta_bound_target*abs(d0)/(abs(derivative_numerator)+delta_bound_target*abs(d1)))
    denominator_lower = abs(d0)-tau*abs(d1)
    m2_change_bound = tau*abs(derivative_numerator)/denominator_lower
    resolvent_neumann_bound = max(tau,m2_change_bound)*norm_G/sigma_A
    real_lower, loss_lower = m20.real-m2_change_bound, m20.imag-m2_change_bound
    def exact_material(t):
        m2 = (n0+n1*t)/(d0+d1*t)
        return np.array([m10+t,m2])
    frozen = dict(contract='NOT_P1_COMPLEX_TWO_CELL_3REAL_1COMPLEXPORT',
        freeze_before_nonbase_nodes=time.time(), free_interval=[-tau,tau],
        fixed_loss1=fixed_loss1, target=[target.real,target.imag],
        proof='explicit Mobius curve, denominator triangle bound, passive-material bound, Neumann resolvent bound')
    (out/'CONTROL_CONTRACT.json').write_text(json.dumps(frozen,indent=2)+'\n')
    nodes=[]
    for t in np.linspace(-tau,tau,9):
        material=exact_material(t)
        u=np.array([material[0].real,material[1].real,material[1].imag])
        q=physical(u,False,'EXPLICIT_CURVE_NODE')
        residual=abs(q['F']-target)*np.sqrt(2)/SIGMA
        nodes.append(dict(t=float(t),u=u.tolist(),material_real=material.real.tolist(),
                          material_imag=material.imag.tolist(),whitened_residual=float(residual)))
    corrected=[]
    tangent = derivative_numerator/d0
    for t in (-tau,tau):
        u=u0+np.array([t,tangent.real*t,tangent.imag*t])
        trace=[]
        for iteration in range(6):
            q=physical(u,True,'ORDINARY_TRANSVERSE_NEWTON')
            rr=np.array([(q['F']-target).real,(q['F']-target).imag])*np.sqrt(2)/SIGMA
            trace.append(dict(iteration=iteration,u=u.tolist(),residual=float(np.linalg.norm(rr))))
            if np.linalg.norm(rr)<1e-10:break
            B=q['J'][:,1:]
            dz=-cho_solve(cho_factor(B.T@B+1e-14*np.eye(2),lower=True),B.T@rr)
            u[1:]+=dz
        corrected.append(dict(t=float(t),free_coordinate_defect=float(abs(u[0]-u0[0]-t)),
            final_residual=trace[-1]['residual'],trace=trace))
    results=dict(status='PASS',actual_discrete_Maxwell=True,
        contract=frozen,rank_at_base=int(np.linalg.matrix_rank(base['J'])),
        transverse_singular_values=np.linalg.svd(base['J'][:,1:],compute_uv=False).tolist(),
        uniform_bounds=dict(denominator_modulus_lower=float(denominator_lower),
            material2_change_modulus_upper=float(m2_change_bound),
            material1_real_lower=float(m10.real-tau),material1_fixed_loss=fixed_loss1,
            material2_real_lower=float(real_lower),material2_loss_lower=float(loss_lower),
            resolvent_neumann_upper=float(resolvent_neumann_bound),
            base_resolvent_sigma_min=float(sigma_A),
            scope='analytic inequalities using stored floating coefficients; not outward-rounded interval arithmetic'),
        node_max_whitened_residual=max(q['whitened_residual'] for q in nodes),
        nodes=nodes,transverse_corrected=corrected,
        cost=dict(RHS=book.rhs,LU=book.lu,tickets=book.tickets),
        theorem_scope='finite two-cell coefficient model; one-real-dimensional passive fixed-response curve on declared interval; NOT P1, NOT Full66, NOT continuous PDE, NOT novelty claim')
    if denominator_lower<=0 or loss_lower<=0 or real_lower<=0 or resolvent_neumann_bound>=1 or results['rank_at_base']!=2 or results['node_max_whitened_residual']>1e-8:
        results['status']='CONTROL_CONDITION_NOT_MET'
    np.savez_compressed(out/'ARRAYS.npz',G=G,E=model.E0,S=model.S,
        base_material=base['m'],base_J=base['J'],target=np.array(target),
        mobius_coefficients=np.array([n0,n1,d0,d1]),tau=np.array(tau))
    (out/'RESULTS.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:results[k] for k in ('status','rank_at_base','uniform_bounds','node_max_whitened_residual','cost')}))


if __name__=='__main__':main()
