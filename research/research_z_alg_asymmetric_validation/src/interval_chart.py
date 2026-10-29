"""Validated local graph for the exact stored binary64 Green-coefficient model.

Hankel coefficient approximation and continuous-PDE discretization are separate.
No claim that this covers the full legal fibre or its conditional probability.
"""
import json,math
import numpy as np
from mpmath import iv
from accounting import P,Book,write
from physics import model_for
from charts import partition
iv.dps=40

def ci(z):return iv.mpc(iv.mpf(float(np.real(z))),iv.mpf(float(np.imag(z))))
def mat(A):return iv.matrix([[ci(z) for z in row] for row in np.asarray(A)])
def upper(x):return math.nextafter(float(x.b),math.inf)
def lower(x):return math.nextafter(float(x.a),-math.inf)
def absup(z):return upper(abs(z))
def absmat(A):return np.asarray([[absup(A[i,j]) for j in range(A.cols)] for i in range(A.rows)])
def rm(A):return iv.matrix([[iv.mpf(float(v)) for v in row] for row in np.asarray(A)])
def isum(values):
    v=iv.mpf(0)
    for a in values:v+=a
    return v

def enclosure(G,S,E0,m,width,R):
    n=len(m);ms=[iv.mpf(float(m[i]))+iv.mpf([-float(width[i]),float(width[i])]) for i in range(n)]
    gg=mat(G);ss=mat(S);ee=mat(E0);rr=mat(R)
    A=iv.eye(n)
    for i in range(n):
        for j in range(n):A[i,j]-=ms[i]*gg[i,j]
    residual=iv.eye(n)-rr*A;M=absmat(residual)
    d=max(upper(isum(iv.mpf(float(x)) for x in row)) for row in M)
    if d>=.8:return None
    Miv=rm(M);RA=absmat(rr);power=Miv*rm(RA);error=power.copy()
    for _ in range(7):power=Miv*power;error+=power
    rnorm=max(upper(isum(iv.mpf(float(x)) for x in row)) for row in RA)
    tail=upper(iv.mpf(rnorm)*iv.mpf(d)**9/(1-iv.mpf(d)))
    inv=iv.matrix(n,n)
    for i in range(n):
        for j in range(n):
            e=upper(error[i,j]+iv.mpf(tail));z=R[i,j]
            inv[i,j]=ci(z)+iv.mpc(iv.mpf([-e,e]),iv.mpf([-e,e]))
    rhs=iv.matrix(n,E0.shape[1])
    for i in range(n):
        for s in range(E0.shape[1]):rhs[i,s]=ms[i]*ee[i,s]
    j=inv*rhs;E=ee+gg*j;y=ss*j;transfer=ss*inv
    J=iv.matrix(2*E0.shape[1],n);f=iv.matrix(2*E0.shape[1],1)
    for s in range(E0.shape[1]):
        f[2*s,0]=y[0,s].real;f[2*s+1,0]=y[0,s].imag
        for i in range(n):
            z=transfer[0,i]*E[i,s];J[2*s,i]=z.real;J[2*s+1,i]=z.imag
    return dict(J=J,f=f,Neumann_norm=d)

def certificate(model,m,J,target,dep,free):
    n=len(m)
    with model.book.op('INTERVAL_INVERSE_PRECONDITIONER',rhs=n,LU=1,validation_inverse_columns=n):
        R=np.linalg.inv(np.eye(n)-m[:,None]*model.G)
    with model.book.op('INTERVAL_REAL_J_PRECONDITIONER',algebraic_real_RHS=4,algebraic_real_factorizations=1):
        Rj=np.linalg.inv(J[:,dep])
    riv=rm(Rj);target=np.asarray(target);point=enclosure(model.G,model.S,model.E0,m,np.zeros(n),R)
    g=point['f']-rm(target[:,None]);center_res=absmat(riv*g)[:,0]
    attempts=[]
    for level in range(13):
        hf=.05/(2**level)
        sensitivity=abs(Rj@J[:,free])@np.full(len(free),hf)
        rd=2*sensitivity+1e-9
        width=np.zeros(n);width[dep]=rd;width[free]=hf
        if np.any(m-width<=.5) or np.any(m+width>=4):attempts.append(dict(level=level,halfwidth=hf,status='LEGAL_BOX_TOO_WIDE'));continue
        with model.book.op('INTERVAL_GRAPH_BOUND',interval_boxes=2,algebraic_certificate=1):
            whole=enclosure(model.G,model.S,model.E0,m,width,R)
            wf=np.zeros(n);wf[free]=hf
            base=enclosure(model.G,model.S,model.E0,m,wf,R)
            if whole is None or base is None:attempts.append(dict(level=level,halfwidth=hf,status='NEUMANN_BOUND_TOO_WIDE'));continue
            D=iv.matrix(4,4);F=iv.matrix(4,len(free))
            for i in range(4):
                for j,k in enumerate(dep):D[i,j]=whole['J'][i,int(k)]
                for j,k in enumerate(free):F[i,j]=base['J'][i,int(k)]
            N=absmat(iv.eye(4)-riv*D);RF=absmat(riv*F)
            ratios=[];qs=[]
            for i in range(4):
                movement=iv.mpf(float(center_res[i]))+isum(iv.mpf(float(RF[i,j]))*iv.mpf(hf) for j in range(len(free)))
                nonlin=isum(iv.mpf(float(N[i,j]))*iv.mpf(float(rd[j])) for j in range(4))
                ratios.append(upper((movement+nonlin)/iv.mpf(float(rd[i]))));qs.append(upper(nonlin/iv.mpf(float(rd[i]))))
            success=max(ratios)<1 and max(qs)<1
            attempts.append(dict(level=level,halfwidth=hf,dependent_radius=rd.tolist(),mapping_inclusion_ratios=ratios,contraction_ratios=qs,Neumann_norm=whole['Neumann_norm'],status='CERTIFIED_FROZEN_COEFFICIENT_LOCAL_GRAPH' if success else 'BOUND_TOO_WIDE'))
            if success:return dict(status=attempts[-1]['status'],free=free.tolist(),dependent=dep.tolist(),material_center=m.tolist(),attempts=attempts,
                coefficient_scope='Exact binary64 constants enclosed with mpmath interval outward rounding; no rigorous Hankel-transcendental coefficient enclosure or continuous PDE certificate',
                implication='For every free-coordinate point in this box, a unique dependent-coordinate root exists in the declared dependent box; no other global branch excluded')
    return dict(status='NOT_CERTIFIED_WITHIN_PREDECLARED_SHRINK_LEVELS',attempts=attempts)

def main():
    book=Book('E1_INTERVAL_GRAPH_CERTIFICATES');model=model_for(book);c=json.loads((P/'EXPERIMENT_CONTRACT.json').read_text())
    rows=[]
    for i,m in enumerate(c['E1']['public_reference_materials']):
        state=model.forward(np.asarray(m),jac=True,label='certificate_center_'+str(i));dep,free=partition(state['J'])
        result=certificate(model,np.asarray(m),state['J'],model.pack(state['y']),dep,free);result['parent']=i;rows.append(result)
    write(P/'audit/INTERVAL_CHART_CERTIFICATES.json',rows);book.end();print(json.dumps(rows,indent=2))
if __name__=='__main__':main()
