import json
import numpy as np
import mpmath as mp
from accounting import P,Book,write

def main():
    c=json.loads((P/'EXPERIMENT_CONTRACT.json').read_text());b=Book('INDEPENDENT_70DPS_FULL_MAXWELL');mp.mp.dps=70
    centers=c['E1']['centers'];sources=c['E1']['sources'];receivers=c['E1']['receivers'];n=6
    v=lambda x:mp.mpf(str(x));k=2*mp.pi;a=v(c['disk_radius'])
    def distance(x,y):return mp.sqrt(sum((v(i)-v(j))**2 for i,j in zip(x,y)))
    def external(d):return mp.j*mp.pi*k*a/2*mp.besselj(1,k*a)*mp.hankel1(0,k*d)
    with b.op('HIGH_PRECISION_HANKEL_SETUP',Green_kernel_evaluations=54,decimal_digits=70):
        G=mp.matrix(n,n);S=mp.matrix(1,n);E0=mp.matrix(n,2)
        for i in range(n):
            for j in range(n):G[i,j]=mp.j*mp.pi*k*a/2*mp.hankel1(1,k*a)-1 if i==j else external(distance(centers[i],centers[j]))
            S[0,i]=external(distance(centers[i],receivers[0]))
            for s in range(2):E0[i,s]=mp.j/4*mp.hankel1(0,k*distance(centers[i],sources[s]))
    rows=[]
    for parent in range(3):
        with np.load(P/f'results/E1/P{parent}/A/REFERENCE.npz') as z:originalG=z['G'].copy();reference=z['y'].copy()
        for role,path in [('reference',P/f'results/E1/P{parent}/A/REFERENCE.npz'),('preselected_radial_0_01',P/f'results/E1/P{parent}/A/RADIAL_0_01/ENDPOINT.npz')]:
            with np.load(path) as z:m=z['m'].copy();oldy=z['y'] if 'y' in z else z['y_full']
            L=mp.eye(n);rhs=mp.matrix(n,2)
            for i in range(n):
                for j in range(n):L[i,j]-=v(m[i])*G[i,j]
                for s in range(2):rhs[i,s]=v(m[i])*E0[i,s]
            currents=mp.matrix(n,2)
            for s in range(2):
                with b.op('HIGH_PRECISION_FULL_SOURCE_SOLVE',rhs=1,LU=1,new_states=1 if s==0 else 0,source=s,decimal_digits=70):
                    solved=mp.lu_solve(L,rhs[:,s])
                    for i in range(n):currents[i,s]=solved[i]
            y=S*currents;Ey=E0+G*currents
            asarray=lambda M:np.array([[complex(M[i,j]) for j in range(M.cols)] for i in range(M.rows)],dtype=np.complex128)
            yd=asarray(y);gd=asarray(G)
            rel=np.linalg.norm(yd-oldy)/np.linalg.norm(oldy);pair=np.linalg.norm(yd-reference)/np.linalg.norm(reference)
            residual=mp.norm(L*currents-rhs)/mp.norm(rhs)
            row=dict(parent=parent,role=role,relative_to_binary64=float(rel),pair_relative_to_reference=float(pair),state_relative_70dps=mp.nstr(residual,60),
                Green_relative=float(np.linalg.norm(gd-originalG)/np.linalg.norm(originalG)),independent_exact_decimal_geometry=True,
                qualification='Numerical 70-digit recomputation, not rigorous transcendental interval or continuum discretization certificate')
            rows.append(row)
            b.save(b.events,m=m,j=asarray(currents),E=asarray(Ey),y=yd,G=gd,S=asarray(S),E0=asarray(E0))
    out=dict(status='PASS' if max(r['relative_to_binary64'] for r in rows)<1e-10 else 'FAIL',rows=rows,new_complex_RHS=12,LU_factorizations=12)
    write(P/'audit/HIGH_PRECISION_FULLWAVE.json',out);b.end(out['status']);print(json.dumps(out,indent=2))
if __name__=='__main__':main()
