"""Material-branch figure, coordinates test, and the data-derived minimax floor."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from flint import arb,acb,ctx
from host_model import reduced_forward
from certify_fold import pulse_operators, reduce_host, positive_ldl_pivots

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'
FIG=ROOT/'figures'


def coordinates_check(H,B,R,z,I):
    # Both t and c are grouped as (all real, all imaginary) for the determinant.
    rest=[i for i in range(len(z)) if i not in I]
    N=np.zeros((3,6),complex)
    N[:,rest]=np.eye(3)
    N[:,I]=-B[rest]@np.linalg.inv(B[I])
    A=np.eye(6)-z[:,None]*H
    V=(R@np.linalg.inv(A)).ravel()
    W=V*z
    cJ=[]
    for q in np.vstack([N,1j*N]):
        cJ.append((q*V.conj()+W*(q@H).conj()).imag)
    DJ=np.array(cJ).T
    E=np.linalg.solve(np.eye(6)-H*z[None,:],B)
    J=(V[:,None]*E).T
    realJ=np.vstack((J.real,J.imag))
    pred=(-1)**3*(abs(np.linalg.det(A))/abs(np.linalg.det(B[I])))**2*np.linalg.det(realJ)
    found=np.linalg.det(DJ)
    h=np.array([.7+1j*.1,-.1+.4j,.3-1j*.2])@N
    qH=h@H
    step=.073
    def C(delta):
        return np.imag((W+delta*h)*(V+delta*qH).conj())
    numeric=(C(step)-2*C(0)+C(-step))/step**2
    formula=2*(h*qH.conj()).imag
    return dict(z=z.tolist(),determinant_measured=found,determinant_formula=pred,
                determinant_relative_error=abs(found-pred)/abs(pred),
                exact_quadratic_central_difference_max_error=max(abs(numeric-formula)))


def main():
    data=np.load(OUT/'host_operators.npz')
    proof=json.loads((OUT/'normal_section.json').read_text())
    count=json.loads((OUT/'point_fibers.json').read_text())
    cert=json.loads((OUT/'fold_certificate.json').read_text())
    mat=[a['material_midpoint'] for a in count['fibers']['minus']['physical_roots']]
    z=np.array(mat)
    z=z[np.argsort(z[:,-1])]
    zc,v=data['zc'],data['null']
    cutoff=np.mean(list(cert['null_hessian_pairing'].values()))
    ctx.prec=128
    separated=[]
    a,b=[r['material_intervals'] for r in count['fibers']['minus']['physical_roots']]
    for x,y in zip(a,b):
        separated.append((arb(x[0]).union(arb(x[1]))-arb(y[0]).union(arb(y[1]))).abs_lower())
    dist=sum(x*x for x in separated).sqrt()
    checks=[coordinates_check(data['H'],data['B'],data['R'],x,count['fibers']['minus']['pivot'])
            for x in (z[0],z[1],np.full(6,2.))]
    ctx.prec=192
    G,E0,S=pulse_operators()
    H,_,_,_,_,_,U=reduce_host(G,E0,S)
    imag=lambda M:(M-M.conjugate().transpose())/acb(0,2)
    diff=imag(H)-U.conjugate().transpose()*imag(G)*U
    assert all(q.contains(0) for q in diff.entries())
    smallest=min(x.lower() for x in positive_ldl_pivots(imag(H)))
    gram_check={'entrywise_congruence_difference_abs_upper':float(np.nextafter(
                    float(max(q.abs_upper() for q in diff.entries())),np.inf)),
                'reduced_imaginary_green_min_LDL_pivot_lower':float(np.nextafter(float(smallest),-np.inf))}
    result=dict(pair_euclidean_lower=float(np.nextafter(float(dist.lower()),-np.inf)),
                pair_rms_lower=float(np.nextafter(float((dist/arb(6).sqrt()).lower()),-np.inf)),
                minimax_rms_lower=float(np.nextafter(float((dist/(2*arb(6).sqrt())).lower()),-np.inf)),
                square_root_prediction=float(np.sqrt(2*data['tau']/abs(cutoff))),
                observed_weak_axis_displacements=((z-zc)@v).tolist(),coordinate_tests=checks,
                host_imaginary_green_congruence=gram_check)
    (OUT/'branch_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    fig,(ax,bx)=plt.subplots(1,2,figsize=(10,4.1),layout='constrained',gridspec_kw={'width_ratios':[1.27,1]})
    t=np.array([q['t'] for q in proof['curve']])
    p=np.array([np.mean(q['p']) for q in proof['curve']])*1e6
    keep=abs(t)<.32
    # Line segments interpolate the verified nodes for display.
    ax.axvspan(-1.,1.,color='#e7edf2')
    ax.plot(p[keep],t[keep],'-',color='#16546d',lw=2,label='Validated section nodes')
    ax.scatter(p[keep],t[keep],s=8,color='#16546d')
    t0=np.linspace(-.32,.32,240)
    ax.plot(.5*cutoff*(t0/v[-1])**2*1e6,t0,':',lw=1.8,color='#71867d',label='Local quadratic prediction')
    for s,color in [(-1.,'#ab4b2c'),(1.,'#787c89')]:
        ax.axvline(s,ls='--',lw=.9,color=color)
    ax.scatter([-1.,-1.],z[:,-1]-zc[-1],s=45,color='#ab4b2c',zorder=5)
    ax.scatter([0.],[0.],color='#122d36',marker='D',s=30,zorder=5)
    ax.text(.54,.16,'0 roots',color='#5b6470',fontsize=10)
    ax.text(-1.06,.263,'2 roots',color='#8a3a20',fontsize=10,ha='center')
    ax.set(xlim=(-1.95,1.25),ylim=(-.32,.32),xlabel='Change in measured data  (×10⁻⁶)',
           ylabel='Slot 6 contrast relative to the fold',title='One turning point on the measured line')
    ax.legend(loc='lower right',frameon=False,fontsize=8)
    x=np.arange(1,7)
    bx.axhline(0,color='0.65',lw=.8)
    for i,color,label in ((0,'#ad5030','Lower slot 6'),(1,'#25677b','Higher slot 6')):
        bx.bar(x+(i-.5)*.35,z[i]-zc,width=.33,color=color,label=label)
    bx.set(xlabel='Unknown material slot',ylabel='Contrast relative to the fold',
           title='Two materials for the same three readings',xticks=x)
    bx.legend(frameon=False,fontsize=9)
    for axis in (ax,bx):
        axis.spines[['top','right']].set_visible(False)
    for ext in ('pdf','png'):
        fig.savefig(FIG/('material_branches.'+ext),dpi=210)
    print(result['pair_rms_lower'],result['minimax_rms_lower'],checks[0]['determinant_relative_error'])

if __name__=='__main__':
    main()
