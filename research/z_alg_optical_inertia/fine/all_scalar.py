"""Two complementary witnesses excluding every global scalar multiplier on the fine grid.

Only auxiliary currents are restricted to block constants to find witnesses. A
counterexample in that subspace also belongs to the full 9216-cell current space.
"""
import json, sys, time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.linalg import eigh, qr
from threadpoolctl import threadpool_limits
from flint import arb,acb,acb_mat,ctx
from fourier_bound import STUDY,PUBLIC,DATA,h
sys.path.insert(0,str(STUDY/'core'))
from certificate import ball,star,frobenius_squared,gershgorin,safe_float
from single_source import fit_proposals

def lo(x):return safe_float(x,upper=False)
def hi(x):return safe_float(x,upper=True)

def compress(spec, S, B):
    # U repeats each coarse current on 3x3 subcells with value 1/3.
    p=len(spec);m=32;axes=np.arange(1-m,m); klen=len(axes)
    fine=[acb.dft([acb(float(z.real),float(z.imag)) for z in row],inverse=True) for row in spec]
    tab=[acb.dft([fine[y][x] for y in range(p)],inverse=True) for x in range(p)]
    g=[]
    for dy in axes:
        row=[]
        for dx in axes:
            row.append(sum((tab[(3*dx+j)%p][(3*dy+i)%p]*((3-abs(j))*(3-abs(i)))
                       for i in range(-2,3) for j in range(-2,3)),acb(0))/9)
        g.append(row)
    sy,sx=np.indices((m,m)).reshape(2,-1)
    index=(sy[:,None]-sy[None,:]+m-1)*klen+(sx[:,None]-sx[None,:]+m-1)
    ft=[x for z in g for x in z]
    G=acb_mat(m*m,m*m,[ft[i] for i in index.ravel()])
    gn=np.array([complex(x) for x in ft])[index]
    def sum_blocks(X):
        r=X.reshape(m,3,m,3,-1).transpose(0,2,4,1,3).reshape(m*m,X.shape[1],9)
        values=[sum((acb(float(z.real),float(z.imag)) for z in nine),acb(0))/3 for nine in r.reshape(-1,9)]
        return acb_mat(r.shape[0],r.shape[1],values),np.array([complex(z) for z in values]).reshape(r.shape[:2])
    S0,s0=sum_blocks(S.T); B0,b0=sum_blocks(B)
    return G,S0.transpose(),B0,gn,s0.T,b0

def prove():
    start=time.perf_counter();ctx.prec,ctx.threads=160,1
    with np.load(PUBLIC/'fine_geometry.npz') as f:spec=f['spectrum'];S=f['PW360_10x20__S'];B=f['PW360_10x20__B']
    G,S,B,g,s,b=compress(spec,S,B); ss=star(S)
    u=.75;t=10000.
    Q=qr(s.conj().T,mode='full')[0][:,len(s):]
    C=np.eye(len(g))-u*h(g)+t*h(g/1j)
    _,v=eigh(h(Q.conj().T@C@Q),subset_by_index=(0,0))
    x=ball(Q@v);v=x-ss*(S*ss).solve(S*x)
    vg=(star(v)*(G*v))[0,0]
    dark=hi(frobenius_squared(v)-arb(u)*vg.real+arb(t)*vg.imag)
    rlo=lo(vg.imag)
    direct=json.loads((DATA/'fine_long_ratio.json').read_text())
    if direct['cap']!=u or direct['ratio']!=t:raise ValueError('Dark proof has different parameters')
    full_dark=direct['verified_form_upper']
    full_rad=direct['radiation_lower']
    report=dict(native_dark_proof='fine_long_ratio.json',native_dark_at_ratio_upper=full_dark,
                native_dark_radiation_lower=full_rad,geometry='96² fine PW360 10x20',cap=u,ratio=t,dark_at_ratio_upper=dark,
                exact_hidden_radiation_lower=rlo,cases={})
    meta=json.loads((PUBLIC/'manifest.json').read_text())['observations']
    for obj in ('01','05'):
        case='optfresh_'+obj+'__pw360_10x20__n02'
        radius=next(z['radius_noise_fro'] for z in meta if z['id']==case)
        with np.load(PUBLIC/(case+'.npz')) as f:y=f['Y']
        p=SimpleNamespace(S=s,B=b,G=g,Ga=h(g/1j),Y=y,n=len(g),s=b.shape[1],r=s.shape[0])
        best,_,records,_=fit_proposals(p,radius,u,t)
        row,X=best
        xx=ball(X);jE=star(xx)*(B+G*xx)
        qi=(star(jE)-jE)/acb(0,2);qr0=star(xx)*xx-arb(u)*(jE+star(jE))/2;D=arb(t)*qi-qr0
        tests={}
        for key,A in [('T',qi),('tT_minus_Qr',D)]:
            V=ball(eigh(h(np.array([[complex(z) for z in rr] for rr in A.tolist()])))[1])
            tests[key]=lo(gershgorin(star(V)*A*V,upper=False))
            tests[key+'_basis']=lo(gershgorin(star(V)*V,upper=False))
        error=hi(frobenius_squared(S*xx-ball(y)).sqrt())
        report['cases'][case]=dict(data_fit_upper=error,radius=radius,strict_signs=tests,proposal=row,
                all_mu_and_all_real_lambda_fail_for_single_source=bool(full_dark<0 and full_rad>0 and min(tests.values())>0 and error<radius))
        print(case,dark,rlo,error,radius,tests,flush=True)
    report['elapsed_seconds']=time.perf_counter()-start
    (DATA/'fine_all_scalar.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':
    with threadpool_limits(1):prove()
