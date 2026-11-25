"""Independent Arb recomputation of selected fixed-trace exclusions."""
import argparse
import json
import os

for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[variable] = "1"

from flint import acb, acb_mat, arb, ctx
import numpy as np
from threadpoolctl import threadpool_limits

from common import FREEZE, HERE, Geometry, clock, elapsed, geom_path, pattern, place, save_json
from operator_bounds import ball, inverse_dft_table
from sector_map import lag_tables


def as_ball(a):
    return acb_mat(*a.shape,[ball(z) for z in a.flat])


def adjoint(a):
    return a.conjugate().transpose()


def identity(n):
    return acb_mat(n,n,[int(i==j) for i,j in np.ndindex(n,n)])


def diag(v):
    return acb_mat(len(v),len(v),[float(v[i]) if i==j else 0
                                 for i,j in np.ndindex(len(v),len(v))])


def norm(a):
    return sum((abs(v).upper()**2 for v in a.entries()),arb(0)).sqrt().upper()


def approx(a):
    return np.array([complex(v) for v in a.entries()]).reshape(a.nrows(),a.ncols())


def spectral(a):
    ap=approx(a)
    values,vectors=np.linalg.eigh((ap+ap.conj().T)/2)
    q=as_ball(vectors)
    d=diag(values)
    orth=norm(adjoint(q)*q-identity(len(values)))
    assert orth<1, repr(orth)
    err=norm(a*q-q*d)+orth*(norm(a)+max(abs(values)))
    err=err.upper()
    negative=sum(arb(float(v))+err<0 for v in values)
    positive=sum(arb(float(v))-err>0 for v in values)
    return values,err,negative,positive


def transformed(g,data,gamma,multiply_g):
    pad=g.symbol.shape[0]
    columns=[]
    for source in range(data.shape[1]):
        vec=data[:,source].reshape(g.side,g.side)
        grid=[[ball(vec[y,x]) if y<g.side and x<g.side else acb(0)
               for x in range(pad)] for y in range(pad)]
        first=[acb.dft(row) for row in grid]
        second=list(zip(*(acb.dft(list(col)) for col in zip(*first))))
        weights=[[v*gamma[y][x]*(ball(g.symbol[y,x]) if multiply_g else acb(1))
                  for x,v in enumerate(row)] for y,row in enumerate(second)]
        spatial=inverse_dft_table(weights)
        columns.append([v for y,row in enumerate(spatial) if y<g.side
                        for x,v in enumerate(row) if x<g.side])
    return acb_mat(g.n,data.shape[1],[z for row in zip(*columns) for z in row])


def block_rows(a,rows):
    return acb_mat(len(rows),a.ncols(),[a[int(i),j] for i in rows for j in range(a.ncols())])


def stack(a,b):
    assert a.ncols()==b.ncols()
    return acb_mat(a.nrows()+b.nrows(),a.ncols(),a.entries()+b.entries())


def block_cols(a,b):
    return adjoint(stack(adjoint(a),adjoint(b)))


class ReplayGeometry:
    def __init__(self,g,b):
        self.g,self.b=g,b
        atom=[[ball(v) for v in row] for row in g.symbol]
        self.gamma=[[1/(1-arb(b)*v.real) for v in row] for row in atom]
        tables=[
            [[acb(v) for v in row] for row in self.gamma],
            [[v*f.conjugate() for v,f in zip(row,coef)] for row,coef in zip(self.gamma,atom)],
            [[v*(f.real**2+f.imag**2) for v,f in zip(row,coef)] for row,coef in zip(self.gamma,atom)]
        ]
        self.blocks=[list(inverse_dft_table(v)) for v in tables]
        self.S=as_ball(g.S)
        self.B=as_ball(g.B)
        self.h=as_ball((b/2)*g.B)
        x0=transformed(g,(b/2)*g.B,self.gamma,False)
        xg=transformed(g,(b/2)*g.B,self.gamma,True)
        s0=transformed(g,g.S.conj().T,self.gamma,False)
        sg=transformed(g,g.S.conj().T,self.gamma,True)
        self.x0,self.xe=x0,xg+self.B
        self.s0,self.sg=s0,sg
        self.P0,self.y0,self.t0=self.S*s0,self.S*x0,adjoint(self.h)*x0

    def covariance(self,cells):
        g=self.g
        k=len(cells)
        out=[]
        for i,j in np.ndindex(2*k,2*k):
            a,b=divmod(i,k),divmod(j,k)
            dif=(g.xy[cells[a[1]]]-g.xy[cells[b[1]]])%g.symbol.shape[0]
            lag=dif[1],dif[0]
            if a[0]==1 and b[0]==0:
                rev=(-dif[1])%g.symbol.shape[0],(-dif[0])%g.symbol.shape[0]
                val=self.blocks[1][rev[0]][rev[1]].conjugate()
            else:
                val=self.blocks[a[0]+b[0]][lag[0]][lag[1]]
            out.append(val)
        raw=acb_mat(2*k,2*k,out)
        return (raw+adjoint(raw))*arb(.5)

    def verify(self,center,shape,lo,hi,broad,w,tau,Y,r,diagnostic_only=False):
        start=clock()
        off,p=pattern(self.g.side,shape)
        cells,p=place(self.g.side,center,off,p)
        ell=p*lo if not broad else np.where(p==1.,lo,0.)
        upper=p*hi if not broad else np.full(len(p),hi)
        k=len(p)
        n_dpos=0
        updates=[]
        for l,u in zip(ell,upper):
            ww=arb(w)
            aa=ww-1
            ab=-(arb(float(l))+arb(float(u)))*ww/2+arb(self.b)/2
            bb=ww*arb(float(l))*arb(float(u))
            determinant=aa*bb-ab*ab
            assert determinant<0 or determinant>0
            n_dpos+=1 if determinant<0 else (2 if aa+bb>0 else 0)
            updates.append((bb/determinant,-ab/determinant,aa/determinant))
        inverse=acb_mat(2*k,2*k,[
            updates[i%k][i//k+j//k] if (i%k)==(j%k) else 0
            for i,j in np.ndindex(2*k,2*k)])
        A=self.covariance(cells)+inverse
        if diagnostic_only:
            return
        raw, err, neg, pos=spectral(A)
        assert neg+pos==len(raw)
        nu=pos-n_dpos
        C=adjoint(stack(block_rows(self.s0,cells),block_rows(self.sg,cells)))
        V=stack(block_rows(self.x0,cells),block_rows(self.xe,cells))
        sol=A.solve(block_cols(adjoint(C),V))
        s=C.nrows()
        XC=acb_mat(2*k,s,[sol[i,j] for i,j in np.ndindex(2*k,s)])
        XV=acb_mat(2*k,V.ncols(),[sol[i,j+s] for i,j in np.ndindex(2*k,V.ncols())])
        P=self.P0-C*XC
        center_data=self.y0-C*XV
        theta=self.t0-adjoint(V)*XV
        P=(P+adjoint(P))*arb(.5)
        theta=(theta+adjoint(theta))*arb(.5)
        port=identity(s)+P/arb(tau)
        pv,pe,pn,pp=spectral(port)
        assert pn+pp==len(pv)
        state=nu-pn
        assert 0<=state<Y.shape[1]
        difference=as_ball(Y)-center_data
        H=adjoint(difference)*port.solve(difference)-arb(tau)*theta
        H=(H+adjoint(H))*arb(.5)
        values,error,_,_=spectral(H)
        floors=[arb(float(v))-error for v in values[:Y.shape[1]-state]]
        gap=sum((v.lower() for v in floors if v>0),arb(0)).sqrt()-arb(r)
        assert gap>0, (center,gap,values,error)
        return dict(centre=center,physical_cells=k,negative_index=int(state),
                    lower_margin=float(np.nextafter(float(gap.lower()),-np.inf)),
                    w=w,tau=tau,trace_eigen_error=float(err),
                    receiver_eigen_error=float(pe),source_eigen_error=float(error),
                    **elapsed(start))


def run(side,shape,obs,mode,samples,interior=False):
    ctx.prec=128
    start=clock()
    g=Geometry(geom_path(side))
    b=FREEZE["model"]["background_upper"]
    glob=f"ipcl-fov-f1__{obs}__weak__7tx9rx__n03"
    geometry=ReplayGeometry(g,b)
    args=(3.,5.,True) if mode=="broad" else (3.75,4.25,False)
    preparation=elapsed(start)
    saved=np.load(HERE/"results"/f"{glob}__b{b:g}__{mode}_sector.npz",allow_pickle=False)
    observe=np.load(HERE/"public"/f"{glob}.npz",allow_pickle=False)
    raw=saved["margin_sector"]
    allowed=np.flatnonzero(raw>0)
    if side==96:
        # A clipped case keeps the interval calculation short while
        # still evaluating native full-domain background covariance.
        off,p=pattern(side,shape)
        support=[len(place(side,c,off,p)[0]) for c in g.xy[allowed]]
        if interior:
            allowed=allowed[np.asarray(support)==len(p)]
        else:
            allowed=allowed[np.argsort(support)[:max(20,samples)]]
    if samples==1:
        ranks=[allowed[np.argmax(raw[allowed])]]
    else:
        ordered=allowed[np.argsort(raw[allowed])]
        ranks=[ordered[int(q)] for q in np.linspace(0,len(ordered)-1,samples)]
    records=[]
    for node in ranks:
        centre=g.xy[node].tolist()
        proposed=float(raw[node])
        checked=geometry.verify(centre,shape,*args,float(saved["weight"][node]),
                                float(saved["sector_scale"][node]),observe["Y"],float(observe["r"]))
        checked["saved_lower_margin"]=proposed
        records.append(checked)
        print("Arb",side,obs,mode,checked,flush=True)
    suffix="_interior" if interior else ""
    save_json(HERE/"results"/f"arb_{obs}_{mode}{suffix}.json",
              dict(precision=ctx.prec,obs_id=glob,records=records,preparation=preparation,
                   total=elapsed(start)))


if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("side",type=int)
    a.add_argument("shape")
    a.add_argument("obs")
    a.add_argument("mode")
    a.add_argument("--samples",type=int,default=1)
    a.add_argument("--interior",action="store_true")
    args=a.parse_args()
    with threadpool_limits(limits=1,user_api="blas"):
        run(args.side,args.shape,args.obs,args.mode,args.samples,args.interior)
