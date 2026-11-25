"""Rebuild selected weak-background PSD proposals from the original Fourier symbol in Arb."""
import argparse
import json
import sys

from flint import acb, acb_mat, arb, ctx
import numpy as np
from scipy.linalg import eigh

from search import HERE, MAPS, Geometry, box, clock, elapsed, geom_path, read_spec, save_json, threadpool_limits
from operator_bounds import ball as scalar_ball, inverse_dft_table
sys.path.insert(0,str(HERE.parent))
from interval_backend import adj,ball,blocks,diagonal,floor,sub,verify_signs


def transpose(table):
    return [list(col) for col in zip(*table,strict=True)]


def fft2(table):
    return transpose([acb.dft(col) for col in transpose([acb.dft(row) for row in table])])


class ExactPadded:
    def __init__(self,g,b):
        self.g=g
        self.b=arb(float(b))
        self.n=g.side
        self.p=g.symbol.shape[0]
        self.symbol=[[scalar_ball(x) for x in row] for row in g.symbol]
        self.gamma=[[1/(1-self.b*entry.real) for entry in row] for row in self.symbol]
        if not all(u > 0 for row in self.gamma for u in row):
            raise ArithmeticError('unresolved positive padded background eigensymbol')
        weighted = [[d*s for d,s in zip(drow,srow,strict=True)] for drow,srow in zip(self.gamma,self.symbol,strict=True)]
        lk=[[[acb(q) for q in row] for row in self.gamma],
            [[q.conjugate() for q in row] for row in weighted],
            [[q*s.conjugate() for q,s in zip(wrow,srow,strict=True)] for wrow,srow in zip(weighted,self.symbol,strict=True)]]
        self.lag=[[tuple(row) for row in inverse_dft_table(q)] for q in lk]
        self.sensor=ball(g.S)
        self.b0=self.b/2*ball(g.B)
        rhs=np.column_stack((g.S.conj().T,b/2*g.B))
        base,green=self.apply(rhs,[self.gamma,weighted])
        pr,ns=g.S.shape[0],g.B.shape[1]
        self.sensorcols=range(pr)
        self.sourcecols=range(pr,pr+ns)
        self.bs=base; self.gs=green
        self.p0=self.sensor*sub(base,range(g.n),self.sensorcols)
        inc=sub(base,range(g.n),self.sourcecols)
        self.y0=self.sensor*inc
        self.t0=adj(self.b0)*inc

    def apply(self,rhs,products):
        results=[[] for _ in products]
        for col in rhs.T:
            cropped=col.reshape(self.n,self.n)
            lattice=[[scalar_ball(v) for v in row] + [acb(0)]*(self.p-self.n) for row in cropped]
            lattice.extend([[acb(0)]*self.p for _ in range(self.p-self.n)])
            transform=fft2(lattice)
            for k,product in enumerate(products):
                table=[[x*y for x,y in zip(fr,pr,strict=True)] for fr,pr in zip(transform,product,strict=True)]
                result=[tuple(row) for row in inverse_dft_table(table)]
                results[k].append([v for row in result[:self.n] for v in row[:self.n]])
        return tuple(acb_mat(rhs.shape[1],self.g.n,[e for col in rows for e in col]).transpose() for rows in results)

    def covariance(self,ids):
        xy=self.g.xy[ids]
        shifts=((xy[:,None,:]-xy[None,:,:])%self.p).reshape(-1,2)
        lags=[acb_mat(len(ids),len(ids),[li[int(y)][int(x)] for x,y in shifts]) for li in self.lag]
        return blocks([[lags[0],lags[1]],[adj(lags[1]),lags[2]]])

    def local(self,node,weight):
        ids,profile,lo,hi=box(node)
        nd=len(ids)
        ratio=arb(float(weight))
        a=ratio-1
        dinv=acb_mat(2*nd,2*nd)
        n_pos=0
        for i,(u,v) in enumerate(zip(lo*profile,hi*profile,strict=True)):
            u,v=arb(float(u)),arb(float(v))
            t=-ratio*(u+v)/2+self.b/2
            f=ratio*u*v
            d=a*f-t*t
            if d < 0:
                n_pos+=1
            elif d > 0 and a+f > 0:
                n_pos+=2
            elif not (d > 0 and a+f < 0):
                raise ArithmeticError('unresolved exact physical update rank')
            dinv[i,i]=f/d
            dinv[i,i+nd]=dinv[i+nd,i]=-t/d
            dinv[i+nd,i+nd]=a/d
        A=self.covariance(ids)+dinv
        full=verify_signs(A)
        if full['unresolved']:
            raise ArithmeticError('unresolved exact local trace signs')
        nu=2*nd-full['negative_at_least']-n_pos
        C=adj(blocks([[sub(self.bs,ids,self.sensorcols)],[sub(self.gs,ids,self.sensorcols)]]))
        eb=ball(self.g.B[ids])
        V=blocks([[sub(self.bs,ids,self.sourcecols)],[sub(self.gs,ids,self.sourcecols)+eb]])
        both=A.solve(blocks([[adj(C),V]]))
        k=C.nrows(); s=V.ncols(); allstate=range(2*nd)
        cs=sub(both,allstate,range(k)); vs=sub(both,allstate,range(k,k+s))
        P=self.p0-C*cs
        yc=self.y0-C*vs
        theta=self.t0-adj(V)*vs
        return P,yc,theta,nu,full


def cutoff(ready,y,r,scale):
    P,yc,theta,nu,phys=ready
    mu=arb(float(scale))
    port=diagonal([1]*P.nrows())+1/mu*P
    counts=verify_signs(port)
    output=dict(local_state_negative_index=nu,strict_port_counts=counts,local_trace_sign_proof=phys)
    if counts['unresolved']:
        return output|{'status':'unresolved_port_signs'}
    actual=nu-counts['negative_at_least']
    output['verified_state_negative_index']=actual
    if actual!=0:
        return output|{'status':'verified_nonzero_state_negative_index'}
    err=ball(y)-yc
    h=adj(err)*port.solve(err)-mu*theta
    s=h.nrows()
    hmid=np.array([complex(x) for x in h.entries()]).reshape(s,s)
    val,Q=eigh((hmid+hmid.conjugate().T)/2,check_finite=False)
    raw=Q[:,val>0]
    if not raw.shape[1]:
        return output|{'status':'no_positive_tail_proposed','verified_margin':-r}
    q=ball(raw)
    gram=adj(q)*q
    denom=-floor(-gram)
    trace=sum(((adj(q)*(h*q))[i,i].real for i in range(q.ncols())),arb(0))
    cert=(trace.lower()/denom.upper())
    margin=arb(0).max(cert).sqrt()-arb(float(r)) if cert>0 else -arb(float(r))
    edge=float(np.nextafter(float(margin.lower()),-np.inf))
    return output|{'status':'EXCLUDED' if margin>0 else 'UNRESOLVED','verified_margin':edge,
                   'raw_arb_midpoint_margin':float(np.sqrt(np.maximum(val,0).sum())-r),
                   'source_positive_tail_frame_columns':q.ncols()}


def run(side, refined=False):
    spec=read_spec()
    if refined and side!=32:
        raise ValueError('only the 32-square broad cases have refined ratio searches')
    result='refined_broad_32.json' if refined else f'result_{side}.json'
    label=f'exact_replay_{"refined_" if refined else ""}{side}.json'
    rows=json.loads((HERE/result).read_text())['nodes']
    g=Geometry(geom_path(side)); started=clock()
    exact=ExactPadded(g,spec['b']); offline=elapsed(started)
    proofs=[]
    for row in rows:
        proposals=[row['positive_continuous_raw_optimum'],row['positive_continuous_verified_best']]
        with np.load(MAPS/'public'/f'{row["oid"]}.npz',allow_pickle=False) as record:
            y,r=record['Y'],float(record['r'])
        prepared={}
        done=set()
        receipt=[]
        for proposal in proposals:
            if proposal is None:
                continue
            w,m=proposal['ratio'],proposal['scale']
            if (w,m) in done:
                continue
            nodeclock=clock(); done.add((w,m))
            try:
                if w not in prepared:
                    prepared[w]=exact.local(row,w)
                proof=cutoff(prepared[w],y,r,m)
            except (ArithmeticError,ValueError,ZeroDivisionError) as exc:
                proof={'status':'UNRESOLVED','reason':str(exc)}
            receipt.append(dict(ratio=w,scale=m,numeric_raw_proposal=proposal['raw_positive_part_margin'],
                                preliminary_verified_margin=proposal['verified_margin'],**proof,cost=elapsed(nodeclock)))
        proofs.append({'oid':row['oid'],'index':row['index'],'family':row['family'],'proposals':receipt})
        save_json(HERE/label,{'precision_bits':ctx.prec,'geometry_arb_preparation_cost':offline,
                                                  'scope':'exact padded Fourier operator contains all zero-extended physical currents',
                                                  'rows':proofs,'elapsed':elapsed(started)})
        print(side,row['index'],receipt,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('side',type=int,choices=(32,96))
    p.add_argument('--refined',action='store_true'); a=p.parse_args()
    ctx.prec,ctx.threads=192,1
    with threadpool_limits(limits=1):run(a.side,a.refined)
