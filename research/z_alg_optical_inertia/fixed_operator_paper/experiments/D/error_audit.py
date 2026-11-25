"""Uniform calibration inequality, two outward cell discretizations, and measured margins."""
import csv,json,sys,time
from pathlib import Path
import numpy as np
from flint import acb,arb,ctx,acb_mat
p=Path(__file__).resolve().parent
study=p.parents[2]
sys.path.insert(0,str(p.parents[1]/'src'))
from fixed_trace import adj,diag
ctx.prec=128

def bounds(b):
    return [float(np.nextafter(float(b.lower()),-np.inf)),float(np.nextafter(float(b.upper()),np.inf))]
def frob(m):return sum((abs(z)**2 for z in m.entries()),arb(0)).sqrt()
def octagon(side):
    # Physical square of dielectric contrast 2 on [-1/8,1/8]^2 within a unit field of view.
    # At both spacings the material boundary coincides with cell edges; eliminate exact vacuum.
    h=arb(1)/side;k=2*arb.pi();u=k*h/arb.pi().sqrt()
    pts=[((arb(i)+arb(.5))/side-arb(.5),(arb(j)+arb(.5))/side-arb(.5)) for i in range(3*side//8,5*side//8) for j in range(3*side//8,5*side//8)]
    def h1(n,x):return acb(x.bessel_j(n),x.bessel_y(n))
    c=acb(0,1)*arb.pi()*u*u.bessel_j(1)/2
    selfterm=acb(0,1)*arb.pi()*u*h1(1,u)/2-1
    def element(a,b):return c*h1(0,k*((a[0]-b[0])**2+(a[1]-b[1])**2).sqrt())
    G=acb_mat([[selfterm if i==j else element(a,b) for j,b in enumerate(pts)] for i,a in enumerate(pts)])
    rec=[(arb(2),arb(0)),(arb(0),arb(2)),(-arb(2),arb(0)),(arb(0),-arb(2))]
    S=acb_mat([[element(r,a) for a in pts] for r in rec])
    B=acb_mat([[acb(0,k*x).exp(),acb(0,k*y).exp()] for x,y in pts])
    return S*(diag([1]*len(pts))-2*G).solve(2*B)
def main():
    started=time.perf_counter()
    a,b=octagon(8),octagon(16)
    discr=frob(a-b)/frob(b)
    # exact common gain bound example. Declared norm bound and residual are measured independently.
    gain=acb(0,arb(1)/32);eps=arb(1)/32
    additive=acb_mat(4,2);additive[0,0]=arb(1)/2048
    y=(1+gain)*b+additive;known_noise=frob(additive)
    cap=(eps*frob(y)+known_noise)/(1-eps)
    assert frob(b-y).upper()<cap.lower()
    source=study/'synthetic/evaluation'
    sampled=[]
    for f in sorted(source.glob('*__*__n02.npz')):
        with np.load(f) as z:
            fine=z['Y_fine_clean'];coarse=z['Y_coarse_proxy']
            sampled.append(dict(id=f.stem,rel=float(np.linalg.norm(coarse-fine)/np.linalg.norm(fine))))
    lab=list(csv.DictReader((study/'ip_guided_v4/routing/measured_port_cases.csv').open()))
    margin=json.loads((study/'validation/transfer/summary.json').read_text())['passive_sector']
    group={}
    for name,ident in [('Fresnel','fresnel'),('Manitoba','manitoba')]:
        rows=[r for r in lab if r['id'].startswith(ident)]
        group[name]=dict(views=len(rows),objects=len(set(r['object_id'] for r in rows)),born_above_proxy=sum(r['historical_exceeds_both_error_proxies']=='True' for r in rows),
                         held_rx=[min(float(r['held_rx_error']) for r in rows),max(float(r['held_rx_error']) for r in rows)],
                         held_tx=[min(float(r['held_tx_error']) for r in rows),max(float(r['held_tx_error']) for r in rows)])
    out=dict(seconds=time.perf_counter()-started,bits=128,interval_cell_kernel=dict(sides=[8,16],active=[4,16],relative_difference=bounds(discr),scope='two exact Richmond discretizations, uniform continuum quadrature error unknown'),
        common_gain=dict(epsilon=str(eps),noise_norm=bounds(known_noise),observed_error=bounds(frob(b-y)),upper=bounds(cap)),
        sample_discrepancy=sampled,laboratory=group,manitoba_cap_one=margin['1.0'],unconditional_loaded_error_radius=None)
    (p/'audit.json').write_text(json.dumps(out,indent=2))
    # Channel summaries only; facility raw fields remain at the original repositories.
    with (p/'ports.csv').open('w') as f:
        keys=('id','object_id','frequency_GHz','held_rx_error','held_tx_error','historical_born_lower_relative','historical_exceeds_both_error_proxies')
        w=csv.DictWriter(f,fieldnames=keys,extrasaction='ignore');w.writeheader();w.writerows(lab)
    print(json.dumps({k:v for k,v in out.items() if k not in ['sample_discrepancy']},indent=2))
if __name__=='__main__':main()
