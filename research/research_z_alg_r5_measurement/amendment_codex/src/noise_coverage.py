"""Uniform interval coefficient/root coverage, without GT or new measurements."""
from complete_quotient import *

def run_noise(case='T0',arm='noisy',eta=1.2e-4):
 mp.iv.prec=128;begin=time.process_time();book=Ledger('parametric_noise_coverage');g=geometry(case);pos,Q=dictionary(case,g);p=Physics(g,pos,book);y=np.load(BASE/'inputs'/f'{case}.npz')['old_'+arm]
 center=chart(p.H,p.B[:,:p.s],p.R,y);uncertain=chart(p.H,p.B[:,:p.s],p.R,y,data_radius=eta);polys,basis,leading=interval_minors(uncertain['coeff']);lead=invertibility([r[len(basis):] for r in polys]);roots,_=enumerate_midpoint(center['coeff']);certs=[];classes=[]
 for x,physical in roots:
  cr=root_certificate(uncertain['coeff'],x,physical,radii=(1e-7,1e-6,1e-5,1e-4,1e-3,.003,.01,.03,.1,.3,1.));certs.append(cr)
  classes.append(classify(uncertain,cr,p.H,p.R,BOX))
 allcert=all(c['status']=='CERTIFIED_COMPLEX_ROOT' for c in certs);distinct=all(disjoint(a['X'][:p.s],b['X'][:p.s]) for a,b in itertools.combinations(certs,2)) if allcert else False
 covered=bool(lead['certified'] and allcert and distinct)
 result=dict(case=case,arm=arm,eta_componentwise=eta,precision_bits=128,leading=lead,all_parametric_complex_roots_certified=allcert,distinct_c_boxes=distinct,noise_feasible_root_coverage=covered,certified_root_count=sum(c['status']=='CERTIFIED_COMPLEX_ROOT' for c in certs),total_roots=len(certs),classifications=classes,certificates=[{k:v for k,v in cr.items() if k!='X'} for cr in certs],scope='Uniform old-data real/imag interval box of nominal rounded discrete model; not a physical model-error calibration.',CPU_seconds=time.process_time()-begin,physics_receipt=book.receipt())
 save(OUT/'results'/f'noise_{case}_{arm}.json',result);print({k:v for k,v in result.items() if k not in ('certificates','classifications','physics_receipt')},flush=True)
 return result
if __name__=='__main__':run_noise(*(sys.argv[1:] or ['T0','noisy']))
