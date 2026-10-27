"""AST-extracted decoder-only QA; never imports the paid physics adapter."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
import ast,json,time,resource,functools
from pathlib import Path
import numpy as np
from scipy.special import ndtr,ndtri,betaincinv
from scipy.stats import norm
from scipy.integrate import quad
from numpy.polynomial.legendre import leggauss
OUT=Path(__file__).resolve().parent;SRC=OUT.parents[1]/'track_a/stratified.py'
def saturation(s):
 t=np.clip((np.asarray(s)+1)/2,0,1);return t*t*t*(10+t*(-15+6*t))
# No source imports, Book/Target construction, forward physics, or GT access.
names={'LO','HI','Z0','A','PH','ULO','UHI','RLO','RHI','ILO','IHI'};funs={'gl','inverse_fraction','fibre','quotient_representative','source_nodes','quotient_nodes'}
node=ast.parse(SRC.read_text());selected=[]
for x in node.body:
 if isinstance(x,ast.Assign) and {n.id for target in x.targets for n in ast.walk(target) if isinstance(n,ast.Name)} <= names:selected.append(x)
 elif isinstance(x,ast.FunctionDef) and x.name in funs:selected.append(x)
ns=dict(np=np,norm=norm,ndtr=ndtr,betaincinv=betaincinv,leggauss=functools.lru_cache(None)(leggauss),saturation=saturation,SIGMA=.00119211798514519)
exec(compile(ast.Module(body=selected,type_ignores=[]),str(SRC),'exec'),ns)
def summarize(nodes):
 V=nodes['V'];w=nodes['prior_weights'];return dict(nodes=len(V),mass=float(w.sum()),branch_mass=[float(w[nodes['layer']==k].sum()) for k in range(3)],phase_gt_half=float(nodes['phase_numerator'].sum()),mean_V=[float((w*V).sum().real),float((w*V).sum().imag)],second_abs_moment=float(np.sum(w*abs(V)**2)))
def main():
 start=time.perf_counter();A=ns['A'];Z0=ns['Z0'];PH=ns['PH'];eq2=quad(lambda z:float(saturation(1.3*z))**2*norm.pdf(z)/Z0,-4,4,points=[-A,A],epsabs=1e-12)[0];rlo,rhi,ilo,ihi=[ns[k] for k in ['RLO','RHI','ILO','IHI']];second=eq2*((rlo*rlo+rlo*rhi+rhi*rhi)+(ilo*ilo+ilo*ihi+ihi*ihi))/3
 analytic=dict(mass=1,branch_mass=[PH,1-2*PH,PH],phase_gt_half=.5,mean_V=[.65,.0525],second_abs_moment=second,normalizer_Z0=Z0,uniform_area_width=(rhi-rlo)*(ihi-ilo))
 source=[]
 for d,nz,nu in [(2,12,12),(2,24,20),(2,40,28),(3,10,8),(3,18,10),(3,40,16)]:source.append(dict(dimension=d,nz=nz,nu=nu,**summarize(ns['source_nodes'](d,nz,nu))))
 quotient=[]
 for nr,nt,nf in [(8,8,32),(16,8,48),(28,12,64),(48,24,96),(80,36,128)]:
  v=summarize(ns['quotient_nodes'](nr,nt,nf));quotient.append(dict(nr=nr,nt=nt,nf=nf,**v));print('quotient',nr,v['mass'],v['branch_mass'],flush=True)
 n=256000;rng=np.random.default_rng(202610270628);u=rng.uniform(ns['ULO'],ns['UHI'],size=(n,3));z=ndtri(u[:,0]);Q=saturation(1.3*z);D=.8+u[:,1]+1j*(.03+.15*u[:,2]);V=Q*D
 mc=dict(n=n,branch_mass=[float(np.mean(z<=-A)),float(np.mean((z>-A)&(z<A))),float(np.mean(z>=A))],phase_gt_half=float(np.mean(Q>.5)),mean_V=[float(V.real.mean()),float(V.imag.mean())],mean_V_standard_error=[float(V.real.std()/np.sqrt(n)),float(V.imag.std()/np.sqrt(n))],second_abs_moment=float(np.mean(abs(V)**2)),second_abs_moment_standard_error=float(np.std(abs(V)**2)/np.sqrt(n)))
 # Pointwise fibre independent q-variable integration away from endpoint.
 fibre=[]
 for R,I in [(.001,.0001),(.02,.002),(.3,.03),(1.2,.0975),(1.7,.17)]:
  d,e=ns['fibre'](R,I,128);lo=max(R/rhi,I/ihi);hi=min(R/rlo,I/ilo,1.);a=float(ns['inverse_fraction'](lo));b=float(ns['inverse_fraction'](hi));f=lambda z:norm.pdf(z)/Z0/float(saturation(1.3*z))**2/((rhi-rlo)*(ihi-ilo));ref=quad(f,a,b,epsabs=1e-9,epsrel=1e-10,points=[0.] if a<0<b else None)[0];er=quad(f,max(a,0),b,epsabs=1e-9,epsrel=1e-10)[0] if b>0 else 0.;fibre.append(dict(R=R,I=I,density=d,reference=ref,relative_error=(d-ref)/ref,phase_numerator=e,phase_reference=er))
 # Representative decoder replay, no Maxwell forward.
 positive=V[Q>0][::1000];errs=[]
 for v in positive:
  q=ns['quotient_representative'](v.real,v.imag);vv=saturation(1.3*q[0])*complex(.8+q[64],.03+.15*q[65]);errs.append(abs(v-vv))
 result=dict(status='DECODER_ALGEBRA_QA_COMPLETE_NO_REFERENCE_CERTIFICATE',source=str(SRC),analytic=analytic,source_quadrature=source,quotient_quadrature=quotient,monte_carlo=mc,fibre_pointwise_checks=fibre,representative_max_abs_error=max(errs),representative_count=len(errs),CPU_threads=1,physics_RHS=0,GT_read=False,Book_or_Target_constructed=False,wall_seconds=time.perf_counter()-start,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 (OUT/'QA.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
