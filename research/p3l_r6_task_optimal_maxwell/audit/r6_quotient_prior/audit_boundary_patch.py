"""Revision2 audit. Reuses previous256k MC receipt; no more random draws."""
import importlib.util,json,time,resource
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('prior_audit',HERE/'audit.py');a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
start=time.perf_counter();source=a.SRC.read_text();(HERE/'SOURCE_USED_BOUNDARY_PATCH.py').write_text(source);rows=[]
for d,nz,nu in [(2,12,12),(2,24,16),(2,32,24),(3,10,8),(3,16,9)]:rows.append(dict(dimension=d,nz=nz,nu=nu,**a.summarize(a.ns['source_nodes'](d,nz,nu))))
qr=[]
for nr,nt,nf in [(12,6,48),(20,8,64),(32,16,96),(64,32,128)]:
 q=a.summarize(a.ns['quotient_nodes'](nr,nt,nf));qr.append(dict(nr=nr,nt=nt,nf=nf,**q));print('quotient',nr,q['mass'],q['branch_mass'],q['phase_gt_half'],flush=True)
old=json.loads((HERE/'QA_initial_unsplit.json').read_text());r=dict(status='BOUNDARY_PATCH_PRIOR_QA_COMPLETE_NO_POSTERIOR_REFERENCE_CERTIFICATE',revision='separate_pure_tensorGL_and_radial_clipping_panels_and_z0_event_split',source_snapshot=str(HERE/'SOURCE_USED_BOUNDARY_PATCH.py'),analytic=old['analytic'],source_quadrature=rows,quotient_quadrature=qr,monte_carlo_reused_from='QA_initial_unsplit.json',additional_MC_samples=0,total_MC_samples=256000,physics_RHS=0,GT_read=False,Book_or_Target_constructed=False,CPU_threads=1,wall_seconds=time.perf_counter()-start,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
(HERE/'QA_boundary_patch.json').write_text(json.dumps(r,indent=2)+'\n')
