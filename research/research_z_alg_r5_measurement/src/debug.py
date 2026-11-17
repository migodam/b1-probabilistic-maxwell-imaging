from common import *
book=Ledger('historical_debug_and_formula_tests');g=base_geometry();positions,Q=dictionary('DEBUG',g);p=Physics(g,positions,book)
old=np.eye(p.L)[:,:p.s];target=np.load(R2/'experiments/general6_public/GENERAL6_00.npz')['y'];ker=p.old_kernel();system=GeneralBilinearDualSystem(ker,target)
r=json.loads((R2/'results/general6_GENERAL6_00_general_quotient20.json').read_text());roots=[np.array(v['m']) for v in r['independently_validated']]
ws=[p.dual(m) for m in roots];cs=[(w-system.w0)@np.linalg.pinv(system.N) for w in ws]
menu=p.menu(Q,roots,roots[0]);rows=[]
for item in menu:
 a=item['a'];b=p.B@a;pred=[(system.w0+c@system.N)@b for c in cs];full=[p.response(m,a,independent=True)[0] for m in roots]
 rows.append(dict(label=item['label'],a=a,pred=pred,full=full,gap=abs(pred[0]-pred[1]),dual_error=max(abs(np.array(pred)-full)),energy=item['energy'],feed_norm=item['feed_norm']))
# Nonzero host, coordinate selection, independent full solve.
m0=np.full(6,1.2);host=ReducedKernel(p.model,m0,active=[0,2,3,5]);z=np.array([.7,2.,-.2,1.1]);m=m0+host.U@z
red=host(z);full=p.model.forward(m,independent=True)['y']
# Actual analytic J versus independent centered finite difference.
y,J=p.response(roots[0],old,True);fd=[]
for i in range(6):
 h=1e-5;d=np.eye(6)[i]*h;fd.append(pack((p.response(roots[0]+d,old,independent=True)-p.response(roots[0]-d,old,independent=True))/(2*h)))
jerr=np.linalg.norm(J-np.array(fd).T)/np.linalg.norm(J)
# Direct completion at deterministic QR selected primitives.
from scipy.linalg import qr
remain=list(range(p.s,p.L));inds=list(range(p.s))
while len(inds)<p.n:
 scores=[]
 for j in remain:scores.append(np.linalg.svd(p.B[:,inds+[j]],compute_uv=False)[-1])
 j=remain.pop(int(np.argmax(scores)));inds.append(j)
y=p.response(roots[0],np.eye(p.L)[:,inds],independent=True);w=np.linalg.solve(p.B[:,inds].T,y);mrec=w/(p.R.reshape(-1)+w@p.H)
# Perturbation formula explicitly changes w0 with t.
dt=np.array([[1e-5+2e-5j,-2e-5+1e-5j,3e-5-1e-5j]])
sys2=GeneralBilinearDualSystem(ker,target+dt)
checks=dict(host_relative=np.linalg.norm(red-full)/np.linalg.norm(full),jacobian_relative=jerr,NB=np.linalg.norm(system.N@ker.q),wB=max(np.linalg.norm(w@ker.q-target) for w in ws),dual_max_error=max(r['dual_error'] for r in rows),old_span_gap=max(r['gap'] for r in rows if r['label'].startswith('old')),largest_accessible_gap=max(r['gap'] for r in rows),direct_material_error=np.linalg.norm(mrec-roots[0]),noise_center_error=np.linalg.norm(sys2.w0-system.w0-(dt@np.linalg.pinv(ker.q)).reshape(-1)))
checks['passed']=checks['host_relative']<1e-12 and jerr<1e-6 and checks['dual_max_error']<1e-11 and checks['old_span_gap']<1e-10 and checks['direct_material_error']<1e-8 and checks['noise_center_error']<1e-12
save(ROOT/'results/debug.json',dict(checks=checks,sources=rows,roots=roots,scope='EXPOSED_HISTORICAL_DEBUG_NOT_NEW_TARGET',receipt=book.receipt()))
print(json.dumps(serial(checks),indent=2))
assert checks['passed']
