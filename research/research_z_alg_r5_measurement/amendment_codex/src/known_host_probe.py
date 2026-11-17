"""One frozen exposed applicability probe: n=12, k=6, known host."""
from complete_quotient import *
from constrained_source import constrained_match

def main():
 begin=time.process_time();mp.iv.prec=128;book=Ledger('H12_known_host');g=geometry('T0');pos,Q=dictionary('T0',g);theta=np.arange(6)*2*np.pi/6+.23;hostpos=.55*np.c_[np.cos(theta),np.sin(theta)];centers=np.vstack([g['centers'],hostpos]);m0=np.r_[np.zeros(6),np.full(6,2.)];z=np.array([19.49822831018447,9.40854042875942,3.8432927569521156,2.490230030345602,19.6893640302288,26.11329805811731])
 save(OUT/'results/H12_contract.json',dict(exposure='EXPOSED_APPLICABILITY_PROBE_NOT_BLIND',n=12,k=6,s=3,centers=centers,known_host=m0,active=list(range(6)),declared_active_material=z,noise=0,host_mismatch='NOT_RUN',source_rule='Constrained pair direction if >=2 certified legal roots; otherwise no ambiguity experiment. No target replacement.',CPU_limit_seconds=600))
 model=Model(book,centers=centers,sources=pos,receivers=g['receivers'],radius=g['radius'],k=g['k']);kernel=ReducedKernel(model,m0,active=np.arange(6));full=model.forward(m0+kernel.U@z,independent=True);yf=full['y'];yr=kernel.evaluate(z);t=(yf-kernel.y0).reshape(-1)[:3];H=kernel.K;B=kernel.q;R=kernel.L;ch=chart(H,B[:,:3],R,t);polys,basis,leading=interval_minors(ch['coeff']);lead=invertibility([r[len(basis):] for r in polys]);roots,diag=enumerate_midpoint(ch['coeff']);certs=[];classes=[];legal=[]
 for x,phys in roots:
  cr=root_certificate(ch['coeff'],x,phys);cl=classify(ch,cr,H,R,BOX);certs.append(cr);classes.append(cl)
  if cl['status']=='CERTIFIED_LEGAL_REAL_MATERIAL':legal.append(np.array([(v['lo']+v['hi'])/2 for v in cl['material']]))
 allcert=all(c['status']=='CERTIFIED_COMPLEX_ROOT' for c in certs);distinct=all(disjoint(a['X'][:3],b['X'][:3]) for a,b in itertools.combinations(certs,2)) if allcert else False;complete=lead['certified'] and allcert and distinct and len(certs)==20 and all(c['status'] in ('CERTIFIED_LEGAL_REAL_MATERIAL','EXCLUDED_NONCONJUGATE','EXCLUDED_OUTSIDE_BOX') for c in classes);pred=[];ws=[]
 for m in legal:
  with book.op('REDUCED_DUAL',reduced_LU=1,reduced_complex_RHS=1):w=np.linalg.solve((np.eye(6)-m[:,None]*H).T,R.T).reshape(-1)*m
  ws.append(w);f=model.forward(m0+kernel.U@m,independent=True)['y'].reshape(-1);pred.append(dict(material=m,old_full_residual=float(np.max(abs(f[:3]-yf.reshape(-1)[:3]))),reduced_full_error=float(np.max(abs(kernel.y0.reshape(-1)+w@B-f))),full_source_response=f))
 source=None
 if len(legal)>=2:
  d=(ws[0]-ws[1])@B;source=constrained_match(d,Q);source['predicted_gap']=float(abs(d@source['a']));source['independent_full_gap']=float(abs((pred[0]['full_source_response']-pred[1]['full_source_response'])@source['a']));source['old_span_max_gap']=float(np.max(abs(d[:3])))
 # Mandatory direct-completion reference for this probe: pivot columns of B until full rank.
 _,_,piv=qr(B,pivoting=True);chosen=list(range(3))
 for i in piv:
  if int(i) not in chosen and np.linalg.matrix_rank(B[:,chosen+[int(i)]])>np.linalg.matrix_rank(B[:,chosen]):chosen.append(int(i))
  if len(chosen)==6:break
 with book.op('DIRECT_HOST_COMPLETION',direct_linear_solve=1):
  w=np.linalg.solve(B[:,chosen].T,(yf-kernel.y0).reshape(-1)[chosen]);recovered=w/(R.reshape(-1)+w@H)
 out=dict(n=12,k=6,s=3,leading=lead,material_box_complete_reduced_rounded_model=bool(complete),legal_count=len(legal),classifications=classes,certificates=[{k:v for k,v in c.items() if k!='X'} for c in certs],full_reduced_error=float(np.max(abs(yf-yr))),nonlinearity_relative_to_first_Born=float(np.linalg.norm(full['j']-(m0+kernel.U@z)[:,None]*model.E0)/np.linalg.norm(full['j'])),predictions=pred,source=source,direct=dict(indices=chosen,new_source_count=len(chosen)-3,material_error=float(np.max(abs(recovered-z))),condition=float(np.linalg.cond(B[:,chosen]))),scope='Certificate applies to rounded reduced H/B/R/t; full n=12 physics independently validated in floating point. Host elimination rounding and Hankel/PDE error are not interval-certified in this probe.',CPU_seconds=time.process_time()-begin,receipt=book.receipt());save(OUT/'results/H12_probe.json',out);print({k:v for k,v in out.items() if k not in ('classifications','certificates','predictions','source','receipt')})
if __name__=='__main__':main()
