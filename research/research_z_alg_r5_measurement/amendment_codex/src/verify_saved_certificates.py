"""Recompute serialized complex-root boxes using real 2N Jacobians.
Independent representation check, not independent external certification software.
"""
from complete_quotient import *
def real_block(J):return [[z.real for z in row]+[-z.imag for z in row] for row in J]+[[z.imag for z in row]+[z.real for z in row] for row in J]
def main():
 start=time.process_time();mp.iv.prec=128;rows=[];book=Ledger('saved_certificate_real_form_audit')
 for case in [f'T{i}' for i in range(6)]:
  g=geometry(case);pos,Q=dictionary(case,g);p=Physics(g,pos,book)
  for arm in ['clean','noisy']:
   doc=json.loads((OUT/'results'/f'complete_{case}_{arm}_128.json').read_text());y=np.load(BASE/'inputs'/f'{case}.npz')['old_'+arm];C=chart(p.H,p.B[:,:p.s],p.R,y)['coeff'];oks=[];bounds=[]
   for cr in doc['certificates']:
    X=[mp.iv.mpc(mp.iv.mpf([v['real']['lo'],v['real']['hi']]),mp.iv.mpf([v['imag']['lo'],v['imag']['hi']])) for v in cr['domain']];xr=[z.real for z in X]+[z.imag for z in X];xc=[ivc(mid(z)) for z in X];point=[z.real for z in xc]+[z.imag for z in xc];F,J=evaluate(C,xc);fr=[z.real for z in F]+[z.imag for z in F];Jr=real_block(J);pre=np.linalg.inv(np.array([[float(v.mid) for v in row] for row in Jr]));_,JX=evaluate(C,X);A=real_block(JX);N=len(xr);M=[[mp.iv.mpf(i==j)-sum((mp.iv.mpf(float(pre[i,k]))*A[k][j] for k in range(N)),mp.iv.mpf(0)) for j in range(N)] for i in range(N)];K=[point[i]-sum((mp.iv.mpf(float(pre[i,k]))*fr[k] for k in range(N)),mp.iv.mpf(0))+sum((M[i][j]*(xr[j]-point[j]) for j in range(N)),mp.iv.mpf(0)) for i in range(N)];norm=max(sum((upper(v) for v in row),mp.iv.mpf(0)) for row in M);inside=all(k.a>x.a and k.b<x.b for k,x in zip(K,xr));oks.append(bool(inside and norm<1));bounds.append(float(norm.b))
   rows.append(dict(case=case,arm=arm,all_serialized_boxes_recertified=all(oks),count=len(oks),max_real_contraction=max(bounds)))
 save(OUT/'results/saved_certificate_audit.json',dict(passed=all(r['all_serialized_boxes_recertified'] for r in rows),rows=rows,CPU_seconds=time.process_time()-start,receipt=book.receipt()));print(rows);assert all(r['all_serialized_boxes_recertified'] for r in rows)
if __name__=='__main__':main()
