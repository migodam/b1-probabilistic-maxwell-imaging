from complete_quotient import *
def readz(z):return complex(z['real'],z['imag'])
def rangez(z):return mp.iv.mpc(mp.iv.mpf([z['real']['lo'],z['real']['hi']]),mp.iv.mpf([z['imag']['lo'],z['imag']['hi']]))
def gap(a,b):
 z=a-b
 def dist(v):return v.a if v.a>0 else (-v.b if v.b<0 else mp.iv.mpf(0))
 q=mp.iv.sqrt(dist(z.real)**2+dist(z.imag)**2);return float(np.nextafter(float(q.a),-np.inf))
def main():
 start=time.process_time();mp.iv.prec=128;book=Ledger('complete_fibre_source_gaps');rows=[]
 for case in [f'T{i}' for i in range(6)]:
  g=geometry(case);pos,Q=dictionary(case,g);p=Physics(g,pos,book)
  for arm in ['clean','noisy']:
   c=json.loads((OUT/'results'/f'complete_{case}_{arm}_128.json').read_text());valid=[cl for cl in c['classifications'] if cl['status']=='CERTIFIED_LEGAL_REAL_MATERIAL']
   if len(valid)<2:continue
   de=json.loads((BASE/'results'/case/arm/'design_freeze.json').read_text());a=np.array([readz(z) for z in de['menu'][de['choices']['maxmin']]['a']]);policies={'original_frozen_maxmin':a}
   if case=='T0' and arm=='noisy':
    designs=json.loads((OUT/'results/constrained_source_design_freeze.json').read_text())['designs'];policies.update({k:np.array([readz(z) for z in v['a']]) for k,v in designs.items()})
   for name,a in policies.items():
    # Exact binary64 source dictionary and weights, interval multiplication before response.
    b=[sum((ivc(p.B[i,j])*ivc(a[j]) for j in range(p.L)),mp.iv.mpc(0)) for i in range(p.n)]
    predictions=[sum((rangez(z)*v for z,v in zip(cl['w'],b)),mp.iv.mpc(0)) for cl in valid];gaps=[gap(u,v) for u,v in itertools.combinations(predictions,2)];rows.append(dict(case=case,arm=arm,policy=name,legal_roots=len(valid),pairwise_gap_lower=gaps,min_gap_lower=min(gaps),covers_all_regular_legal_materials_for_exact_old_data=c['material_box_complete'],new_measurement_error_half_gap=min(gaps)/2,certified_predictions=[interval_json(z) for z in predictions],scope='Exact old data fixed at this center; new measurement absolute complex error < half min gap; nominal rounded model; does NOT cover uncertainty in old data or physical model.'))
 save(OUT/'results/complete_fibre_source_gaps.json',dict(rows=rows,CPU_seconds=time.process_time()-start,receipt=book.receipt()));print([{k:v for k,v in r.items() if k not in ('certified_predictions','scope')} for r in rows])
if __name__=='__main__':main()
