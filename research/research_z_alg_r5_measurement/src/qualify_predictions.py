"""Post-discovery audit with NO GT. No new roots, sources or measurements."""
from common import *
import certification as c,mpmath as mp,itertools

def decode(x):
 if isinstance(x,dict) and set(x)=={'real','imag'}:return complex(x['real'],x['imag'])
 if isinstance(x,dict):return {k:decode(v) for k,v in x.items()}
 if isinstance(x,list):return [decode(v) for v in x]
 return x
book=Ledger('interval_source_prediction_and_menu_audit');output={}
for case in ['T0','T1','T2','T3','T4','T5']:
 g=geometry(case);pos,Q=dictionary(case,g);p=Physics(g,pos,book);output[case]={}
 for arm in ('clean','noisy'):
  path=ROOT/'results'/case/arm;q=decode(json.loads((path/'Q.json').read_text()));design=decode(json.loads((path/'design_freeze.json').read_text()));preds=[];boxes=[]
  for idx,r in enumerate(q['certificates']):
   cr=r.get('certificate',{})
   if cr.get('status')!='CERTIFIED_LOCAL_REAL_ROOT':continue
   with c._interval_precision(128):
    safe=[{'lo':float(np.nextafter(v['lo'],np.inf)),'hi':float(np.nextafter(v['hi'],-np.inf))} for v in cr['intervals']['X']]
    if not all(v['lo']<k['lo'] and k['hi']<v['hi'] for v,k in zip(safe,cr['intervals']['K'])):raise AssertionError('Serialized inner-box inclusion not verified')
    X=[mp.iv.mpf([v['lo'],v['hi']]) for v in safe];boxes.append((idx,X))
    with book.op('INTERVAL_PREDICTION',interval_matrix_inversions=1,interval_complex_RHS=p.L):
     F,_,_,_=c._interval_forward_jacobian(X,p.H.tolist(),p.B.tolist(),p.R.tolist(),[0j]*p.L)
     yy=[mp.iv.mpc(F[2*j],F[2*j+1]) for j in range(p.L)]
     row=[]
     for src in design['menu']:
      y=sum((v*c._iv_complex(a) for v,a in zip(yy,src['a'])),mp.iv.mpc(0));row.append({'real':c._report_interval(y.real),'imag':c._report_interval(y.imag)})
     preds.append(dict(root_index=idx,certified_safe_binary64_box=safe,serialization_inclusion_check=True,response_boxes=row))
  pairs=[]
  for i,j in itertools.combinations(range(len(preds)),2):
   idx,X=boxes[i];jdx,Y=boxes[j];disjoint=any(x.b<y.a or y.b<x.a for x,y in zip(X,Y));gaps=[]
   with c._interval_precision(128):
    for a,b in zip(preds[i]['response_boxes'],preds[j]['response_boxes']):
     # Each outward report encloses exact interval. Distance between rectangles lower bound.
     da=[]
     for k in ('real','imag'):
      ar=mp.iv.mpf([a[k]['lo'],a[k]['hi']]);br=mp.iv.mpf([b[k]['lo'],b[k]['hi']]);da.append(c._abs_lower(ar-br))
     gap=mp.iv.sqrt(mp.iv.mpf(da[0])**2+mp.iv.mpf(da[1])**2)
     gaps.append(float(gap.a))
   chosen=design['choices']['maxmin'];pairs.append(dict(root_i=idx,root_j=jdx,material_boxes_disjoint=disjoint,source_gap_lower_bounds=gaps,chosen_gap_lower=gaps[chosen],chosen_noise_model_margin=gaps[chosen]-2*(EPS+MODEL_ALLOWANCE)))
  # Reconstruct every omitted pair direction for an audit; no new allowed source added.
  omitted=[];roots=[np.array(m) for m in q['accepted_roots']]
  for i,j in itertools.combinations(range(len(roots)),2):
   d=(p.dual(roots[i])-p.dual(roots[j]))@p.B;a=np.linalg.pinv(Q,rcond=1e-10)@d.conj();power=float(np.real(np.vdot(a,Q@a)));feed=float(np.linalg.norm(a)/np.sqrt(power)) if power>0 else None
   omitted.append(dict(pair=[i,j],unconstrained_power=power,normalized_feed_norm=feed,admissible=bool(power>1e-12 and feed<=4),reason='FEED_NORM_LIMIT' if power>1e-12 and feed>4 else 'LOW_OR_INVALID_POWER' if power<=1e-12 else 'ADMISSIBLE'))
  label='CERTIFIED_AMBIGUOUS' if any(t['material_boxes_disjoint'] for t in pairs) else 'CANDIDATE_ONLY' if roots else 'NO_FEASIBLE_CANDIDATE'
  output[case][arm]=dict(status=label,status_scope='At least two distinct real roots of rounded discrete old-data central equation, if certified ambiguous; no exhaustiveness.',response_boxes=preds,pairs=pairs,omitted_or_attempted_pair_directions=omitted,global_material_box_coverage='UNKNOWN',noise_feasible_set_coverage='UNKNOWN',interval_prediction_scope='Exact rounded G/E/S, exact binary64 source weights, certified local material boxes; no uncertain-data tube or physical-kernel error certificate.')
save(ROOT/'results/qualification.json',output);save(ROOT/'results/qualification_receipt.json',book.receipt());print({k:{a:v['status'] for a,v in d.items()} for k,d in output.items()})
