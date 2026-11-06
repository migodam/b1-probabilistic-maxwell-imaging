from fresnel_validate import *
from scipy.linalg import svd

def main():
 cpu=time.process_time();rows=[];data={};refs={}
 for name,cs in [('dielTM_dec8f',[np.array([0.,.03])]),('twodielTM_8f',[np.array([0.,-.045]),np.array([0.,.045])]),('FoamTwinDielTM',None)]:
  foam=name=='FoamTwinDielTM';frequency=8 if foam else 2;d=prepare(ROOT/'data_local/fresnel'/f'{name}.exp',frequency=frequency,side=60 if foam else 24,width=.2 if foam else .16,matrix_free=foam);data[name]=d
  if not foam:refs[name]=reference(d,cs)[0]
  nr,ns=d['Y'].shape;k=2*np.pi*frequency*1e9/299792458.;ra=2*np.pi*np.arange(nr)/nr;rx=(1.67 if foam else .76)*np.c_[np.cos(ra),np.sin(ra)];sa=2*np.pi*np.arange(ns)/ns;src=(1.67 if foam else .72)*np.c_[np.cos(sa),np.sin(sa)];cal=[]
  for s in range(ns):
   axis=-src[s]/np.linalg.norm(src[s]);normal=np.array([-axis[1],axis[0]])
   def B(points):
    v=points-src[s];dist=np.linalg.norm(v,axis=1);ang=np.arctan2(v@normal,v@axis);return np.stack([.25j*hankel1(abs(m),k*dist)*np.exp(1j*m*ang) for m in range(-2,3)],axis=1)
   fit=d['calibration'][s]['fit_channels'];A=B(rx);scale=np.linalg.norm(A[fit],axis=0);Z=A/scale;u,q,vh=svd(Z[fit],full_matrices=False);keep=q>1e-4*q[0];inverse=(vh[keep].conj().T/q[keep])@u[:,keep].conj().T;coef=inverse@d['I'][fit,s]
   C=(B(d['centers'])/scale)@inverse;E=C@d['I'][fit,s]
   cal.append(dict(source=s+1,retained_rank=int(keep.sum()),retained_min_relative_sigma=float(q[keep][-1]/q[0]),
    normalized_coefficient_norm=float(np.linalg.norm(coef)),normalized_fit_data_norm=float(np.linalg.norm(d['I'][fit,s])),
    relative_coefficient_amplification=float(np.linalg.norm(coef)/np.linalg.norm(d['I'][fit,s])),
    interior_relative_calibration_condition=float(np.linalg.norm(C,2)*np.linalg.norm(d['I'][fit,s])/np.linalg.norm(E))))
  rows.append(dict(case=name,source_diagnostics=cal))
 # Additional-known-reference gain calibration is a diagnostic separate from
 # all inverse arms. Fit single known nominal material using train channels;
 # transfer exactly same single complex coefficient to twin, without twin fit.
 single=data['dielTM_dec8f'];train=single['train_mask'];p=refs['dielTM_dec8f'];gain=np.vdot(p[train],single['Y'][train])/np.vdot(p[train],p[train]);transfer=[]
 for name in refs:
  d=data[name]
  p=gain*refs[name];mask=d['withheld_mask'];transfer.append(dict(case=name,heldout_relative=float(np.linalg.norm((p-d['Y'])[mask])/np.linalg.norm(d['Y'][mask]))))
 record=dict(illumination_stability=rows,known_reference_transfer=dict(gain_real=float(gain.real),gain_imag=float(gain.imag),fit_dataset='single known nominal eps3 radius15mm +y30mm TRAIN ONLY',validation=transfer,used_in_main_inversion=False,qualification='Additional known-material and nominal-pose prior; material/pose uncertainty can enter gain; not unique antenna calibration'),cpu_seconds=time.process_time()-cpu)
 (ROOT/'RESULTS/fresnel/ILLUMINATION_STABILITY_AND_REFERENCE_TRANSFER.json').write_text(json.dumps(record,indent=2))
 print(json.dumps(dict(relative_condition_median=[float(np.median([r['interior_relative_calibration_condition'] for r in x['source_diagnostics']])) for x in rows],known_reference_transfer=record['known_reference_transfer'])),flush=True)
if __name__=='__main__':main()
