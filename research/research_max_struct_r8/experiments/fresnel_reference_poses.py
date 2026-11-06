from fresnel_validate import *
def main():
 start=time.process_time();rows=[]
 for name in ['dielTM_dec8f','twodielTM_8f']:
  d=prepare(ROOT/'data_local/fresnel'/f'{name}.exp');choices=[]
  for deg in ([0,90,180,270] if name.startswith('diel') else [0,90]):
   axis=np.array([np.cos(deg*np.pi/180),np.sin(deg*np.pi/180)]);cs=[-.03*axis] if name.startswith('diel') else [-.045*axis,.045*axis]
   pred,rec=reference(d,cs);train=d['train_mask'];gain=np.vdot(pred[train],d['Y'][train])/np.vdot(pred[train],pred[train])
   residual=float(np.linalg.norm((pred-d['Y'])[train])/np.linalg.norm(d['Y'][train]));gained=float(np.linalg.norm((gain*pred-d['Y'])[train])/np.linalg.norm(d['Y'][train]))
   choices.append(dict(degrees=deg,train_residual=residual,diagnostic_gain_train_residual=gained,gain=[float(gain.real),float(gain.imag)],prediction=pred))
  best=min(choices,key=lambda r:r['diagnostic_gain_train_residual']);mask=d['withheld_mask'];gain=complex(*best['gain']);best['selected_pose_heldout_residual']=float(np.linalg.norm((gain*best['prediction']-d['Y'])[mask])/np.linalg.norm(d['Y'][mask]))
  for row in choices:row.pop('prediction')
  rows.append(dict(case=name,selection='TRAINING ONLY finite ambiguity of unlabeled primary schematic axes; gain/pose diagnostic never used by inverse arms',choices=choices,selected_degrees=best['degrees']));print(json.dumps(rows[-1]),flush=True)
 (ROOT/'RESULTS/fresnel/REFERENCE_POSE_DIAGNOSTIC.json').write_text(json.dumps(dict(cases=rows,cpu_seconds=time.process_time()-start),indent=2))
if __name__=='__main__':main()
