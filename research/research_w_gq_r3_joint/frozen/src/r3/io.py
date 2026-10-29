"""Online artifacts contain measurements/candidates only; no truth reader."""
from pathlib import Path
import json,time
import numpy as np
from .accounting import plain
class Trace:
    def __init__(self,out):self.out=Path(out);self.out.mkdir(parents=True,exist_ok=True);self.count=0
    def record(self,z,pred,objective,**meta):
        self.count+=1;name=f'candidate_{self.count:05d}.npz';np.savez_compressed(self.out/name,z=z,data=pred.data,**({'jacobian':pred.jacobian} if pred.jacobian is not None else {}))
        row=dict(time_unix=time.time(),candidate=self.count,file=name,objective=float(objective),**meta)
        with (self.out/'TRACE.jsonl').open('a') as f:f.write(json.dumps(row,default=plain,allow_nan=False)+'\n')
        return self.count
    def event(self,**meta):
        with (self.out/'STEPS.jsonl').open('a') as f:f.write(json.dumps(dict(time_unix=time.time(),**meta),default=plain,allow_nan=False)+'\n')
def save_endpoint(model,z,pred,out,summary):
    out=Path(out);maps=model.material_maps(z);np.savez_compressed(out/'endpoint.npz',z=z,data=pred.data,**maps)
    np.savez_compressed(out/'predicted_measurement.npz',data=pred.data,y_pred=pred.data,frequencies=np.asarray(model.config['frequencies']))
    if pred.jacobian is not None:
        from .inverse import prior
        rp,jp=prior(model,z);np.savez_compressed(out/'endpoint_exact_jacobian.npz',z=z,data=pred.data,exact_jacobian=pred.jacobian,parameter_scales=np.asarray(model.spec['scale']),residual=pred.data-model.observation,measurement=model.observation,prior_residual=rp,prior_jacobian=jp,sigma=model.config['sigma'])
    ix=(maps['pixel_x']>=-1.2)&(maps['pixel_x']<=1.2);iy=(maps['pixel_y']>=1.4)&(maps['pixel_y']<=2.8);material=np.stack([maps['image_eps'][np.ix_(iy,ix)],maps['image_loss'][np.ix_(iy,ix)]],axis=-1);np.save(out/'target_material.npy',material)
    reconstruction=dict(summary=summary,config=model.config,parameter_names=model.spec['names'],z=z,target_material_format='array(y,x,channels): real epsilon, loss ell at reference f1',target_pixel_x=maps['pixel_x'][ix],target_pixel_y=maps['pixel_y'][iy]);(out/'reconstruction.json').write_text(json.dumps(reconstruction,default=plain,indent=2)+'\n')
    import csv
    decisions={}
    if (out/'STEPS.jsonl').exists():
        for line in (out/'STEPS.jsonl').read_text().splitlines():
            row=json.loads(line)
            if row.get('event')=='TR_DECISION':decisions[row['candidate']]=row
    fields=['candidate','time_unix','stage','kind','objective','accepted','actual_decrease','predicted_decrease','rho','projected_gradient','file']
    with (out/'optimization_trace.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        for line in (out/'TRACE.jsonl').read_text().splitlines():
            row=json.loads(line);row.update({k:v for k,v in decisions.get(row['candidate'],{}).items() if k in fields});writer.writerow({k:row.get(k,'') for k in fields})
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for name,box in [('reconstruction',(-1.5,1.5,.4,2.8)),('target',(-1.2,1.2,1.4,2.8))]:
        fig,axes=plt.subplots(1,2,figsize=(9,4),constrained_layout=True)
        for ax,key,title in zip(axes,['image_eps','image_loss'],['Reconstructed real permittivity','Reconstructed loss parameter']):
            im=ax.imshow(maps[key],origin='lower',extent=[maps['pixel_x'][0],maps['pixel_x'][-1],maps['pixel_y'][0],maps['pixel_y'][-1]],aspect='equal',cmap='magma');ax.set(xlim=box[:2],ylim=box[2:],xlabel='x / wavelength',ylabel='y / wavelength',title=title);fig.colorbar(im,ax=ax,shrink=.8)
        fig.savefig(out/(name+'.png'),dpi=180,bbox_inches='tight');plt.close(fig)
    (out/'SUMMARY.json').write_text(json.dumps(summary,default=plain,indent=2)+'\n')
    (out/'ONLINE_SEALED.json').write_text(json.dumps(dict(time_unix=time.time(),selected_by='minimum full measurement likelihood plus frozen proper prior',truth_access=False,parameter_names=model.spec['names'],summary=summary),default=plain,indent=2)+'\n')
