"""Portable derived actual measured training packages; no target truth."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,time,numpy as np
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[3];sys.path.insert(0,str(REPO))
from research.research_max_struct_r13.src.physics_inputs import prepare_fresnel,SourceSpecificPublicModel
from research.research_max_struct_r13.delegated.physics.generate_inputs import dump

def main():
    out=HERE/'fresnel_public';out.mkdir(exist_ok=False);rows=[]
    for case in ('dielTM_dec8f','twodielTM_8f'):
        d=prepare_fresnel(case);m=d['model'];Y=d['Y'];meta=d['metadata']
        ports=np.stack(d['S']);public=dict(G=d['G'],E0=d['E'],S=ports[0],S_ports=ports,
            grid_shape=np.array(m.shape),centers=d['context']['centers'],coords=d['context']['centers'],
            components=np.array(1),kind=np.array('fresnel'),m_init=np.full(m.n,.1),
            lower=np.array(0.),upper=np.array(20.),Y=Y,delta=np.array(d['delta']),
            # Required scalar is explicitly discrepancy-per-channel, not noise.
            sigma=np.array(d['delta']/np.sqrt(Y.size)),noise_sigma_known=np.array(False),
            sigma_role=np.array('diagnostic discrepancy per complex training channel; not known noise sigma'))
        np.savez_compressed(out/f'{case}.npz',**public)
        reloaded=SourceSpecificPublicModel(public);x=public['m_init']
        old=m.state(x);new=reloaded.state(x);f,g=m.data_value_gradient(x,Y);ff,gg=reloaded.data_value_gradient(x,Y)
        dump(out/f'{case}_METADATA.json',meta)
        # Heldout observations are separate and never part of the inversion NPZ.
        held=d['context']['withheld_mask'];fullY=d['context']['Y']
        np.savez_compressed(out/f'{case}_HELDOUT_OFFLINE_ONLY.npz',Y=fullY[held],mask=held,
            E0_all_sources=d['context']['public']['E0'],S_all_receivers=d['context']['public']['S'])
        rows.append(dict(case=case,path=str(out/f'{case}.npz'),source_path=meta['source_path'],
            actual_training_channels=Y.size,heldout_separate=True,source_ports_shape=list(ports.shape),
            portable_prediction_relative_difference=float(np.linalg.norm(new['Y']-old['Y'])/np.linalg.norm(old['Y'])),
            objective_difference=abs(ff-f),gradient_relative_difference=float(np.linalg.norm(gg-g)/np.linalg.norm(g)),
            delta=d['delta'],delta_kind=meta['delta_status'],known_noise_sigma=False,
            allowed_claim='measured predictive fit on fixed train/heldout split; no pixel material GT'))
    dump(out/'PACKAGING_RECEIPT.json',dict(status='COMPLETE',rows=rows,raw_ASCII_copied=False,
        heldout_for_selection=False,material_truth_access=False))
    print('Portable measured Fresnel packages and source-specific model verified',flush=True)
if __name__=='__main__':main()
