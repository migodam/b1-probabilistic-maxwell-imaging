"""Offline-only scientific scores. Truth opens only after the online seal."""
from pathlib import Path
import argparse,json,sys,time,csv
import numpy as np
from generate import ROOT,plain,write_json,maps

def score_endpoint(endpoint,truth,y):
    a=np.load(endpoint);z=np.asarray(a['z']);zt=truth['z'];estimated=maps(z);actual=maps(zt)
    m=estimated['target_mask'];mt=actual['target_mask'];union=m|mt
    pred=np.asarray(a['data']);clean=truth['clean'];sigma=.01
    c=estimated['target_material']-1;ct=actual['target_material']-1
    # Sorted semiaxes remove the exact ellipse axis/angle parametrization symmetry.
    axes=np.sort(z[14:16]);taxes=np.sort(zt[14:16]);angle_error=None
    if abs(taxes[1]-taxes[0])>1e-5:
        aa=z[18]+(np.pi/2 if z[14]<z[15] else 0);at=zt[18]+(np.pi/2 if zt[14]<zt[15] else 0)
        angle_error=float(abs((aa-at+np.pi/2)%np.pi-np.pi/2))
    white=(pred-y)/sigma;cleanwhite=(pred-clean)/sigma
    return dict(endpoint=str(endpoint.relative_to(ROOT)) if endpoint.is_relative_to(ROOT) else str(endpoint),
        position_error=float(np.linalg.norm(z[12:14]-zt[12:14])),x_error=float(z[12]-zt[12]),y_error=float(z[13]-zt[13]),
        semiaxis_absolute_errors=np.abs(axes-taxes).tolist(),semiaxis_relative_errors=(np.abs(axes-taxes)/taxes).tolist(),semiaxis_RMSE=float(np.sqrt(np.mean((axes-taxes)**2))),
        epsilon_real_absolute_error=float(abs(z[16]-zt[16])),epsilon_real_relative_error=float(abs(z[16]-zt[16])/(1+zt[16])),contrast_relative_error=float(abs(z[16]-zt[16])/zt[16]),
        loss_absolute_error=float(abs(z[17]-zt[17])),epsilon_complex_relative_error_f1=float(abs((z[16]-zt[16])+1j*(z[17]-zt[17]))/abs(1+zt[16]+1j*zt[17])),target_complex_intrinsic_relative_error=float(abs((z[16]-zt[16])+1j*(z[17]-zt[17]))/abs(zt[16]+1j*zt[17])),
        shape_IoU=float(np.sum(m&mt)/max(np.sum(union),1)),shape_mask_relative_L2=float(np.sqrt(np.sum(m!=mt)/max(np.sum(mt),1))),target_material_contrast_relative_L2=float(np.linalg.norm(c-ct)/max(np.linalg.norm(ct),1e-300)),
        target_material_RMSE_union=float(np.sqrt(np.mean(abs(c[union]-ct[union])**2))) if union.any() else 0.,
        wall_epsilon_RMSE=float(np.sqrt(np.mean((z[:6]-zt[:6])**2))),wall_loss_RMSE=float(np.sqrt(np.mean((z[6:12]-zt[6:12])**2))),
        observation_white_norm=float(np.linalg.norm(white)),observation_white_RMS=float(np.linalg.norm(white)/np.sqrt(2*white.size)),
        clean_prediction_white_norm=float(np.linalg.norm(cleanwhite)),clean_prediction_white_RMS=float(np.linalg.norm(cleanwhite)/np.sqrt(2*cleanwhite.size)),
        angle_error_radians=angle_error,estimated_parameters=z.tolist(),truth_parameters=zt.tolist(),
        gates=dict(position_le_0p1=bool(np.linalg.norm(z[12:14]-zt[12:14])<=.1),shape_IoU_ge_0p7=bool(np.sum(m&mt)/max(np.sum(union),1)>=.7),epsilon_relative_le_0p2=bool(abs(z[16]-zt[16])/(1+zt[16])<=.2)))

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--case',required=True);p.add_argument('--repeat',type=int,default=0);p.add_argument('--allow-dev-exposure',action='store_true');a=p.parse_args();out=a.output.resolve()
    if not (out/'ONLINE_SEALED.json').exists():raise RuntimeError('Refusing private truth read before ONLINE_SEALED.json')
    if not (out/'endpoint.npz').exists():raise FileNotFoundError(out/'endpoint.npz')
    seal=json.loads((out/'ONLINE_SEALED.json').read_text())
    if a.case.startswith('blind') and not any(p.exists() for p in [ROOT/'METHODS_FROZEN.json',ROOT/'configs/METHODS_FROZEN.json']):
        raise RuntimeError('Blind truth remains closed until definitive METHODS_FROZEN.json exists')
    truth=np.load(ROOT/'inputs/private'/a.case/'truth.npz');y=np.load(ROOT/'inputs/online'/f'{a.case}_r{a.repeat}'/'observation.npz')['y']
    selected=score_endpoint(out/'endpoint.npz',truth,y)
    summary=seal.get('summary',{})
    selected_rec=next((r for r in summary.get('starts',[]) if r.get('stage')==summary.get('selected_start')),{})
    selected.update(final_exact_projected_gradient=selected_rec.get('final_exact_projected_gradient'),stop=selected_rec.get('stop'))
    selected['joint_imaging_gate']='PASS' if all(selected['gates'].values()) else 'FAIL'
    selected.update(method=summary.get('method',out.name),objective=summary.get('objective'),converged=summary.get('converged'),selected_start=summary.get('selected_start'))
    branches=[]
    for endpoint in sorted(set(out.glob('branches/**/endpoint.npz')) | set(out.glob('START_*_ENDPOINT.npz'))):
        branch=score_endpoint(endpoint,truth,y)
        if endpoint.name.startswith('START_'):
            k=int(endpoint.name.split('_')[1]);rec=next((v for v in summary.get('starts',[]) if v.get('stage')==f'joint_start_{k}'),{})
            branch.update(objective=rec.get('objective'),stop=rec.get('stop'),evaluations=rec.get('evaluations'),accepted_steps=rec.get('accepted_steps'),final_exact_projected_gradient=rec.get('final_exact_projected_gradient'))
        branches.append(branch)
    result=dict(case=a.case,repeat=a.repeat,method=selected['method'],cohort=('ORACLE_DIAGNOSTIC' if seal.get('diagnostic_only') else 'DEVELOPMENT_EXPOSED' if a.case=='dev' else 'INDEPENDENT_FROZEN_BLIND'),status='POST_SEAL_SCORED',scored_unix=time.time(),online_seal=seal,selected=selected,branches=branches,
        selection_rule='online endpoint fixed before scoring; private errors never select a branch',scope='independent off-grid draw within declared six-cell wall and ellipse family; no arbitrary support or calibrated posterior claim')
    write_json(out/'SCORES.json',result)
    em=maps(np.load(out/'endpoint.npz')['z']);np.save(out/'evaluation_target_material.npy',em['target_material']);np.save(out/'estimated_wall.npy',em['wall_eps']+1j*em['wall_loss']);np.save(out/'target_mask.npy',em['target_mask'])
    np.savez_compressed(out/'material_maps.npz',**em)
    write_json(out/'EVALUATION_MAPS_FORMAT.json',dict(runtime_target_material='target_material.npy is preserved in runtime channel format declared by reconstruction.json',evaluation_target_material='evaluation_target_material.npy: complex epsilon at f1 on full evaluation grid',estimated_wall='estimated_wall.npy: complex wall epsilon at f1 on same evaluation grid',axes='pixel_x and pixel_y in material_maps.npz',scoring_grid=[401,301],inference_unchanged=True))
    print(json.dumps(dict(case=a.case,repeat=a.repeat,selected={k:selected[k] for k in ['position_error','shape_IoU','epsilon_real_relative_error','loss_absolute_error','wall_epsilon_RMSE','observation_white_RMS']},branches=len(branches))))
if __name__=='__main__':main()
