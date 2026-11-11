"""Plot only already sealed/scored material arrays; zero FW or optimization."""
from pathlib import Path
from datetime import datetime,timezone
import json,time

HERE=Path(__file__).resolve().parent


def read(path):return json.loads(Path(path).read_text())


def main():
    cpu,wall=time.process_time(),time.perf_counter()
    prior_receipt=read(HERE/'RAW_VS_CORRECTOR_RECEIPT.json') if (HERE/'RAW_VS_CORRECTOR_RECEIPT.json').exists() else None
    # The plane/order/scale are fixed for every case before loading its labels.
    config=dict(status='FIXED_DISPLAY_OF_SAVED_SCORED_MATERIALS',
        created_utc=datetime.now(timezone.utc).isoformat(),cases=['random_high8','random_mixed2_8','random_high8_crossgrid8'],
        columns=['GT volume mean','DR raw original-Phi selection','DR corrected','CIE corrected'],
        slice_axis='z',slice_index=3,grid_side=6,slice_coordinate=.05,width=.6,
        color_limits=[0,8],selection='raw is saved minoriginalPhi proposal, neverGTbest',
        new_forward=0,new_optimization=0)
    (HERE/'RAW_VS_CORRECTOR_CONFIG.json').write_text(json.dumps(config,indent=2)+'\n')
    seal=read(HERE/'ALL12_PREDICTIONS_SEALED.json')
    assert seal['status']=='ALL12_PREDICTIONS_AND_PREFIX_CHOICES_SEALED'
    scored=read(HERE/'OFFLINE_SCORES.json')
    assert scored['status']=='OFFLINE_SCORED_AFTER_ALL12_PUBLIC_SEAL'
    decomposition=read(HERE/'RAW_CORRECTOR_DECOMPOSITION.json')
    table={(r['case_id'],r['method']):r for r in decomposition['rows']}
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(3,4,figsize=(15.5,11.8),layout='constrained')
    provenance=[]
    names=['high χ8 matched','mixed χ2/8 matched','high χ8, 8³ → 6³ pair']
    for row,(cid,label) in enumerate(zip(config['cases'],names)):
        base=HERE/'RESULTS'/cid
        with np.load(HERE/'sealed_truth'/f'{cid}.npz',allow_pickle=False) as z:truth=z['chi'].copy()
        with np.load(base/'DR/prediction.npz',allow_pickle=False) as z:raw=z['raw'].copy();dr=z['prediction'].copy()
        with np.load(base/'CIE_sourcewise/prediction.npz',allow_pickle=False) as z:cie=z['prediction'].copy()
        drinfo=table[(cid,'DR')];cieinfo=table[(cid,'CIE_sourcewise')]
        with np.load(base/'DR/paid_materials.npz',allow_pickle=False) as z:
            assert np.array_equal(z['m'][drinfo['raw']['selected_history_index']],raw)
        history=read(base/'DR/history.json')
        proposal=[r for r in history if r['phase']!='corrector']
        assert min(r['phi'] for r in proposal)==drinfo['raw']['Phi']
        panels=[(truth,None,'GT volume mean'),
                (raw,drinfo['raw'],f"DR raw: min-Φ proposal\nprojection iteration {drinfo['raw']['selected_projection_iteration']}"),
                (dr,drinfo['corrected'],'DR + same-Φ corrector'),
                (cie,cieinfo['corrected'],'CIE + same-Φ corrector')]
        for col,(material,info,title) in enumerate(panels):
            ax=axes[row,col]
            im=ax.imshow(np.asarray(material).real.reshape(6,6,6)[:,:,3].T,
                         origin='lower',cmap='viridis',vmin=0,vmax=8)
            ax.set_xticks([]);ax.set_yticks([]);ax.set_title(title,fontsize=11)
            if col==0:
                ax.set_ylabel(label,fontsize=11)
                ax.set_xlabel('Reference, same χ color scale\nfixed upper-central xy plane z=+0.05',fontsize=9)
            else:
                noise=info['noise'] if col==1 else info
                passed=noise['raw_CN99_pass']
                ratio=noise['residual_over_CN99_delta']
                status='IN' if passed else 'OUT'
                ax.set_xlabel(f"relL2={info['relative_material_l2']:.4f}   "
                              f"D={info['data_term']:.3g}   R={info['prior']:.3g}\n"
                              f"physical CN99: {status}  (||r||/δ={ratio:.3f})",
                              fontsize=9,color='#16733b' if passed else '#a32121')
        provenance.append(dict(case_id=cid,raw_prediction=str(base/'DR/prediction.npz'),
            raw_history_index=drinfo['raw']['selected_history_index'],
            raw_projection_iteration=drinfo['raw']['selected_projection_iteration'],
            raw_Phi=drinfo['raw']['Phi'],raw_fields_or_material_optimized_again=False))
    fig.colorbar(im,ax=list(axes.ravel()),shrink=.7,pad=.012,label='χ, common [0,8] scale')
    fig.suptitle('Saved high-contrast true3D materials: raw DR versus the same-objective corrector',fontsize=14,y=.995)
    fig.get_layout_engine().set(rect=(0,.065,1,.88))
    fig.text(.5,.014,'D = inherited source-normalized data term; R = unchanged quadratic prior.  '
             'Raw is selected by original Φ, never by GT.\n'
             'CN99 is the physical raw-noise ball; the cross-grid pair has additional model mismatch.  '
             'Saved arrays only: 0 new Maxwell solves / 0 optimization.',
             ha='center',va='bottom',fontsize=9)
    fig.savefig(HERE/'RAW_VS_CORRECTOR.png',dpi=160);plt.close(fig)
    receipt=dict(status='COMPLETE_SAVED_ARRAY_FIGURE',created_utc=datetime.now(timezone.utc).isoformat(),
        cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
        clock_scope='postprocessing start beforeNumPy/Matplotlib imports throughfigureoutput',
        new_forward=0,new_Maxwell_LU=0,new_source_RHS=0,new_adjoint_RHS=0,new_optimization=0,
        truth_access='offline previouslyscored material reference only; noselection',
        cases=provenance,config='RAW_VS_CORRECTOR_CONFIG.json',output='RAW_VS_CORRECTOR.png')
    (HERE/'RAW_VS_CORRECTOR_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    costfile=HERE/'RAW_VS_CORRECTOR_REVISION_COSTS.json'
    costlog=read(costfile) if costfile.exists() else dict(runs=[])
    if prior_receipt is not None and not costlog['runs']:
        costlog['runs'].append(dict(prior_receipt,role='initial figure; title overlap found in visual review'))
    costlog['runs'].append(dict(receipt,role='current figure; reserved title and footer layout'))
    costfile.write_text(json.dumps(costlog,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ('status','cpu_seconds','wall_seconds','new_forward')}),flush=True)


if __name__=='__main__':main()
