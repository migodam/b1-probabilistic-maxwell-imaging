#!/usr/bin/env python3
"""Render approved saved R15 selections/diagnostic JSON only; no rescoring/physics."""
from __future__ import annotations
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import argparse,json,sys,time,shutil
from pathlib import Path
from datetime import datetime,timezone
START_CPU=time.process_time();START_WALL=time.perf_counter()
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(Path(__file__).resolve().parent))
import analyze_campaign as campaign
base=campaign.base


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=ROOT);p.add_argument('--seal',type=Path)
    a=p.parse_args();root=a.root.resolve();seal=(a.seal or root/'predictions/PREDICTION_SEAL_MAIN_APPROVED_V1.json').resolve()
    campaign.authorize(root,seal)  # approved24x11 predictions before any GT load
    analysis=root/'analysis/twofold_exposed24_v1';main_receipt=base.read(analysis/'RECEIPT.json')
    if Path(main_receipt['seal']).resolve()!=seal:raise ValueError('Approved seal differs from scored selection source')
    out=root/'FIGURES/final';out.mkdir(parents=True,exist_ok=True)
    os.environ['MPLCONFIGDIR']=str(out/'.matplotlib')
    record=dict(started_utc=datetime.now(timezone.utc).isoformat(),status='RUNNING',approved_seal=str(seal),
        scoring=False,new_forward=0,new_LU=0,new_RHS=0,new_SVD=0,optimization=False,selection_modified=False,
        sources=[],artifacts=[],issues=[],scope='Presentation of already-scored public primary/raw/corrected selections and saved offline diagnostics; scientific acceptance remains parent')
    try:
        import numpy as np
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'pdf.fonttype':42})
        allviews_path=analysis/'ALL_VIEWS.json';allviews=base.read(allviews_path)
        record['sources'].append(str(allviews_path))
        lookup={(r['case_id'],r['solver'],r['view']):r for r in allviews if r['view'] in base.VIEWS}
        material={};truth={};vmax={};selected_sources=[];gt_sources=[]
        def load_selected(cid,method,view):
            key=(cid,method,view)
            if key in material:return material[key]
            row=lookup.get(key);array=None
            try:
                if row is None:raise ValueError('Selected row absent')
                path=Path(row['source_file'])
                if not path.resolve().is_relative_to(root/'RESULTS/twofold_exposed24_v1'):raise ValueError('Selected source outside approved campaign')
                with np.load(path,allow_pickle=False) as z:array=z['m'].copy().reshape(24,24)
                selected_sources.append(dict(case_id=cid,method=method,view=view,candidate=row.get('selection_candidate'),source_view=row.get('source_view'),path=str(path),bytes=path.stat().st_size))
            except Exception as exc:record['issues'].append(dict(case_id=cid,method=method,view=view,error=repr(exc)))
            material[key]=array;return array
        for cid in campaign.CASES:
            path=root.parent/'research_max_struct_r14/sealed_truth/independent24'/(cid+'.npz')
            # Material pixels only, authorized by the existing approved main seal.
            with np.load(path,allow_pickle=False) as z:truth[cid]=z['material_coarse'].copy().reshape(24,24)
            gt_sources.append(dict(case_id=cid,path=str(path),key='material_coarse'))
            arrays=[truth[cid]]+[load_selected(cid,m,'primary') for m in campaign.METHODS]
            vmax[cid]=max(float(np.max(x)) for x in arrays if x is not None)
        record['selected_material_sources']=selected_sources;record['authorized_GT_material_sources']=gt_sources
        def save(fig,name,caption):
            for ext in ('png','pdf'):
                target=out/(name+'.'+ext);fig.savefig(target,dpi=160,bbox_inches='tight',pad_inches=.10)
                record['artifacts'].append(dict(path=str(target.relative_to(root)),bytes=target.stat().st_size,caption=caption))
            plt.close(fig)
        def show(ax,array,maximum):
            ax.set_xticks([]);ax.set_yticks([])
            if array is None:ax.text(.5,.5,'MISSING',ha='center',va='center',transform=ax.transAxes);return None
            return ax.imshow(array,origin='lower',vmin=.02,vmax=maximum,cmap='viridis',interpolation='nearest')
        def fmt(value):return f'{value:.5g}' if base.finite(value) else 'NA'
        def metric_text(cid,method,view):
            row=lookup.get((cid,method,view),{})
            return 'NRMSE '+fmt(row.get('NRMSE'))+'   IoU '+fmt(row.get('IoU'))+'\nresidual/δ '+fmt(row.get('residual_over_delta'))
        # All24 cases/all11 selected primary methods; top=.95 removes old blank band.
        for page,start in enumerate((0,12),1):
            cases=campaign.CASES[start:start+12];fig=plt.figure(figsize=(24,25.2))
            grid=fig.add_gridspec(12,13,width_ratios=[1]*12+[.07],left=.045,right=.985,bottom=.025,top=.95,hspace=.13,wspace=.065)
            for i,cid in enumerate(cases):
                im=None
                for j,method in enumerate(('GT',)+campaign.METHODS):
                    ax=fig.add_subplot(grid[i,j]);image=show(ax,truth[cid] if method=='GT' else load_selected(cid,method,'primary'),vmax[cid])
                    if image is not None:im=image
                    if i==0:ax.set_title(method,fontsize=8,pad=8)
                    if j==0:ax.set_ylabel(cid,fontsize=9)
                cb=fig.add_subplot(grid[i,12]);fig.colorbar(im,cax=cb);cb.tick_params(labelsize=7)
            fig.suptitle(f'R15 exposed24 — all eleven archived primary selections (page {page}/2)',y=.982,fontsize=15)
            fig.text(.5,.007,'Shared color scale within each case; archived paid min-Phi primary, no new selection or scoring.',ha='center',fontsize=9)
            save(fig,f'all24_primary_page{page}','Complete24x11 primary material gallery; top=.95; row scales identical to original main gallery maxima over GT+all11 primary arrays.')
        focus=('r14_new_0009','r14_new_0020','r14_new_0021','r14_new_0023')
        focus_methods=('g19','g19_s19','twofold_hop','twofold_dr','dr','cie')
        fig=plt.figure(figsize=(19,12.5));grid=fig.add_gridspec(4,8,width_ratios=[1]*7+[.07],left=.055,right=.985,bottom=.07,top=.92,hspace=.44,wspace=.10)
        for i,cid in enumerate(focus):
            im=None
            for j,m in enumerate(('GT',)+focus_methods):
                ax=fig.add_subplot(grid[i,j]);image=show(ax,truth[cid] if m=='GT' else load_selected(cid,m,'primary'),vmax[cid])
                if image is not None:im=image
                if i==0:ax.set_title(m,fontsize=10,pad=9)
                if j==0:ax.set_ylabel(cid,fontsize=10);ax.set_xlabel('GT inverse-grid reference',fontsize=8)
                else:ax.set_xlabel(metric_text(cid,m,'primary'),fontsize=8,labelpad=6)
            cb=fig.add_subplot(grid[i,7]);fig.colorbar(im,cax=cb);cb.tick_params(labelsize=8)
        fig.suptitle('Four fixed exposed cases — material metrics and physical residual shown separately',y=.976,fontsize=15)
        fig.text(.5,.018,'NRMSE, IoU and residual/δ are copied from ALL_VIEWS; primary selectors unchanged.\nEach row shares the full eleven-arm primary gallery scale. GT is a coarse-grid material reference, not a noise certificate.',ha='center',fontsize=9)
        save(fig,'fixed4_primary_material_residual_tradeoffs','Fixed cases0009/0020/0021/0023, full IDs; GT+G19/S19/hop/cyclicDR/standaloneDR/CIE. No unified winner or new metric.')
        # Original raw/corrected comparison, including recorded raw-initial fallback.
        corrector_methods=('g19','g19_s19','twofold_hop','twofold_dr')
        corr_cases=('r14_new_0013','r14_new_0019')
        columns=[('GT',None)]+[(m,v) for m in corrector_methods for v in ('raw','corrected')]
        fig=plt.figure(figsize=(23,7.3));grid=fig.add_gridspec(2,10,width_ratios=[1]*9+[.07],left=.045,right=.985,bottom=.13,top=.84,hspace=.48,wspace=.10)
        for i,cid in enumerate(corr_cases):
            arrays=[truth[cid]]+[load_selected(cid,m,v) for m,v in columns[1:]];maximum=max(float(np.max(x)) for x in arrays if x is not None);im=None
            for j,((m,v),array) in enumerate(zip(columns,arrays)):
                ax=fig.add_subplot(grid[i,j]);image=show(ax,array,maximum)
                if image is not None:im=image
                if i==0:ax.set_title('GT' if m=='GT' else m+'\n'+v,fontsize=9,pad=8)
                if j==0:ax.set_ylabel(cid,fontsize=10);ax.set_xlabel('GT inverse-grid reference',fontsize=8)
                else:ax.set_xlabel(metric_text(cid,m,v),fontsize=8,labelpad=6)
            cb=fig.add_subplot(grid[i,9]);fig.colorbar(im,cax=cb);cb.tick_params(labelsize=8)
        fig.suptitle('Fixed corrector contrasts: archived raw and corrected material views',y=.976,fontsize=15)
        fig.text(.5,.026,'All displayed metrics are original ALL_VIEWS values; unchanged NRMSE/IoU success rule.\nHere cyclic DR means twofold_dr. Same scale across each nine-panel row; no rescoring or candidate replacement.',ha='center',fontsize=9)
        save(fig,'corrector_fixed0013_0019','GT+four methods raw/corrected; existing raw-source/initial fallback and corrected behavior preserved. Continuous metrics copied unchanged.')
        comparison_path=root/'audit/readout_counterfactual_v1/COMPARISON.json';comparison=base.read(comparison_path);record['sources'].append(str(comparison_path))
        if len(comparison)!=16:record['issues'].append(dict(issue='EXPECTED16_COUNTERFACTUAL_ROWS',actual=len(comparison)))
        x=np.arange(len(comparison));linear=[r.get('linear_predicted_delta_Phi',np.nan) for r in comparison];actual=[r.get('actual_delta_Phi',np.nan) for r in comparison]
        cos=[r.get('cosine_to_negative_physical_gradient',np.nan) for r in comparison]
        linear=[v if base.finite(v) else np.nan for v in linear];actual=[v if base.finite(v) else np.nan for v in actual];cos=[v if base.finite(v) else np.nan for v in cos]
        labels=[r['case_id']+'\n'+r['readout'] for r in comparison]
        fig,axes=plt.subplots(2,1,figsize=(19,10.5),sharex=True,gridspec_kw={'height_ratios':[1.3,1]})
        axes[0].bar(x-.18,linear,width=.36,label='linear predicted ΔPhi = grad(Phi) · Δm',color='#2166ac')
        axes[0].bar(x+.18,actual,width=.36,label='actual archived ΔPhi',color='#b35806');axes[0].axhline(0,color='black',lw=.7);axes[0].legend(fontsize=10);axes[0].set_ylabel('Signed ΔPhi');axes[0].grid(axis='y',alpha=.2)
        axes[1].scatter(x,cos,color='#1b7837',s=34);axes[1].axhline(0,color='gray',lw=.7);axes[1].set_ylim(-1.08,1.08);axes[1].set_ylabel('cos(Δm, −grad physical Phi)');axes[1].grid(axis='y',alpha=.2)
        axes[1].set_xticks(x,labels,rotation=55,ha='right',fontsize=8)
        fig.suptitle('All sixteen fixed-state readout counterfactuals: signed linear prediction, actual change and gradient cosine',fontsize=14)
        fig.subplots_adjust(top=.92,bottom=.23,hspace=.16)
        fig.text(.5,.022,'OFFLINE_DIAGNOSTIC; fixed checkpoint J/Q/PR retained. Rows and actual stage/beta come from COMPARISON.json.\nNo GT selection, no new readout/forward evaluation; null cosine remains undefined.',ha='center',fontsize=9)
        save(fig,'readout16_linear_actual_phi_gradient_cosine','All16 saved counterfactual rows, signed Phi changes; cosine to negative original physical gradient; no GT-based selection.')
        # Preserve exactly the three original remote PNGs as archived presentations.
        dynamics=root/'FIGURES/dynamics';dynamics.mkdir(parents=True,exist_ok=True)
        for name in ('samecheckpoint_phi_gap_material.png','finite_response_gains_box_stage.png','cycle_material_bands_projection_loss.png'):
            source=root/'analysis/fullstate_v1/figures'/name;target=dynamics/('original_remote_'+name)
            shutil.copyfile(source,target)
            record['artifacts'].append(dict(path=str(target.relative_to(root)),bytes=target.stat().st_size,source=str(source),kind='Original remote PNG preserved unchanged; third retains known legend bug for provenance'))
        cycle_path=root/'analysis/fullstate_v1/MAIN_CYCLES.json';cycles=base.read(cycle_path);record['sources'].append(str(cycle_path))
        methods=('twofold_hop','twofold_dr','twofold_ap','twofold_full_dr','twofold_reset');bands=('U11','U10','U01','U00');colors=['#2166ac','#e08214','#1b7837','#762a83']
        fig,axes=plt.subplots(1,3,figsize=(20,7.2));undefined={}
        for mi,method in enumerate(methods):
            group=[r for r in cycles if r['method']==method and r.get('available')]
            undefined[method]=sum(not base.finite(r.get('U11_material_fraction')) for r in group)
            for r in group:
                if base.finite(r.get('material_outside_real_U01_fraction')):axes[0].scatter(mi,r['material_outside_real_U01_fraction'],s=12,alpha=.45,color='#2166ac')
                if base.finite(r.get('current_projection_relative_loss')):axes[1].scatter(mi,r['current_projection_relative_loss'],s=12,alpha=.45,color='#b35806')
            for bi,band in enumerate(bands):
                values=[r[band+'_material_fraction'] for r in group if base.finite(r.get(band+'_material_fraction'))]
                if values:axes[2].scatter([mi+.13*(bi-1.5)]*len(values),values,s=12,alpha=.4,color=colors[bi])
            if group and undefined[method]==len(group):
                for ax in (axes[0],axes[2]):ax.text(mi,.97,'Δm=0\nundefined',ha='center',va='top',fontsize=8,color='#555')
        # Explicit legend handles remain correct when the first method has no Δm.
        fig.legend(handles=[Line2D([0],[0],marker='o',linestyle='none',color=c,label=b) for c,b in zip(colors,bands)],loc='upper center',bbox_to_anchor=(.835,.935),ncol=4,frameon=False,fontsize=9)
        for ax in axes:
            ax.set_xlim(-.45,4.45)
            ax.set_xticks(range(len(methods)),methods,rotation=32,ha='right',fontsize=9);ax.grid(axis='y',alpha=.2)
        axes[0].set_ylim(-.03,1.04);axes[2].set_ylim(-.03,1.04)
        axes[0].set_title('Cycle Δm energy outside real U01');axes[0].set_ylabel('Fraction; zero Δm is undefined')
        axes[1].set_title('Saved complex-current injection projection loss');axes[1].set_ylabel('Relative norm')
        axes[2].set_title('Cycle Δm energy by real-material sector');axes[2].set_ylabel('Fraction; no zero-filling')
        fig.suptitle('First four archived cycle windows only — actual stage0, beta(current,material)=(1,1)',fontsize=14,y=.98)
        fig.subplots_adjust(top=.87,bottom=.27,wspace=.27)
        caption='Real-material band fractions describe Δm; routing H is a complex-current subspace and is not a material-step constraint.\nHop/reset zero Δm fractions stay undefined. These early four beta=(1,1) windows do not establish later-stage behavior.'
        fig.text(.5,.025,caption,ha='center',fontsize=10)
        save(fig,'cycle_material_bands_projection_loss_fixed',caption)
        record['cycle_fraction_undefined_counts']=undefined
        record['cycle_rows']=len(cycles)
        record['cycle_actual_before_stage_metadata']=list({json.dumps(r.get('before_stage'),sort_keys=True) for r in cycles})
        record['selected_material_sources']=selected_sources
        record['status']='COMPLETE_WITH_RETAINED_ISSUES' if record['issues'] else 'COMPLETE'
    except BaseException as exc:
        record.update(status='FAILED',error=repr(exc));raise
    finally:
        record.update(CPU_seconds=time.process_time()-START_CPU,wall_seconds=time.perf_counter()-START_WALL,whole_process_CPU_seconds_at_receipt=time.process_time(),finished_utc=datetime.now(timezone.utc).isoformat())
        (out/'RECEIPT.json').write_text(json.dumps(record,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
        with (out/'RENDER_ATTEMPTS.jsonl').open('a') as f:f.write(json.dumps(record,ensure_ascii=False,allow_nan=False)+'\n')
    print(json.dumps({k:record[k] for k in ('status','CPU_seconds','wall_seconds')},indent=2))


if __name__=='__main__':main()
