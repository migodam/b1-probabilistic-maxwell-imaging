"""Render every sealed physics_v1 vector case/method as exact coarse voxels.

Parent-authorized terminal physics seal is validated BEFORE NumPy, selected
materials or GT are loaded. No forward solve, interpolation, surface smoothing,
new support threshold, result selection or additional experiment is performed.
Each page has at most four methods plus GT, a true3D voxel row, and three exact
orthogonal central-slice rows with the same physical chi color scale [0,8].
"""
from __future__ import annotations
import os
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[_key]='1'
import argparse,csv,json,time,traceback
from pathlib import Path
from datetime import datetime,timezone
try:
    from .seal_results import read,validate_seal
except ImportError:
    from seal_results import read,validate_seal


def utc():return datetime.now(timezone.utc).isoformat()


def dump(path,value):
    def encode(x):
        if isinstance(x,Path):return str(x)
        if hasattr(x,'item'):return x.item()
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(value,indent=2,default=encode,allow_nan=False)+'\n')


def table(out,name,rows):
    dump(out/(name+'.json'),rows)
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with (out/(name+'.csv')).open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=keys);writer.writeheader()
        for row in rows:
            writer.writerow({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in row.items()})


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--seal',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--methods-per-page',type=int,default=4,choices=(1,2,3,4))
    p.add_argument('--dpi',type=int,default=160)
    a=p.parse_args();root=a.root.resolve();out=a.out.resolve()
    if not 100<=a.dpi<=300:raise ValueError('dpi must be between100 and300')
    start,wall=time.process_time(),time.perf_counter()
    # All private reads and numerical imports are strictly below this guard.
    seal=read(a.seal);validate_seal(seal,require_authorization=True)
    if Path(seal['root']).resolve()!=root:raise RuntimeError('Reviewed seal root differs from requested root')
    cfg=seal['manifest_snapshot']
    if cfg.get('campaign')!='physics_v1' or seal.get('campaign')!='physics_v1':
        raise RuntimeError('Only the original authorized physics_v1 manifest is accepted')
    if len(cfg['tasks'])!=60:raise RuntimeError('Expected original complete60-task physics manifest')
    if {t['task_id'] for t in cfg['tasks']}!={t['task_id'] for t in seal['task_accounting']}:
        raise RuntimeError('Reviewed seal does not account for every original task')
    tasks=[t for t in cfg['tasks'] if t.get('family')=='vector3d']
    expected_cases={'slotted_high_frame','unequal_coupled_steps','mixed_crossed_bracket','slotted_high_frame_crossgrid8'}
    expected_methods={'cie','dual_cie','global_dual','dr','ap','event_integer','direct'}
    if len(tasks)!=28 or {t['case_id'] for t in tasks}!=expected_cases:
        raise RuntimeError('Original four-vector-condition manifest is required, with no case filtering')
    cases={}
    for task in tasks:cases.setdefault(task['case_id'],[]).append(task)
    for cid,case in cases.items():
        if len(case)!=7 or {t['method'] for t in case}!=expected_methods:
            raise RuntimeError('All seven original methods are required for every vector condition')
    if any(out.is_relative_to(Path(d['path']).resolve()) for d in seal['directories']):
        raise ValueError('Plot output must be outside sealed result and queue directories')
    out.mkdir(parents=True,exist_ok=False);figdir=out/'FIGURES';figdir.mkdir()
    authorization_time=utc();figures=[];panels=[];issues=[];access=[]
    attempts=0;imports_CPU=None;imports_wall=None;fatal_error=None
    try:
        imp_cpu,imp_wall=time.process_time(),time.perf_counter()
        import numpy as np
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.colors import Normalize
        from matplotlib.cm import ScalarMappable
        imports_CPU=time.process_time()-imp_cpu;imports_wall=time.perf_counter()-imp_wall
        side=6;width=.24;support_threshold=1.;mid=3
        norm=Normalize(vmin=0.,vmax=8.,clip=True);cmap=plt.get_cmap('viridis')
        camera=dict(elevation=25.,azimuth=-55.);views=('initial','raw','primary','corrected')
        coordinate_edges=np.linspace(-width/2,width/2,side+1)
        xyz=np.meshgrid(coordinate_edges,coordinate_edges,coordinate_edges,indexing='ij')
        plane_coordinate=-width/2+(mid+.5)*width/side
        def material(path,key):
            if not path.is_file():return None,'MISSING_SAVED_ARRAY'
            try:
                with np.load(path,allow_pickle=False) as z:value=z[key].copy()
                if value.size!=216 or np.iscomplexobj(value) or not np.isfinite(value).all():
                    return None,'INVALID_REAL_6_CUBED_ARRAY'
                return value.reshape(side,side,side),'AVAILABLE'
            except Exception as exc:return None,repr(exc)
        def axis3(ax):
            ax.set(xlim=(-.12,.12),ylim=(-.12,.12),zlim=(-.12,.12))
            ax.set_box_aspect((1,1,1));ax.view_init(elev=camera['elevation'],azim=camera['azimuth'])
            ax.set_xticks([-.12,0,.12]);ax.set_yticks([-.12,0,.12]);ax.set_zticks([-.12,0,.12])
            ax.set_xlabel('x / lambda0',fontsize=7,labelpad=0);ax.set_ylabel('y / lambda0',fontsize=7,labelpad=0)
            ax.set_zlabel('z / lambda0',fontsize=7,labelpad=0);ax.tick_params(labelsize=6,pad=0)
        for cid,case_tasks in cases.items():
            source_case=case_tasks[0]
            with np.load(root/source_case['input_file'].replace('\\','/'),allow_pickle=False) as public:
                if tuple(map(int,public['grid_shape']))!=(6,6,6) or float(public['width'])!=.24:
                    raise ValueError('Original exact6-cubed/.24lambda vector grid required')
                if float(public['lower'])!=0 or float(public['upper'])!=8:
                    raise ValueError('Original shared physical contrast display/box [0,8] required')
            spec=source_case.get('truth_spec',cfg.get('truth_specs',{}).get(cid,{}))
            if not spec.get('path') or spec.get('key')!='chi':
                raise RuntimeError('Original coarse vector truth path/key missing from reviewed manifest')
            gtpath=root/spec['path'];access.append(dict(kind='offline_GT',case_id=cid,
                path=str(gtpath.resolve()),access_at_utc=utc(),authorized_at_utc=authorization_time))
            gt,gt_status=material(gtpath,spec['key'])
            if gt is None:issues.append(dict(case_id=cid,role='GT',error=gt_status))
            loaded={};metadata={};taskmap={t['method']:t for t in case_tasks}
            methods=[t['method'] for t in case_tasks]
            for method in methods:
                folder=root/'RESULTS'/cfg['campaign']/taskmap[method]['task_id']
                receipt=read(folder/'RECEIPT.json') if (folder/'RECEIPT.json').is_file() else{}
                for view in views:
                    path=folder/(view+'.npz');volume,status=material(path,'m')
                    loaded[(method,view)]=(volume,status,path)
                    saved=read(folder/(view+'.json')) if (folder/(view+'.json')).is_file() else{}
                    metadata[(method,view)]=dict(saved_selection_metadata=saved,run_status=receipt.get('status','MISSING_RECEIPT'))
                    if volume is None:issues.append(dict(case_id=cid,method=method,view=view,error=status))
            cross=bool(spec.get('crossgrid_pair_not_independent',False))
            physical_label='cross-grid8^3 -> inverse6^3' if cross else'matched6^3 generation/inversion'
            for view in views:
                for page,offset in enumerate(range(0,len(methods),a.methods_per_page),1):
                    subset=methods[offset:offset+a.methods_per_page]
                    columns=([('GT',gt,gt_status,gtpath)]+
                        [(method,*loaded[(method,view)]) for method in subset])
                    ident=f'V3D__{cid}__{view}__page{page:02d}'
                    png=figdir/(ident+'.png');pdf=figdir/(ident+'.pdf')
                    attempts+=1;cpu0,wall0=time.process_time(),time.perf_counter();fig=None
                    record=dict(figure_id=ident,case_id=cid,object_id=spec.get('object_id',cid),
                        view=view,methods=subset,page=page,maximum_methods_per_page=a.methods_per_page,
                        png_download_path=str(png.relative_to(out)),pdf_download_path=str(pdf.relative_to(out)),
                        png_absolute_path=str(png),pdf_absolute_path=str(pdf),
                        source_manifest=str(seal['manifest_path']),physical_condition=physical_label,
                        camera=camera,contrast_range=[0,8],support_rule='chi > 1',
                        exact_material_grid=[6,6,6],interpolation='NONE',slice_index=mid,
                        slice_plane_coordinate_lambda=plane_coordinate)
                    try:
                        fig=plt.figure(figsize=(3.15*len(columns),11.7))
                        grid=fig.add_gridspec(4,len(columns),height_ratios=[1.5,1,1,1])
                        for col,(label,volume,status,path) in enumerate(columns):
                            title='Offline coarse GT' if label=='GT' else label
                            if label=='global_dual':title+='\nglobal-PR Dual-CIE'
                            ax=fig.add_subplot(grid[0,col],projection='3d');axis3(ax)
                            panel=dict(figure_id=ident,case_id=cid,method=label,view='GT' if label=='GT' else view,
                                source_array=str(path.resolve()),availability=status,
                                voxel_support_rule='chi > 1',voxel_color_range=[0,8],material_resampled=False)
                            if label!='GT':panel.update(metadata[(label,view)])
                            if volume is not None:
                                filled=volume>support_threshold
                                ax.voxels(*xyz,filled,facecolors=cmap(norm(volume)),edgecolor=(.12,.12,.12,.65),
                                    linewidth=.23,shade=False)
                                panel.update(material_min=float(volume.min()),material_max=float(volume.max()),
                                    support_cells=int(filled.sum()),below_box_cells=int((volume<0).sum()),
                                    above_box_cells=int((volume>8).sum()))
                                title+=f'\nchi [{volume.min():.3g},{volume.max():.3g}]; support {filled.sum()}/216'
                                if not filled.any():ax.text2D(.5,.43,'EMPTY chi>1 SUPPORT',transform=ax.transAxes,ha='center',fontsize=7)
                                if (volume<0).any() or (volume>8).any():title+='\nOUT-OF-BOX CELLS; fixed color saturation'
                            else:
                                ax.text2D(.5,.45,'UNAVAILABLE\n'+status,transform=ax.transAxes,ha='center',fontsize=7,wrap=True)
                            ax.set_title(title,fontsize=8,pad=2)
                            slices=(None,None,None) if volume is None else(volume[mid,:,:].T,volume[:,mid,:].T,volume[:,:,mid].T)
                            for row,(image,plane,axes) in enumerate(zip(slices,('YZ: x','XZ: y','XY: z'),(('y','z'),('x','z'),('x','y'))),1):
                                aa=fig.add_subplot(grid[row,col]);aa.set_aspect('equal')
                                aa.set_title(f'{plane}={plane_coordinate:+.3f} lambda0; exact index {mid}',fontsize=7)
                                if image is None:
                                    aa.text(.5,.5,'UNAVAILABLE',transform=aa.transAxes,ha='center',fontsize=8)
                                    aa.set(xlim=(-.12,.12),ylim=(-.12,.12))
                                else:aa.imshow(image,origin='lower',extent=(-.12,.12,-.12,.12),
                                    cmap=cmap,norm=norm,interpolation='nearest')
                                aa.set_xticks([-.12,0,.12]);aa.set_yticks([-.12,0,.12]);aa.tick_params(labelsize=6)
                                aa.set_xlabel(axes[0]+' / lambda0',fontsize=7);aa.set_ylabel(axes[1]+' / lambda0',fontsize=7)
                            panels.append(panel)
                        fig.suptitle(cid+' | '+view+' | '+physical_label+
                            '\nSynthetic full-vector CM+RR DDA; finite coarse-grid evidence, continuum accuracy unqualified',fontsize=11,y=.985)
                        fig.subplots_adjust(left=.05,right=.90,bottom=.055,top=.88,hspace=.6,wspace=.34)
                        cax=fig.add_axes([.925,.18,.014,.63]);fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=cax,
                            label='physical susceptibility chi; fixed [0,8]')
                        fig.text(.5,.014,'Exact6^3 cells, no interpolation or refined-looking surface; chi>1 voxels show support; exact central sections show interior cuts. Same camera/scale for every panel.',
                            ha='center',fontsize=7,wrap=True)
                        fig.savefig(png,dpi=a.dpi);fig.savefig(pdf)
                        record.update(status='COMPLETE',png_bytes=png.stat().st_size,pdf_bytes=pdf.stat().st_size)
                    except Exception as exc:
                        record.update(status='PLOT_FAILED',error=repr(exc));issues.append(dict(figure_id=ident,error=repr(exc)))
                    finally:
                        if fig is not None:plt.close(fig)
                        record.update(CPU_seconds=time.process_time()-cpu0,wall_seconds=time.perf_counter()-wall0)
                        figures.append(record)
        validate_seal(seal,require_authorization=True)
    except BaseException as exc:
        fatal_error=repr(exc);issues.append(dict(stage='FATAL',error=fatal_error,traceback=traceback.format_exc()))
    finally:
        table(out,'FIGURE_MANIFEST',figures);table(out,'FIGURE_PANEL_TABLE',panels)
        dump(out/'ISSUES.json',issues);dump(out/'PRIVATE_ACCESS_LOG.json',access)
        dump(out/'PLOT_COST_RECEIPT.json',dict(status='FAILED' if fatal_error else'COMPLETE_WITH_MISSING_OR_FAILED_PANELS' if issues else'COMPLETE',
            created_at_utc=utc(),authorized_seal=str(a.seal.resolve()),authorized_at_utc=authorization_time,
            campaign='physics_v1',vector_conditions=4,independent_material_objects=3,methods_per_condition=7,
            requested_views=['initial','raw','primary','corrected'],expected_method_view_panels=112,
            expected_figures=4*4*((7+a.methods_per_page-1)//a.methods_per_page),
            current_run_figure_attempts=attempts,completed_figures=sum(f['status']=='COMPLETE' for f in figures),
            current_run_failed_figure_attempts=sum(f['status']!='COMPLETE' for f in figures),
            prior_failed_plot_attempt_count=None,prior_attempt_scope='unknown; never inferred as zero',
            import_CPU_seconds=imports_CPU,import_wall_seconds=imports_wall,
            figure_CPU_seconds_sum=sum(f['CPU_seconds'] for f in figures),
            figure_wall_seconds_sum=sum(f['wall_seconds'] for f in figures),
            total_plot_CPU_seconds=time.process_time()-start,total_plot_wall_seconds=time.perf_counter()-wall,
            cost_scope='post-seal plot/input-load/export only; figure/import times overlap total; no online cost attribution',
            new_forward_solves=0,new_LU=0,new_RHS=0,material_grid=[6,6,6],resampling=False,
            material_truth_access='authorized original coarse chi only; no generator geometry or heldout measured archive',
            scientific_scope='visual display of all saved conditions/methods/views; no new scientific acceptance or continuum claim',
            fatal_error=fatal_error))
    if fatal_error:raise RuntimeError(fatal_error)
    print(json.dumps(dict(status='VECTOR_VOLUME_GALLERY_WRITTEN',figures=len(figures),output=str(out))))


if __name__=='__main__':main()
