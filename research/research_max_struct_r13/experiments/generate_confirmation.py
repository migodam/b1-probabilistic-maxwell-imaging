"""Generate independent R13 public observations; never optimize or score reconstructions.

Material seeds/families, clean observations and cell-average truth go exclusively
under sealed_truth/fresh48. Public inverse inputs contain all576unknowns.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import argparse,json,sys,time,traceback
from pathlib import Path


def main():
    wall0,cpu0=time.perf_counter(),time.process_time()
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);args=ap.parse_args()
    root=args.root.resolve();repo=root.parents[1];sys.path.insert(0,str(repo))
    import numpy as np
    from scipy.stats import chi2
    from research.research_fbg_r3_experiments.src.physics import make_grid, MaxwellModel
    publicdir=root/'data_public/fresh48';sealdir=root/'sealed_truth/fresh48';audit=root/'delegated/evidence/generation'
    audit.mkdir(parents=True,exist_ok=True)
    if publicdir.exists() or sealdir.exists():raise FileExistsError('Fresh48 destination exists; never silently overwrite sealed observations/truth.')
    publicdir.mkdir(parents=True);sealdir.mkdir(parents=True)
    def dump(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
    publicmanifest={'schema_version':'r13_fresh48_v1','status':'GENERATING','material_parameter':'real_scalar_contrast_chi','unknowns_per_case':576,'inverse_grid':[24,24],'roi_width_wavelengths':.96,'k':float(2*np.pi),'lower':.02,'upper':20.,'initial_material':.1,'noise_model':'independent circular complex Gaussian; real and imag variance sigma²/2; sigma=relative_RMS*RMS(clean observation)','delta_definition':'sigma*sqrt(chi2_ppf(.99,2*N_complex)/2); noise-only radius, does not include cross-grid model mismatch','case_files':[]}
    families=['multi_island_ellipses','continuous_random_field','random_level_patchwork','perforated_slab_closed_holes','open_notched_porous_body','thin_bridge_network','eccentric_multilayers','mixed_scale_disconnected']
    # The material-family assignment is deliberately absent from public records.
    assignments=[(f,v) for f in range(8) for v in range(6)]
    permutation=np.random.default_rng(131299).permutation(48)
    assignment={int(permutation[j]):fv for j,fv in enumerate(assignments)}
    privatecases=[]
    for idx in range(48):
        f,v=assignment[idx]
        cross=(v==0 or (f%2==0 and v==3)) # 8+4 objects; no duplicated-noise copies.
        privatecases.append({'case_id':f'r13_fresh_{idx:04d}','material_seed':131300+idx,'family':families[f],'family_index':f,'family_variant':v,'forward_side':36 if cross else 24,'inverse_side':24,'acquisition':['full12','full6','half12'][v//2],'relative_noise_rms':[.01,.02][v%2],'noise_seed':231300+idx,'quadrature_reference_side':216})
    private={'schema_version':'r13_fresh48_generation_v1','status':'PREDECLARED_BEFORE_GENERATION','family_count':8,'objects_per_family':6,'objects':48,'cross_grid_objects':12,'material_seed_rule':'131300+opaque_index','assignment_permutation_seed':131299,'truth_semantics':'square-cell contrast average of one common216² continuous-geometry quadrature; disk-pulse VIE remains inherited model','test_truth_policy':'generator-only; solver loads only data_public; score only after predictions sealed','cases':privatecases}
    dump(sealdir/'GENERATOR_MANIFEST.json',private);dump(publicdir/'MANIFEST.json',publicmanifest)
    # Only two physical grids are assembled; acquisition variants select actual
    # inherited source/receiver rows/columns rather than inventing new kernels.
    setup=[];grids={}
    for side in [24,36]:
        w,c=time.perf_counter(),time.process_time();p=make_grid(side,sources=12,receivers=48,k=2*np.pi)
        grids[side]=p;setup.append({'category':'operator_setup','side':side,'Green_setups':1,'LU':0,'complex_RHS':0,'cpu_seconds':time.process_time()-c,'wall_seconds':time.perf_counter()-w})
    axis=(np.arange(216)+.5)*.96/216-.48;yy,xx=np.meshgrid(axis,axis,indexing='ij')
    def ellipse(cx,cy,rx,ry,angle=0.):
        co,si=np.cos(angle),np.sin(angle);x=co*(xx-cx)+si*(yy-cy);y=-si*(xx-cx)+co*(yy-cy)
        return (x/rx)**2+(y/ry)**2<=1
    def segment(ax,ay,bx,by,width):
        dx,dy=bx-ax,by-ay;t=np.clip(((xx-ax)*dx+(yy-ay)*dy)/(dx*dx+dy*dy),0,1)
        return (xx-ax-t*dx)**2+(yy-ay-t*dy)**2<=width*width
    def material(case):
        rng=np.random.default_rng(case['material_seed']);f=case['family_index'];v=case['family_variant']
        high=float(rng.uniform(12,18));low=float(rng.uniform(2,5));angle=float(rng.uniform(-np.pi,np.pi));shift=rng.uniform(-.05,.05,2)
        h=low if v==0 else high;m=np.full(xx.shape,.1);detail={'high':high,'low':low,'rotation_radians':angle,'shift':shift.tolist()}
        def put(mask,val):m[mask]=val
        if f==0:
            count=int(rng.integers(3,7));objs=[]
            for j in range(count):
                cx,cy=rng.uniform(-.29,.29,2);rx,ry=rng.uniform(.055,.13,2);a=float(rng.uniform(-np.pi,np.pi));val=h if j%2==0 else low
                put(ellipse(cx,cy,rx,ry,a),val);objs.append([cx,cy,rx,ry,a,val])
            detail['ellipses']=objs
        elif f==1:
            z=np.zeros(xx.shape);bumps=[]
            for j in range(8):
                cx,cy=rng.uniform(-.3,.3,2);width=float(rng.uniform(.07,.19));amp=float(rng.uniform(-1,1));z+=amp*np.exp(-((xx-cx)**2+(yy-cy)**2)/(2*width*width));bumps.append([cx,cy,width,amp])
            z=(z-z.min())/(z.max()-z.min());taper=np.clip((.40-np.sqrt((xx-shift[0])**2+(yy-shift[1])**2))/.09,0,1)
            m=.1+taper*(low+(h-low)*z-.1);detail['bumps']=bumps
        elif f==2:
            centers=rng.uniform(-.32,.32,(7,2));d=(xx[:,:,None]-centers[:,0])**2+(yy[:,:,None]-centers[:,1])**2;label=np.argmin(d,axis=2)
            vals=rng.uniform(low,h,7);vals[::3]=h;mask=ellipse(*shift,.37,.29,angle)
            put(mask,vals[label][mask]);detail.update(centers=centers.tolist(),levels=vals.tolist())
        elif f==3:
            co,si=np.cos(angle),np.sin(angle);x=co*(xx-shift[0])+si*(yy-shift[1]);y=-si*(xx-shift[0])+co*(yy-shift[1]);body=(abs(x)<.34)&(abs(y)<.27)
            put(body,h);holes=[]
            for cx,cy,r in [(-.15,-.1,.065),(.12,.11,.085),(.16,-.12,.045)]:
                cx+=float(rng.uniform(-.025,.025));cy+=float(rng.uniform(-.025,.025));hole=(x-cx)**2+(y-cy)**2<r*r;put(hole,.1);holes.append([cx,cy,r])
            put(body&(x<-.25),low);detail['holes_in_body_coordinates']=holes
        elif f==4:
            co,si=np.cos(angle),np.sin(angle);x=co*(xx-shift[0])+si*(yy-shift[1]);y=-si*(xx-shift[0])+co*(yy-shift[1]);body=(abs(x)<.34)&(abs(y)<.29)
            # A C-shaped rectangular notch and a separate closed oval hole.
            notch=(x>-.10)&(abs(y)<rng.uniform(.065,.10));closed=((x+.21)/.055)**2+((y-.14)/.06)**2<1
            put(body&~notch&~closed,h);put(body&(y<-.19),low)
        elif f==5:
            pts=np.array([[-.27,-.21],[.25,-.18],[-.03,.26]])+shift
            widths=[float(rng.uniform(.018,.033)),float(rng.uniform(.025,.045))]
            for j in [0,1]:put(segment(*pts[j],*pts[2],widths[j]),low if j==0 else h)
            for j,(cx,cy) in enumerate(pts):put(ellipse(cx,cy,float(rng.uniform(.07,.105)),float(rng.uniform(.055,.095)),angle+j),h if j!=1 else low)
            detail.update(junctions=pts.tolist(),bridge_halfwidths=widths)
        elif f==6:
            put(ellipse(*shift,.35,.285,angle),h)
            center=shift+rng.uniform(-.06,.06,2);put(ellipse(*center,.25,.18,angle+.3),low)
            inner=center+rng.uniform(-.035,.035,2);put(ellipse(*inner,.13,.09,angle-.4),h*.65)
            put(ellipse(*(inner+np.array([.02,-.01])),.055,.04,angle),.1)
            detail.update(layer_centers=[shift.tolist(),center.tolist(),inner.tolist()])
        else:
            centers=np.array([[-.24,-.20],[.19,.20],[-.21,.21],[.22,-.21]])+rng.uniform(-.035,.035,(4,2))
            for j,(cx,cy) in enumerate(centers):
                rx,ry=(.12,.085) if j<2 else(.05,.075);mask=ellipse(cx,cy,rx,ry,angle+j*.4)
                val=low if j%2 else h
                # One smooth internal gradient gives a mixed piecewise/continuous object.
                if j==1:put(mask,(val+2*(xx-cx))[mask])
                else:put(mask,val)
            detail['island_centers']=centers.tolist()
        return np.clip(m,.02,20),detail
    def average(a,side):
        q=216//side;return a.reshape(side,q,side,q).mean(axis=(1,3)).ravel()
    def acquisition(side,kind):
        p=grids[side];si=np.arange(0,12,2) if kind=='full6' else np.arange(12);ri=np.arange(24) if kind=='half12' else np.arange(48)
        return dict(G=p['G'],E0=p['E0'][:,si],S=p['S'][ri],coords=p['centers'],centers=p['centers'],sources=p['sources'][si],receivers=p['receivers'][ri],grid_shape=np.array([side,side]),k=float(p['k']),radius=float(p['radius']),cell_radius=float(p['radius'])),si,ri
    costs=[];sealedmetrics=[]
    try:
        for case in privatecases:
            c0,w0=time.process_time(),time.perf_counter();cid=case['case_id'];ref,params=material(case);coarse=average(ref,24);fine=average(ref,case['forward_side'])
            forward,si,ri=acquisition(case['forward_side'],case['acquisition']);model=MaxwellModel(forward);state=model.state(fine);clean=state['Y'].copy();counts=dict(model.cost)
            sigma=case['relative_noise_rms']*float(np.linalg.norm(clean)/np.sqrt(clean.size))
            noise_rng=np.random.default_rng(case['noise_seed']);noise=sigma/np.sqrt(2)*(noise_rng.normal(size=clean.shape)+1j*noise_rng.normal(size=clean.shape));Y=clean+noise
            delta=float(sigma*np.sqrt(chi2.ppf(.99,2*Y.size)/2));inverse,_,_=acquisition(24,case['acquisition'])
            mismatch=None;coarse_prediction=clean
            if case['forward_side']!=24:
                cm=MaxwellModel(inverse);coarse_state=cm.state(coarse);coarse_prediction=coarse_state['Y'];mismatch=float(np.linalg.norm(coarse_prediction-clean)/np.linalg.norm(clean))
                for k,val in cm.cost.items():counts[k]=counts.get(k,0)+val
            # Full individual input files match run_experiment.py's existing API.
            np.savez_compressed(publicdir/(cid+'.npz'),**inverse,Y=Y,sigma=sigma,delta=delta,m_init=np.full(576,.1),components=1,lower=.02,upper=20.,inverse_side=24,forward_side=case['forward_side'],noise_relative_rms=case['relative_noise_rms'])
            np.savez_compressed(sealdir/(cid+'.npz'),material_coarse=coarse,material_forward=fine,Y_clean=clean,Y_coarse_truth=coarse_prediction,noise=noise,forward_side=case['forward_side'])
            params['case_id']=cid;dump(sealdir/(cid+'_GEOMETRY.json'),params)
            sealedmetrics.append({'case_id':cid,'crossgrid':case['forward_side']!=24,'clean_coarse_forward_relative_mismatch':mismatch,'realized_noise_over_delta':float(np.linalg.norm(noise)/delta),'generation_forward_relative_residual':float(state['state_relative_residual'])})
            rec={'case_id':cid,'input_file':str((publicdir/(cid+'.npz')).relative_to(root)),'inverse_side':24,'forward_side':case['forward_side'],'sources':len(si),'receivers':len(ri),'receiver_aperture':'half_circle' if case['acquisition']=='half12' else'full_circle','relative_noise_RMS':case['relative_noise_rms'],'sigma':sigma,'delta_noise_only':delta,'model_mismatch_status':'SEPARATE_CROSSGRID_DIAGNOSTIC_NOT_INCLUDED_IN_NOISE_RADIUS' if case['forward_side']!=24 else'MATCHED_DISCRETE_MODEL','material_truth_in_public':False}
            publicmanifest['case_files'].append(rec)
            costs.append({'case_id':cid,'counts':counts,'cpu_seconds_including_serialization':time.process_time()-c0,'wall_seconds_including_serialization':time.perf_counter()-w0,'status':'GENERATED','solver_run':False})
            dump(audit/'CASE_GENERATION_COSTS.json',costs);dump(publicdir/'MANIFEST.json',publicmanifest)
            print(cid,'generated',flush=True)
        private['status']='GENERATION_COMPLETE_TRUTH_SEALED';dump(sealdir/'GENERATOR_MANIFEST.json',private);dump(sealdir/'GENERATION_DIAGNOSTICS.json',sealedmetrics)
        publicmanifest['status']='READY_FOR_FROZEN_METHOD_DISPATCH';dump(publicdir/'MANIFEST.json',publicmanifest)
        totals={'LU':sum(r['counts'].get('LU',0) for r in costs),'complex_RHS':sum(r['counts'].get('actual_complex_RHS',0) for r in costs),'setup':setup,'process_CPU_seconds':time.process_time()-cpu0,'process_wall_seconds':time.perf_counter()-wall0,'case_CPU_sum':sum(r['cpu_seconds_including_serialization'] for r in costs),'case_wall_sum':sum(r['wall_seconds_including_serialization'] for r in costs),'memory_peak_bytes':None,'status':'READY','objects':48,'cross_grid_objects':12,'case_manifest_public':str(publicdir/'MANIFEST.json'),'case_manifest_sealed':str(sealdir/'GENERATOR_MANIFEST.json'),'noise_only_delta_not_model_error_certificate':True,'inversion_or_reconstruction_scoring_performed':False}
        dump(audit/'GENERATION_RECEIPT.json',totals);print(json.dumps(totals),flush=True)
    except BaseException:
        dump(audit/'GENERATION_FAILURE.json',{'traceback':traceback.format_exc(),'process_CPU_seconds':time.process_time()-cpu0,'process_wall_seconds':time.perf_counter()-wall0,'completed_objects':len(costs),'setup':setup,'last_case':case,'completed_costs':costs});raise
if __name__=='__main__':main()
