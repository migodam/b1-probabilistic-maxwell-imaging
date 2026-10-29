"""Independent R3 offline generator. Public inputs never contain clean data or truth."""
from pathlib import Path
import argparse, json, sys, time
import os
for _thread_key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS']:os.environ[_thread_key]='1'
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'frozen/src' if (ROOT/'frozen/src').exists() else ROOT/'src'))
CASES={'dev':33102901,'blind_normal':33102911,'blind_cavity':33102912,'blind_strong':33102913}
NAMES=[*[f'wall_eps{i}' for i in range(6)],*[f'wall_loss{i}' for i in range(6)],'target_x','target_y','target_rx','target_ry','target_contrast','target_loss','target_angle']
def plain(a):
    if hasattr(a,'tolist'):return a.tolist()
    if hasattr(a,'item'):return a.item()
    if isinstance(a,Path):return str(a)
    raise TypeError(type(a).__name__)
def write_json(path,value):
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value,indent=2,default=plain,allow_nan=False)+'\n')
def sample(case):
    """Only frozen independent draws; shape strata are explicit conditioning."""
    r=np.random.default_rng(CASES[case]);e=r.uniform(2,4.5,6);l=r.uniform(.02,.12,6)
    x=r.uniform(-.8,.8);y=r.uniform(1.8,2.4);rx=r.uniform(.085,.155);ry=r.uniform(.085,.155)
    contrast=r.uniform(6.5,8) if case=='blind_strong' else r.uniform(2,5)
    loss=r.uniform(.005,.20);angle=r.uniform(-np.pi/2,np.pi/2)
    if case in ('dev','blind_normal'):ry=rx;angle=0.
    cavity=int(r.integers(6)) if case=='blind_cavity' else None
    if cavity is not None:e[cavity]=1.;l[cavity]=0.
    return np.r_[e,l,x,y,rx,ry,contrast,loss,angle],dict(case=case,seed=CASES[case],shape_stratum='circle' if case in ('dev','blind_normal') else 'ellipse',cavity_cell=cavity)
def maps(z):
    x=np.linspace(-1.6,1.6,401);y=np.linspace(.4,2.8,301);xx,yy=np.meshgrid(x,y)
    dx=xx-z[12];dy=yy-z[13];c=np.cos(z[18]);s=np.sin(z[18]);u=c*dx+s*dy;v=-s*dx+c*dy
    target=(u/z[14])**2+(v/z[15])**2<=1
    ix=np.searchsorted(np.array([-1.5,-.5,.5,1.5]),xx,side='right')-1
    iy=np.searchsorted(np.array([.5,.875,1.25]),yy,side='right')-1
    wall=(xx>=-1.5)&(xx<=1.5)&(yy>=.5)&(yy<=1.25)
    ind=np.clip(iy,0,1)*3+np.clip(ix,0,2)
    we=np.where(wall,z[:6][ind],1.);wl=np.where(wall,z[6:12][ind],0.)
    te=np.where(target,1+z[16],1.);tl=np.where(target,z[17],0.)
    return dict(pixel_x=x,pixel_y=y,target_material=te+1j*tl,target_mask=target,wall_eps=we,wall_loss=wl)
def physics(refine,*,h=None,targetstep=None):
    from r3.physics import JointModel,default_config
    from r3.accounting import BudgetBook
    cfg=default_config();cfg.update(wall_refine=refine)
    if h is not None:cfg['h']=h
    if targetstep is not None:cfg['targetstep']=targetstep
    return JointModel(cfg,BudgetBook(root=ROOT)),cfg

def generate(case,repeats,refine):
    freeze=ROOT/'configs/PROTOCOL.json'
    if not freeze.exists():
        # Root may use this conventional name instead; no implicit authorization.
        freeze=ROOT/'CONTRACT.json'
    if not freeze.exists():raise RuntimeError('Root must write frozen configs/PROTOCOL.json or CONTRACT.json before generation')
    private=ROOT/'inputs/private'/case;private.mkdir(parents=True,exist_ok=True)
    z,meta=sample(case)
    if (private/'truth.npz').exists():
        saved=np.load(private/'truth.npz');clean=saved['clean'];cfg=json.loads((private/'truth.json').read_text())['truth_config']
        if not np.array_equal(saved['z'],z):raise RuntimeError('frozen generator changed')
    else:
        model,cfg=physics(refine);p=model.forward(z,jacobian=False);clean=np.asarray(p.data)
        np.savez_compressed(private/'truth.npz',z=z,clean=clean,**maps(z))
        np.savez_compressed(private/'oracle_wall.npz',eps=z[:6],loss=z[6:12])
        write_json(private/'truth.json',dict(**meta,parameter_names=NAMES,z=z,truth_config=cfg,generated_unix=time.time(),online_exposure='dev exposed after online seal; blind cases private until own seal',truth_mesh_refinement=refine))
    parent=list(CASES).index(case)
    for repeat in range(repeats):
        destination=ROOT/'inputs/online'/f'{case}_r{repeat}';destination.mkdir(parents=True,exist_ok=True)
        if (destination/'observation.npz').exists():raise FileExistsError(destination/'observation.npz')
        seed=33103001+parent*10+repeat;r=np.random.default_rng(seed)
        y=clean+.01*(r.normal(size=clean.shape)+1j*r.normal(size=clean.shape))
        np.savez_compressed(destination/'observation.npz',y=y)
        # Public config is deliberately inverse resolution. No sampled coefficient can enter.
        from r3.physics import default_config
        public=default_config()
        write_json(destination/'config.json',public)
        write_json(destination/'observation.json',dict(case=case,repeat=repeat,noise_seed=seed,sigma_real_imag=.01,data_shape=list(y.shape),packing='frequency,source,receiver; whitened all real then all imaginary',truth_access='OFFLINE_ONLY',measurement_file='observation.npz',config_file='config.json'))
    print(json.dumps(dict(case=case,repeats=repeats,status='GENERATED',truth_private=True,online_keys=['y'],refinement=refine)))

def validate(output,refine,case=None):
    output=Path(output)
    if not (output/'ONLINE_SEALED.json').exists():raise RuntimeError('endpoint must be sealed')
    if (output/'FINER_VALIDATION.npz').exists():raise FileExistsError('validation already exists')
    a=np.load(output/'endpoint.npz');model,cfg=physics(refine,h=.14,targetstep=.075);p=model.forward(a['z'],jacobian=False,label='Independent finer endpoint validation')
    d=np.asarray(p.data)-a['data'];np.savez_compressed(output/'FINER_VALIDATION.npz',data=p.data,difference=d)
    report=dict(status='DONE',refinement=refine,config=cfg,raw_complex_difference_norm=float(np.linalg.norm(d)),raw_complex_RMS_per_channel=float(np.linalg.norm(d)/np.sqrt(d.size)),white_difference_norm=float(np.linalg.norm(d)/.01),white_difference_RMS=float(np.linalg.norm(d)/(.01*np.sqrt(d.size*2))),relative_difference=float(np.linalg.norm(d)/max(np.linalg.norm(p.data),1e-300)),scope='independent finer discrete forward check; no continuum certification')
    if case is not None:
        if case.startswith('blind') and not any(p.exists() for p in [ROOT/'METHODS_FROZEN.json',ROOT/'configs/METHODS_FROZEN.json']):raise RuntimeError('Blind validation truth closed until METHODS_FROZEN.json')
        truth=np.load(ROOT/'inputs/private'/case/'truth.npz');pt=model.forward(truth['z'],jacobian=False,label='Independent finer truth validation');dt=pt.data-truth['clean'];de=p.data-pt.data
        np.savez_compressed(output/'FINER_TRUTH_VALIDATION.npz',truth_data=pt.data,truth_grid_difference=dt,endpoint_prediction_difference=de)
        obs=np.load(ROOT/'inputs/online'/f'{case}_r0'/'observation.npz')['y'];noise=obs-truth['clean']
        from r3.inverse import objective
        original_summary=json.loads((output/'SUMMARY.json').read_text())
        fine_objective=objective(model,obs,a['z'],p)
        report.update(endpoint_original_full_objective=original_summary['objective'],endpoint_fine_full_objective=float(fine_objective),full_objective_change=float(fine_objective-original_summary['objective']),endpoint_original_data_chi2=float(np.sum(abs(a['data']-obs)**2)/.01**2),endpoint_fine_data_chi2=float(np.sum(abs(p.data-obs)**2)/.01**2),sigma_real_imag=.01,n_complex_channels=int(obs.size),n_real_coordinates=int(2*obs.size),expected_noise_white_norm=float(np.sqrt(2*obs.size)),actual_noise_complex_norm=float(np.linalg.norm(noise)),actual_noise_white_norm=float(np.linalg.norm(noise)/.01),actual_noise_white_RMS=float(np.linalg.norm(noise)/(.01*np.sqrt(2*obs.size))),full_received_total_signal_complex_norm=float(np.linalg.norm(pt.data)),full_signal_to_actual_noise_norm=float(np.linalg.norm(pt.data)/np.linalg.norm(noise)),signal_scope='full received total response (fixed incident plus scattering); not target-only SNR',truth_generator_raw_complex_difference_norm=float(np.linalg.norm(dt)),truth_generator_raw_complex_RMS_per_channel=float(np.linalg.norm(dt)/np.sqrt(dt.size)),fine_endpoint_diagnostics=p.diagnostics,fine_truth_diagnostics=pt.diagnostics,no_new_numeric_acceptance_gate=True,no_automatic_PASS=True)
        report.update(case=case,truth_generator_white_difference_norm=float(np.linalg.norm(dt)/.01),truth_generator_white_difference_RMS=float(np.linalg.norm(dt)/(.01*np.sqrt(dt.size*2))),finer_endpoint_vs_truth_white_norm=float(np.linalg.norm(de)/.01),finer_endpoint_vs_truth_white_RMS=float(np.linalg.norm(de)/(.01*np.sqrt(de.size*2))))
    write_json(output/'FINER_VALIDATION.json',report)

def main():
    p=argparse.ArgumentParser();p.add_argument('--case',choices=CASES);p.add_argument('--repeats',type=int,default=1);p.add_argument('--refine',type=int,default=2);p.add_argument('--validate-endpoint',type=Path);p.add_argument('--validation-case',choices=CASES)
    a=p.parse_args()
    if a.validate_endpoint:validate(a.validate_endpoint,a.refine,a.validation_case)
    elif a.case:generate(a.case,a.repeats,a.refine)
    else:p.error('--case or --validate-endpoint is required')
if __name__=='__main__':main()
