"""Independent write-once seal and post-closure scoring; no forward solver.

`seal` never opens sealed_truth. `score --authorize-sealed` requires that seal.
The inverse runner does not import this module. No hashes are computed.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,json,csv,time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
ARMS=['native','global','profile','direct']

def utc():return datetime.now(timezone.utc).isoformat()
def read(p):return json.loads(p.read_text())
def write_new(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x')as f:json.dump(x,f,indent=2,allow_nan=False)

def seal():
    started,cpu=time.perf_counter(),time.process_time();entries=[]
    manifest=read(ROOT/'data_public/manifest.json')
    for record in manifest['records']:
        ident=record['observation_id']
        for arm in ARMS:
            path=ROOT/f'results/blind/{ident}/{arm}/RECEIPT.json';rec=read(path)
            if rec['status']!='SEALED_PREDICTION':raise RuntimeError(f'Unclosed {ident}/{arm}')
            pred=ROOT/rec['prediction']
            with np.load(pred,allow_pickle=False)as d:
                if d['m'].shape!=(576,)or d['Y_prediction'].shape!=(48,12):raise ValueError('Wrong prediction dimensions')
                if not all(np.isfinite(d[k]).all()for k in d.files):raise ValueError('Nonfinite prediction')
            entries.append(dict(observation_id=ident,arm=arm,receipt=str(path.relative_to(ROOT)),prediction=rec['prediction'],
                finished_utc=rec['finished_utc'],status=rec['status'],prediction_bytes=pred.stat().st_size,
                prediction_mtime_ns=pred.stat().st_mtime_ns,receipt_bytes=path.stat().st_size,receipt_mtime_ns=path.stat().st_mtime_ns))
    if len(entries)!=12:raise RuntimeError('Expected3 cases x4 arms')
    result=dict(status='PREDICTIONS_SEALED',sealed_utc=utc(),entries=entries,expected_arms=ARMS,expected_cases=3,
      independent_check='All12 completed receipts and finite576-cell predictions present before any scoring truth access',
      truth_opened=False,code_hash=None,config_hash=None,hash_reason='User forbids hash checks; path/size/mtime and closure timestamps retained',
      source='experiments/score_blind.py::seal',new_physical_LU=0,new_forward_RHS=0,
      wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu,
      timer_scope='Public manifest/12 receipts/prediction validity read; own seal JSON write excluded')
    write_new(ROOT/'audit/PREDICTIONS_SEALED.json',result)
    print(json.dumps(dict(status=result['status'],count=12,sealed_utc=result['sealed_utc'],truth_opened=False)))

def otsu(m):
    if np.ptp(m)<=1e-12:return np.zeros(m.size,bool),float(m[0])
    hist,edges=np.histogram(m,bins=256);centers=(edges[:-1]+edges[1:])/2
    weight=np.cumsum(hist);moment=np.cumsum(hist*centers);den=weight*(weight[-1]-weight)
    between=np.full(256,-np.inf);ok=den>0
    between[ok]=(moment[-1]*weight[ok]-moment[ok]*weight[-1])**2/den[ok]
    threshold=float(centers[np.argmax(between)])
    return m>threshold,threshold

def iou(p,t):
    n=np.count_nonzero(p|t);return float(np.count_nonzero(p&t)/n)if n else 1.

def prior(m):
    z=m.reshape(24,24);d=np.r_[np.diff(z,axis=0).ravel(),np.diff(z,axis=1).ravel()]/19.5
    return float(1e-4*np.mean(d*d)+1e-4*np.mean(np.sqrt(d*d+.01**2)-.01))

def score():
    started,cpu=time.perf_counter(),time.process_time()
    closed=read(ROOT/'audit/PREDICTIONS_SEALED.json')
    if closed['status']!='PREDICTIONS_SEALED'or len(closed['entries'])!=12:raise RuntimeError('Complete pre-truth seal required')
    for e in closed['entries']:
        p=ROOT/e['prediction'];r=ROOT/e['receipt']
        if (p.stat().st_size,p.stat().st_mtime_ns)!=(e['prediction_bytes'],e['prediction_mtime_ns']):raise RuntimeError('Prediction changed after seal')
        if (r.stat().st_size,r.stat().st_mtime_ns)!=(e['receipt_bytes'],e['receipt_mtime_ns']):raise RuntimeError('Receipt changed after seal')
    # FIRST PRIVATE TRUTH ACCESS OCCURS ONLY HERE, after seal verification.
    opened_utc=utc();mapping=read(ROOT/'sealed_truth/mapping_PRIVATE.json')
    meta={x['observation_id']:x for x in mapping['records']}
    rows=[];materials={};truths={};gt_diagnostics=[]
    for ident,info in meta.items():
        with np.load(ROOT/'sealed_truth'/info['truth_file'],allow_pickle=False)as d:
            truth=d['m'].copy();support=d['mask'].ravel().astype(bool);noise=d['noise'].copy()
        with np.load(ROOT/f'data_public/{ident}.npz',allow_pickle=False)as d:
            Y=d['Y'].copy();sigma=float(d['sigma']);radius=float(d['noise_radius']);coords=d['coords'].copy()
        truths[ident]=truth
        selfmask,selfthreshold=otsu(truth)
        gt_diagnostics.append(dict(observation_id=ident,family=info['family'],self_Otsu_IoU=iou(selfmask,support),
          self_Otsu_threshold=selfthreshold,realized_GT_whitened_residual=float(np.linalg.norm(noise)/sigma),
          radius99=radius,realized_GT_in_radius=bool(np.linalg.norm(noise)/sigma<=radius)))
        for arm in ARMS:
            rec=read(ROOT/f'results/blind/{ident}/{arm}/RECEIPT.json')
            with np.load(ROOT/rec['prediction'],allow_pickle=False)as d:m=d['m'].copy();yp=d['Y_prediction'].copy()
            materials[(ident,arm)]=m
            mask,th=otsu(m);nrmse=float(np.linalg.norm(m-truth)/np.sqrt(truth.size)/12.)
            res=float(np.linalg.norm(yp-Y));relative=res/float(np.linalg.norm(Y));white=res/sigma
            physical_phi=.5*relative**2+prior(m)
            lo,hi=.02,20.;violation=float(max(0.,lo-m.min(),m.max()-hi));legal=violation<=1e-12
            datagate=white<=radius;image=nrmse<=.15 and iou(mask,support)>=.7
            centerdistance=float(np.linalg.norm(coords[mask].mean(0)-coords[support].mean(0)))if mask.any()and support.any()else None
            row=dict(observation_id=ident,family=info['family'],arm=arm,stage='POST_CLOSURE_SCORED',amplitude=12.,host=.1,
              amplitude_normalized_nrmse=nrmse,RMSE=12*nrmse,otsu_iou=iou(mask,support),otsu_threshold=th,
              fixed_midpoint_iou=iou(m>6.1,support),foreground_area_error_fraction=float((mask.sum()-support.sum())/support.sum()),
              support_centroid_distance=centerdistance,support_centroid_null_reason='Empty predicted or truth support'if centerdistance is None else None,
              relative_residual=relative,whitened_residual=white,noise_radius99=radius,noise_gate=bool(datagate),
              historical_noisy_relative_gate=bool(relative<=.003),material_gate=bool(nrmse<=.15),image_gate=bool(image),
              material_and_data_gate=bool(nrmse<=.15 and datagate and legal),joint_gate=bool(image and datagate and legal),
              physical_objective=physical_phi,physical_objective_receipt_delta=physical_phi-rec['physical_objective'],
              relative_residual_receipt_delta=relative-rec['physical_relative_residual'],
              whitened_residual_receipt_delta=white-rec['physical_whitened_residual'],material_constraint_violation=violation,
              material_legal=bool(legal),material_min=float(m.min()),material_max=float(m.max()),
              wall_seconds=rec['wall_seconds'],cpu_seconds=rec['cpu_seconds'],physical_cost=rec['physical_model_cost'],
              receipt=f'results/blind/{ident}/{arm}/RECEIPT.json',prediction=rec['prediction'])
            rows.append(row)
    summary=[]
    for arm in ARMS:
        rr=[r for r in rows if r['arm']==arm]
        summary.append(dict(arm=arm,cases=3,noise_pass=sum(r['noise_gate']for r in rr),
          material_pass=sum(r['material_gate']for r in rr),image_pass=sum(r['image_gate']for r in rr),
          joint_pass=sum(r['joint_gate']for r in rr),material_and_data_pass=sum(r['material_and_data_gate']for r in rr),
          median_nrmse=float(np.median([r['amplitude_normalized_nrmse']for r in rr])),median_iou=float(np.median([r['otsu_iou']for r in rr])),
          wall_seconds=sum(r['wall_seconds']for r in rr),cpu_seconds=sum(r['cpu_seconds']for r in rr),
          LU=sum(r['physical_cost'].get('LU',0)for r in rr),actual_complex_RHS=sum(r['physical_cost'].get('actual_complex_RHS',0)for r in rr)))
    result=dict(status='POST_CLOSURE_SCORED',seal='audit/PREDICTIONS_SEALED.json',sealed_utc=closed['sealed_utc'],
      first_truth_opened_utc=opened_utc,source='experiments/score_blind.py',rows=rows,summaries=summary,gt_gate_diagnostics=gt_diagnostics,
      scientific_scope='Three specific fresh objects; no population inference, universal method advantage or independent-root claim',
      selection='Sealed online original-physical-objective incumbent, no GT reselection; no best@K reported',
      gates=dict(NRMSE=.15,OtsuIoU=.7,data='99% marginal per-observation complex Gaussian whitened radius',legality='public box with1e-12 numeric tolerance'),
      new_physical_LU=0,new_forward_RHS=0,code_hash=None,config_hash=None,hash_reason='User forbids hash checks',
      wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu,timer_scope='Seal validation,private/public reads,metrics; output/plot serialization separately below')
    write_new(ROOT/'results/BLIND_SCORES.json',result)
    with (ROOT/'results/BLIND_SCORES.csv').open('x',newline='')as f:
        fields=list(rows[0]);w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for r in rows:w.writerow({k:json.dumps(v)if isinstance(v,(dict,list))else v for k,v in r.items()})
    plot_started,plot_cpu=time.perf_counter(),time.process_time()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    order=sorted(meta,key=lambda i:['ring','curved_bar','two_islands'].index(meta[i]['family']))
    fig,axes=plt.subplots(3,5,figsize=(15,9.5),constrained_layout=True)
    for ri,ident in enumerate(order):
        info=meta[ident]
        for ci,arm in enumerate(['truth']+ARMS):
            ax=axes[ri,ci];z=truths[ident]if arm=='truth'else materials[(ident,arm)]
            im=ax.imshow(z.reshape(24,24),origin='lower',vmin=.02,vmax=20,cmap='viridis',interpolation='nearest')
            ax.set_xticks([]);ax.set_yticks([])
            if ri==0:ax.set_title({'truth':'Held-out truth','native':'Native sourcewise','global':'Global PR','profile':'Known profile','direct':'Direct physical'}[arm],fontsize=12)
            if arm=='truth':ax.set_ylabel(info['family']+'\n'+ident[-8:],fontsize=11);ax.set_xlabel('A=12; host=.1',fontsize=10)
            else:
                row=next(r for r in rows if r['observation_id']==ident and r['arm']==arm)
                gate='JOINT PASS'if row['joint_gate']else('DATA PASS / IMAGE FAIL'if row['noise_gate']else'NOISE GATE FAIL')
                ax.set_xlabel(f"N={row['amplitude_normalized_nrmse']:.3f}; IoU={row['otsu_iou']:.3f}\nr={row['relative_residual']:.4g}; {gate}",fontsize=9)
    fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.88,label='Shared real contrast m; fixed [.02,20]')
    fig.suptitle('Three fresh scalar-TM objects; sealed online selections; no GT-ranked candidates',fontsize=14)
    dest=ROOT/'figures';dest.mkdir(exist_ok=True)
    fig.savefig(dest/'BLIND_RECONSTRUCTIONS.png',dpi=180);fig.savefig(dest/'BLIND_RECONSTRUCTIONS.pdf');plt.close(fig)
    write_new(ROOT/'audit/BLIND_PLOT_RECEIPT.json',dict(status='COMPLETED',source='experiments/score_blind.py',new_physical_LU=0,new_forward_RHS=0,
      wall_seconds=time.perf_counter()-plot_started,cpu_seconds=time.process_time()-plot_cpu,
      input='Stored sealed predictions + post-closure truth',origin='lower; identical raw array order for all panels',
      color_limits=[.02,20.],outputs=['figures/BLIND_RECONSTRUCTIONS.png','figures/BLIND_RECONSTRUCTIONS.pdf'],
      timer_scope='matplotlib imports,layout,two renders; own receipt write excluded'))
    print(json.dumps(dict(status=result['status'],summaries=summary,first_truth_opened_utc=opened_utc)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('operation',choices=['seal','score']);p.add_argument('--authorize-sealed',action='store_true');a=p.parse_args()
    if a.operation=='seal':seal()
    elif a.authorize_sealed:score()
    else:p.error('Score requires parent closure authorization --authorize-sealed')
