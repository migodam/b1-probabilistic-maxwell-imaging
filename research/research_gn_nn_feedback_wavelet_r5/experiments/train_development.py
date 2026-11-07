"""Finite R5 TRAIN/VAL pilot; no old TEST access or new-blind claim."""
from pathlib import Path
import sys,argparse,json,time,random
import numpy as np
import torch
from torch.utils.data import DataLoader
from threadpoolctl import threadpool_limits

STUDY=Path(__file__).resolve().parents[1]
REPO=STUDY.parents[1]
sys.path.insert(0,str(REPO))
from research.research_gn_nn_feedback_wavelet_r5.src.data import cohort,DevelopmentAnchors,move,run_model,R4
from research.research_gn_nn_feedback_wavelet_r5.src.closed_loop import FeedbackCNN,FixedOperators
from research.research_gn_nn_feedback_wavelet_r5.src.coupled_loss import CoupledTrainingLoss

ARMS={'pixel_closed_coupled':('pixel',True,100.),
      'haar_closed_coupled':('haar',True,100.),
      'db2_closed_coupled':('db2',True,100.),
      'haar_closed_no_fullfield':('haar',True,0.),
      'haar_free_coupled':('haar',False,100.)}

def save(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')

def epoch(model,ops,criterion,loader,stats,device,rng,optimizer=None):
    training=optimizer is not None;model.train(training)
    sums={};cases=0;started=time.perf_counter();grad_norm=None
    with torch.set_grad_enabled(training):
        for cpu_batch in loader:
            batch=move(cpu_batch,device)
            steps=rng.choice((2,4,6)) if training else 4
            if training:optimizer.zero_grad(set_to_none=True)
            out=run_model(model,ops,batch,stats,steps)
            loss,metrics=criterion(out,batch['m_gt'],batch['public_y_pair'],batch['public_noise'],batch['j0'],evaluate_field=not training)
            if training:
                loss.backward()
                grad_norm=float(torch.nn.utils.clip_grad_norm_(model.parameters(),10.).detach())
                optimizer.step()
            n=len(batch['m0']);cases+=n
            err=(out['m_proposal']-batch['m_gt']).square().mean(1).sqrt()
            metrics['material_NRMSE']=(err/batch['amplitude']).mean().detach()
            metrics['material_RMSE']=err.mean().detach()
            metrics['accept_fraction']=out['trace']['accepted'].float().mean().detach()
            metrics['max_energy_increase']=out['trace']['energy_change'].max().detach()
            # All checkpoint selection uses the SAME composite validation score,
            # including actual F, also for the no-full-field training ablation.
            metrics['selection_score']=metrics['material_loss']+100*metrics['full_data_hinge']+.01*metrics['constitutive_normalized']+.01*metrics['fixed_point_tail']
            for key,value in metrics.items():
                if torch.is_tensor(value)and value.numel()==1:
                    sums[key]=sums.get(key,0.)+float(value)*n
            if 'final_field_qualified'in metrics:
                sums['physical_qualified']=sums.get('physical_qualified',0.)+float(metrics['final_field_qualified'].sum())
    return {**{k:v/cases for k,v in sums.items()},'cases':cases,'seconds':time.perf_counter()-started,'last_gradient_norm_before_clipping':grad_norm}

def main():
    p=argparse.ArgumentParser();p.add_argument('--arms',nargs='+',choices=list(ARMS),default=list(ARMS));p.add_argument('--epochs',type=int,default=20);p.add_argument('--batch',type=int,default=8);p.add_argument('--device',default='mps');p.add_argument('--seed',type=int,default=7);p.add_argument('--patience',type=int,default=8);p.add_argument('--wiring-only',action='store_true');p.add_argument('--anchor-mode',choices=['strict','noiseaware'],default='strict');p.add_argument('--noise-ball',action='store_true');a=p.parse_args()
    torch.set_num_threads(1);device=torch.device(a.device)
    anchor_root=STUDY/'inputs/noiseaware_anchors' if a.anchor_mode=='noiseaware' else R4/'anchors'
    train,val,stats,gnorm,snorm=cohort(anchor_root)
    train_set=DevelopmentAnchors(train,stats,gnorm,snorm,anchor_root);val_set=DevelopmentAnchors(val,stats,gnorm,snorm,anchor_root)
    cohort_file='COHORT_NOISEAWARE.json' if a.anchor_mode=='noiseaware' else 'COHORT.json'
    save(STUDY/'audit'/cohort_file,{'scope':'DEVELOPMENT TRAIN/VAL, no TEST opened','anchor_mode':a.anchor_mode,'anchor_root':str(anchor_root),'train_mothers':len({r['parent_id']for r in train}),'train_observations':len(train),'val_mothers':len({r['parent_id']for r in val}),'val_observations':len(val),'train_parent_ids':sorted({r['parent_id']for r in train}),'val_parent_ids':sorted({r['parent_id']for r in val}),'stats':stats,'Gnorm_cached':gnorm,'Snorm_cached':snorm,'step_size':'public-anchor spectral upper bound; actual per-step Armijo','cache_seconds':[train_set.cache_receipt['load_seconds'],val_set.cache_receipt['load_seconds']]})
    ops=FixedOperators(R4/'data/operators.npz').to(device)
    for arm in a.arms:
        ops.call_counts={key:0 for key in ops.call_counts}
        representation,closed,weight=ARMS[arm]
        random.seed(a.seed);np.random.seed(a.seed);torch.manual_seed(a.seed)
        rng=random.Random(a.seed)
        model=FeedbackCNN(representation,closed,noise_ball=a.noise_ball).to(device)
        criterion=CoupledTrainingLoss(R4/'data/operators.npz',ops=ops,lambda_field=weight)
        optimizer=torch.optim.AdamW(model.parameters(),lr=.001,weight_decay=.0001)
        gen=torch.Generator().manual_seed(a.seed)
        train_loader=DataLoader(train_set,batch_size=a.batch,shuffle=True,generator=gen)
        val_loader=DataLoader(val_set,batch_size=a.batch,shuffle=False)
        models_root='models_noiseaware' if a.anchor_mode=='noiseaware' else 'models'
        run_name=arm.replace('_closed_coupled','_closed_noiseball_coupled') if a.noise_ball else arm
        outdir=STUDY/models_root/run_name;outdir.mkdir(parents=True,exist_ok=True)
        config={**vars(a),'arm':run_name,'representation':representation,'closed':closed,'lambda_field':weight,'params':model.parameter_counts(),'TRAIN':128,'VAL':78,'test_status':'NOT_RUN','forward_loss':'full nonlinear Maxwell complex128 retained LU and analytic VJP at every selected inverse stage','inverse_layer_LU':0,'train_depths':[2,4,6],'val_depth':4,'new_weights':'FROM_SCRATCH_NO_R4_INITIALIZATION'}
        save(outdir/'config.json',config)
        history=[];best=float('inf');best_epoch=-1;started=time.perf_counter()
        if a.wiring_only:
            # One actual TRAIN batch, including real Maxwell forward/VJP and
            # learned inverse path. No random algebra/error-floor sweep.
            small=DataLoader(torch.utils.data.Subset(train_set,list(range(min(a.batch,len(train_set))))),batch_size=a.batch)
            tr=epoch(model,ops,criterion,small,stats,device,rng,optimizer)
            save(STUDY/'audit/ACTUAL_TRAIN_WIRING.json',{'arm':arm,'metrics':tr,'physics':criterion.receipt(),'inverse_counts':ops.call_counts,'params':model.parameter_counts()})
            print(json.dumps({'wiring':arm,**tr}),flush=True);break
        for ep in range(1,a.epochs+1):
            tr=epoch(model,ops,criterion,train_loader,stats,device,rng,optimizer)
            va=epoch(model,ops,criterion,val_loader,stats,device,rng)
            row={'epoch':ep,'TRAIN':tr,'VAL':va};history.append(row)
            score=va['selection_score']
            if score<best:
                best,best_epoch=score,ep
                torch.save({'format':'gn_nn_r5_closed_cnn_v1','state_dict':model.state_dict(),'stats':stats,'config':config,'epoch':ep,'VAL_selection_score':score},outdir/'best.pt')
            save(outdir/'history.json',history)
            save(outdir/'physics_cost.json',criterion.receipt())
            print(json.dumps({'arm':arm,'epoch':ep,'TRAIN_NRMSE':tr['material_NRMSE'],'VAL_NRMSE':va['material_NRMSE'],'VAL_data':va.get('relative_public_residual'),'VAL_qualified':va.get('physical_qualified'),'score':score,'seconds':time.perf_counter()-started}),flush=True)
            if ep-best_epoch>=a.patience:break
        save(outdir/'SUMMARY.json',{'status':'TRAINED_DEVELOPMENT_ONLY','epochs':len(history),'best_epoch':best_epoch,'selection_score':best,'wall_seconds':time.perf_counter()-started,'physics_counts':criterion.bridge.counts,'fixed_operator_calls':ops.call_counts,'new_test':'NOT_RUN'})
    print('R5 finite training request complete',flush=True)

if __name__=='__main__':
    with threadpool_limits(limits=1):main()
