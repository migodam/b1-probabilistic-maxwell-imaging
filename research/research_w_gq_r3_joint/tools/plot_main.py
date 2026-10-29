"""Main r0 scientific figures from sealed predictions and post-seal scores only."""
from pathlib import Path
import os,json,argparse
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from generate import ROOT,maps
FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
CASES=['blind_normal','blind_cavity','blind_strong'];METHODS=['FULL_GN','FIXED_CURRENT_CHART','ADAPTIVE_CURRENT_CHART'];COLORS=['#1f77b4','#d97911','#238b45'];LABELS=['Ordinary joint GN','Fixed forcing-current\ntangent chart','Adaptive forcing-current\ntangent chart'];CASE_LABELS=['Normal wall / circle','Boundary-connected air defect\n(historical ID: blind_cavity)','Strong target / ellipse']
plt.rcParams.update({'font.size':10,'axes.titlesize':11,'savefig.dpi':170})
def read(p):return json.loads(p.read_text())
def output(c,m):return ROOT/'results/online'/f'{c}_r0'/m
def load():return {(c,m):read(output(c,m)/'SCORES.json') for c in CASES for m in METHODS if (output(c,m)/'SCORES.json').exists()}
def save(fig,name):fig.savefig(FIG/(name+'.png'),bbox_inches='tight');fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
def draw_image(ax,z,title):
    q=maps(z);im=ax.imshow(q['target_material'].real,origin='lower',extent=[q['pixel_x'][0],q['pixel_x'][-1],q['pixel_y'][0],q['pixel_y'][-1]],aspect='equal',vmin=1,vmax=9,cmap='viridis',interpolation='nearest');ax.set(xlim=(-1,1),ylim=(1.6,2.65),xlabel=r'$x/\lambda_0$',title=title);return im

def material_overview(scores):
    fig,axes=plt.subplots(3,4,figsize=(14,8.4),layout='constrained')
    for i,c in enumerate(CASES):
        truth=np.load(ROOT/'inputs/private'/c/'truth.npz')['z'];im=draw_image(axes[i,0],truth,'Offline truth' if i==0 else '')
        for j,m in enumerate(METHODS,1):
            q=scores[c,m]['selected'];a=np.load(output(c,m)/'endpoint.npz');im=draw_image(axes[i,j],a['z'],(LABELS[j-1]+'\n' if i==0 else '')+f"position {q['position_error']:.3f}; IoU {q['shape_IoU']:.3f}")
        axes[i,0].set_ylabel(CASE_LABELS[i]+'\n'+r'$y/\lambda_0$')
    fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.75,label='Target real permittivity at f=1');fig.suptitle('Three frozen independent blind worlds, noise repeat 0\nFull nonlinear Maxwell F in every arm; forcing-current charts approximate material tangent J\nEndpoints selected by full measurement objective before offline truth scoring',fontsize=13);save(fig,'R0_TARGET_MATERIAL_OVERVIEW')

def wall_overviews(scores):
    for field,limit,label,name in [('wall_eps',(1,6),'Wall real permittivity','R0_WALL_EPS_OVERVIEW'),('wall_loss',(0,.25),'Wall loss coefficient at f=1','R0_WALL_LOSS_OVERVIEW')]:
        fig,axes=plt.subplots(3,4,figsize=(15,7.2),layout='constrained')
        for i,c in enumerate(CASES):
            truth=np.load(ROOT/'inputs/private'/c/'truth.npz')['z'];sets=[truth]+[np.load(output(c,m)/'endpoint.npz')['z'] for m in METHODS]
            for j,z in enumerate(sets):
                q=maps(z);im=axes[i,j].imshow(q[field],origin='lower',extent=[q['pixel_x'][0],q['pixel_x'][-1],q['pixel_y'][0],q['pixel_y'][-1]],aspect='equal',vmin=limit[0],vmax=limit[1],cmap='viridis',interpolation='nearest');axes[i,j].set(xlim=(-1.5,1.5),ylim=(.5,1.25),xlabel=r'$x/\lambda_0$')
                if j==0:title='Offline true wall' if i==0 else 'Offline true wall'
                else:
                    s=scores[c,METHODS[j-1]]['selected'];key='wall_epsilon_RMSE' if field=='wall_eps' else 'wall_loss_RMSE';title=(LABELS[j-1]+'\n' if i==0 else '')+f"wall RMSE {s[key]:.4f}; target pos {s['position_error']:.3f}"
                axes[i,j].set_title(title,fontsize=10)
            axes[i,0].set_ylabel(CASE_LABELS[i]+'\n'+r'$y/\lambda_0$')
        fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.8,label=label);fig.suptitle('Joint unknown-wall / target inversion: three blind r0 parents\nAll six macro cells unknown online; boundary-connected air defect is not an enclosed cavity',fontsize=13);save(fig,name)

def quality_cost(scores):
    agg=read(ROOT/'results/aggregation/AGGREGATE.json');rows={(r['parent'],r['method']):r for r in agg['rows'] if r['noise']==0 and r['status']=='POST_SEAL_SCORED'}
    fig,ax=plt.subplots(2,3,figsize=(14,8),layout='constrained');xs=np.arange(3)
    for j,m in enumerate(METHODS):
        vals=[scores[c,m]['selected'] for c in CASES]
        for a,key,title,gate in [(ax[0,0],'position_error','Position error / reference wavelength',.1),(ax[0,1],'shape_IoU','Shape IoU',.7),(ax[0,2],'epsilon_real_relative_error','Target real permittivity relative error',.2)]:
            a.plot(xs,[v[key] for v in vals],marker=['o','s','^'][j],color=COLORS[j],label=LABELS[j].replace('\n',' '));a.axhline(gate,color='k',ls='--',lw=.7);a.set_title(title);a.set_xticks(xs,['normal','boundary air defect','strong'],rotation=15)
        rr=[rows[c,m] for c in CASES];rhs=sum(r['rhs'] for r in rr);fe=sum(r['full_fe'] for r in rr);ax[1,0].bar(j,rhs,color=COLORS[j]);ax[1,0].text(j,rhs,f'{rhs:,} RHS\n{fe:,} FE',ha='center',va='bottom',fontsize=9)
        wall=sum(r['wall_seconds'] for r in rr);ax[1,1].bar(j,wall,color=COLORS[j]);ax[1,1].text(j,wall,f'{wall:.1f}',ha='center',va='bottom')
        for i,c in enumerate(CASES):
            cpu=rows[c,m]['CPU_seconds'];x=i+(j-1)*.23
            if cpu is None:ax[1,2].plot(x,0,marker='x',color=COLORS[j]);ax[1,2].annotate('NOT RECORDED',(x,0),xytext=(0,10),textcoords='offset points',ha='center',va='bottom',fontsize=7,rotation=90,color=COLORS[j])
            else:ax[1,2].bar(x,cpu,width=.21,color=COLORS[j])
    ax[0,0].legend(fontsize=8,loc='upper left');ax[1,0].set(title='Cold method costs across three matched r0 cases',ylabel='Reserved complex RHS (FE counts annotated)',xticks=range(3),xticklabels=['GN','Fixed tangent','Adaptive tangent']);ax[1,0].margins(y=.2)
    ax[1,1].set(title='Actual end-to-end wall seconds, all four starts',ylabel='seconds',xticks=range(3),xticklabels=['GN','Fixed tangent','Adaptive tangent']);ax[1,1].margins(y=.15)
    ax[1,2].set(title='CPU seconds by parent; gaps are not zero',ylabel='seconds',xticks=xs,xticklabels=['normal','boundary air defect','strong']);fig.suptitle('Quality and paid cost must be read together\nFixed-chart early stalls are not a matched-quality acceleration; two missing CPU receipts remain missing',fontsize=13);save(fig,'R0_QUALITY_COST')

def branches(scores):
    fig,ax=plt.subplots(1,3,figsize=(14,4.6),layout='constrained')
    for i,c in enumerate(CASES):
        for j,m in enumerate(METHODS):
            s=scores[c,m];summary=read(output(c,m)/'SUMMARY.json');recs={r['stage']:r for r in summary['starts']}
            for b in s['branches']:
                k=int(Path(b['endpoint']).name.split('_')[1]);rec=recs.get(f'joint_start_{k}',{});x=rec.get('objective',b.get('objective'));ax[i].scatter(x,b['position_error'],s=45,color=COLORS[j],marker=['o','s','^'][j],alpha=.75,label=LABELS[j].replace('\n',' ') if k==0 else None)
            q=s['selected'];ax[i].scatter(summary['objective'],q['position_error'],s=160,color=COLORS[j],marker='*',edgecolor='black',linewidth=.7)
        ax[i].axhline(.1,color='k',ls='--',lw=.8);ax[i].set(yscale='log',xlabel='Full measurement objective + frozen prior',ylabel=r'Position error / $\lambda_0$ (log scale)',title=CASE_LABELS[i])
    ax[0].legend(fontsize=8,loc='best');fig.suptitle('All four retained starts per method; stars are measurement-selected endpoints\nOffline errors never selected a branch; branch-level raw values are in BRANCH_RESULTS.csv',fontsize=13);save(fig,'R0_ALL_BRANCHES')

def iterations():
    fig,axes=plt.subplots(3,3,figsize=(14,10.5),layout='constrained')
    for i,c in enumerate(CASES):
        for j,m in enumerate(METHODS):
            a=axes[i,j];d=output(c,m);tr=[json.loads(l) for l in (d/'TRACE.jsonl').read_text().splitlines() if l];steps=[json.loads(l) for l in (d/'STEPS.jsonl').read_text().splitlines() if l];accepted={r['candidate'] for r in steps if r.get('event')=='TR_DECISION' and r.get('accepted')};summary=read(d/'SUMMARY.json')
            for k in range(4):
                rr=[r for r in tr if r.get('stage')==f'joint_start_{k}'];x=np.arange(1,len(rr)+1);ys=[r['objective'] for r in rr];a.plot(x,ys,lw=.9,alpha=.75,label=f'start {k}');sel=[v for v,r in enumerate(rr) if r['candidate'] in accepted];a.scatter(x[sel],[ys[v] for v in sel],s=8)
            a.axhline(summary['objective'],color='black',ls='--',lw=.7);a.set(yscale='log',xlabel='Full-objective evaluations within start',ylabel='Full objective (log scale)',title=(LABELS[j]+'\n' if i==0 else '')+CASE_LABELS[i]);a.grid(alpha=.18)
    axes[0,0].legend(fontsize=8);fig.suptitle('Actual joint-start trajectories, including rejected full-F trials\nSmall dots mark accepted steps; wall initialization and 18 coarse seeds are paid in cost totals\nFinal exact-J stationarity and stop reasons are reported separately',fontsize=13);save(fig,'R0_OBJECTIVE_TRAJECTORIES')

def controls():
    labels=['Offline truth','Prior parameter-mean point','Wall-first fixed','Known-wall oracle (diagnostic)'];methods=['PRIOR_ONLY','WALL_FIRST_FIXED','KNOWN_WALL_ORACLE'];fig,axes=plt.subplots(3,4,figsize=(14,8.5),layout='constrained')
    for i,c in enumerate(CASES):
        im=draw_image(axes[i,0],np.load(ROOT/'inputs/private'/c/'truth.npz')['z'],labels[0] if i==0 else '')
        for j,m in enumerate(methods,1):
            d=output(c,m)
            if (d/'SCORES.json').exists():
                q=read(d/'SCORES.json')['selected'];im=draw_image(axes[i,j],np.load(d/'endpoint.npz')['z'],(labels[j]+'\n' if i==0 else '')+f"position {q['position_error']:.3f}; IoU {q['shape_IoU']:.3f}")
            else:axes[i,j].set_facecolor('#eeeeee');axes[i,j].text(.5,.5,'NOT YET SCORED\nNo result imputed',ha='center',va='center',transform=axes[i,j].transAxes);axes[i,j].set_title(labels[j] if i==0 else '')
        axes[i,0].set_ylabel(CASE_LABELS[i]+'\n'+r'$y/\lambda_0$')
    fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.75,label='Target real permittivity at f=1');fig.suptitle('Matched ordinary controls, noise repeat 0\nPrior point is not a Bayes-optimal material image; known-wall arm reads only wall truth after nonoracle seals',fontsize=13);save(fig,'R0_MATCHED_CONTROLS')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--only',default='material,walls,quality,branches,iterations,controls');args=ap.parse_args();todo=args.only.split(',');s=load()
    if len(s)!=9:raise RuntimeError('All nine r0 score files required')
    if 'material' in todo:material_overview(s)
    if 'walls' in todo:wall_overviews(s)
    if 'quality' in todo:quality_cost(s)
    if 'branches' in todo:branches(s)
    if 'iterations' in todo:iterations()
    if 'controls' in todo:controls()
    (FIG/'R0_FIGURE_MANIFEST.json').write_text(json.dumps(dict(figures=['R0_TARGET_MATERIAL_OVERVIEW','R0_WALL_EPS_OVERVIEW','R0_WALL_LOSS_OVERVIEW','R0_QUALITY_COST','R0_ALL_BRANCHES','R0_OBJECTIVE_TRAJECTORIES','R0_MATCHED_CONTROLS'],formats=['png','pdf'],source='sealed endpoint arrays, post-seal scores, raw trace/steps, runtime receipts via aggregation',forward='full nonlinear Maxwell in all methods; charts only approximate material J',n_independent_parents=3,noise_repeat=0,oracle='diagnostic, wall coefficients only',cavity_scope='boundary-connected air defect, historical ID blind_cavity'),indent=2)+'\n')
    print('Five main figure families saved in figures/, PNG+PDF')
if __name__=='__main__':main()
