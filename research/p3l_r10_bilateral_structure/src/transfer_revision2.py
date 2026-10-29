"""Separately frozen event/domain repair. Never overwrites V1 numerical outcomes."""
from transfer_campaign import *
from scipy.special import ndtr

CONFIG=P/'config/transfer_revision2_frozen_v2.json'

def partition(c):
    mode=np.array(c['mode']);sd=np.array(c['marginal_SD']);cl,ch=np.array(c['core_box']);axes=[]
    for k in range(2):
        points=np.r_[LO[k],HI[k],cl[k],ch[k],mode[k]-2*sd[k],mode[k],mode[k]+2*sd[k],[1.3,.105][k]]
        axes.append(np.unique(np.clip(points,LO[k],HI[k])))
    return axes,cl,ch

def freeze2():
    parents=json.loads((P/'config/transfer_frozen.json').read_text())['conditional_parents'];counts={}
    for parent in parents:
        r=json.loads((OUT/f'{parent}_conditional.json').read_text());axes,cl,ch=partition(r);core=tail=0
        for i in range(len(axes[0])-1):
            for j in range(len(axes[1])-1):
                center=np.array([(axes[0][i]+axes[0][i+1])/2,(axes[1][j]+axes[1][j+1])/2]);inside=np.all(center>=cl)&np.all(center<=ch)
                if inside:core+=1
                else:tail+=1
        # GH6 inner real, GL6 outer imag on every event/peak partition interval.
        inner_nodes=6*6*sum(cl[1]<=(axes[1][j]+axes[1][j+1])/2<=ch[1] for j in range(len(axes[1])-1))
        counts[parent]=dict(core_cells=core,tail_cells=tail,original_nodes=core*(4**2+6**2)+tail*(1**2+2**2),independent_nodes_max=inner_nodes)
    total=sum(v['original_nodes']+v['independent_nodes_max'] for v in counts.values());prior=40608
    c=dict(id='R10_CONDITIONAL_EVENT_SPLIT_REVISION2',frozen_utc=time.time(),approved_by_root=True,parents=parents,
        parent_V1='transfer_frozen.json; V1 status and results retained',reason='repair narrowpeak quadrature and hard event boundary integration; independent final comparisons may remain unqualified',
        source_V1_RHS256=prior,RHS256_cap=70000,LU256_cap=3000,wall_seconds_cap=900,RSS_bytes_cap=6000000000,
        extra_RHS256_cap=29392,root_cap_interpretation='root formally reallocated total70000 at17:52UTC; enforce29392 remaining inclusive all attempts',
        original_rule='event-split fullrectangle composite GL4/6 core and GL1/2 tail; peak boundaries mode+-2SD and+-8SD',
        independent_rule='triangular Gaussian covariance transform; original imaginary outerintervals split at event and peak boundaries, GL6 on LOCAL coreonly; conditional real GH6; legal prior rejection; real event empirical indicator is localcrosscheckonly, main rule eventdomain is exactlysplit; localcheck doesnotcertify full tails',
        counts=counts,planned_nodes_max=total,planned_RHS256_max=8*total,
        tolerances=dict(relative_evidence=2e-4,amplitude_mean=2e-5,event_probability=5e-4),
        fullsupport_certificate=False,GT_access=False)
    if CONFIG.exists():raise RuntimeError('revision2 freeze already exists')
    write(CONFIG,c);print(json.dumps(plain(c),indent=2))


class RevisionBook(Book):
    def close(self,status,error=None):
        row=dict(status=status,error=error,RHS256=self.rhs,LU256=self.lu,RHS576=0,additional_RHS256=self.rhs-self.c['source_V1_RHS256'],wall_seconds=time.perf_counter()-self.start,CPU_seconds=time.process_time()-self.cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        self.event('END_REVISION2',**{k:v for k,v in row.items() if k not in ('RHS256','LU256','RHS576')});write(OUT/'RECEIPT_revision2.json',row)

def load_cached_shape(parent,book):
    shape=Shape.__new__(Shape);shape.parent=parent;shape.book=book
    with np.load(OUT/f'{parent}_setup.npz') as z:
        for k in ('B','L','R','T','Q','LQ','QR','Fh'):setattr(shape,k,z[k].copy())
        shape.y=z['observation'].copy();shape.fraction=z['fraction'].copy()
    shape.I=np.eye(256,dtype=complex)
    return shape

def run_parent(shape,c):
    parent=shape.parent;old=json.loads((OUT/f'{parent}_conditional.json').read_text());axes,cl,ch=partition(old);mode=np.array(old['mode']);H=np.array(old['Hessian']);cov=np.linalg.inv(H);shift=old['profile_phi'];saved=[]
    with np.load(OUT/f'{parent}_setup.npz') as z:fraction=z['fraction'].copy()
    with np.load(R9/'results/stable_roi/SETUP.npz') as z:roi=z['ROI'].copy()
    amean=float(fraction[roi].mean())
    def node(v,w,method,level,region):
        if len(saved)>=c['counts'][parent]['original_nodes']+c['counts'][parent]['independent_nodes_max']:raise RuntimeError('REVISION2_NODE_CAP')
        f,_,_,_=shape.forward(complex(*v));phi=float(.5*np.linalg.norm(white(f-shape.y))**2);saved.append((v[0],v[1],w,phi,method,level,region));return phi
    def moments(nodes,weights,phis,regions,events):
        v=np.asarray(nodes);w=np.asarray(weights);phi=np.asarray(phis);lw=np.log(w)+shift-phi-np.log(.15);z=logsumexp(lw);pn=np.exp(lw-z);mean=pn@v
        return dict(log_evidence=float(z-shift),scaled_evidence=float(np.exp(z)),amplitude_mean=mean,amplitude_variance=pn@((v-mean)**2),event_probability=pn@np.asarray(events),material_ROI_mean=[HOST.real+amean*mean[0],HOST.imag+amean*mean[1]],tail_probability=float(np.sum(pn[np.asarray(regions)!='core'])))
    refs=[]
    for level,(nc,nt) in enumerate([(4,1),(6,2)]):
        nodes=[];weights=[];phis=[];regions=[];events=[]
        for ix in range(len(axes[0])-1):
            for iy in range(len(axes[1])-1):
                aa=np.array([axes[0][ix],axes[1][iy]]);bb=np.array([axes[0][ix+1],axes[1][iy+1]]);center=(aa+bb)/2;core=bool(np.all(center>=cl)&np.all(center<=ch));region='core' if core else 'tail';n=nc if core else nt;x,w=leggauss(n)
                for i in range(n):
                    for j in range(n):
                        v=aa+(bb-aa)*(np.array([x[i],x[j]])+1)/2;ww=w[i]*w[j]*np.prod(bb-aa)/4;nodes.append(v);weights.append(ww);regions.append(region);events.append([center[0]>1.3,center[1]>.105]);phis.append(node(v,ww,'event_split_original_GL',level,region))
        row=moments(nodes,weights,phis,regions,events);row.update(core_order=nc,tail_order=nt,node_count=len(nodes));refs.append(row)
    # Positive triangular covariance transform with imaginary coordinate preserved.
    alpha=cov[0,1]/cov[1,1];conditional_sd=float(np.sqrt(cov[0,0]-cov[0,1]**2/cov[1,1]));zz,ww=hermgauss(6);x,w=leggauss(6);nodes=[];weights=[];phis=[];regions=[];events=[]
    for j in range(len(axes[1])-1):
        aa,bb=axes[1][j:j+2];center=(aa+bb)/2
        if not cl[1]<=center<=ch[1]:continue
        for k in range(6):
            vi=aa+(bb-aa)*(x[k]+1)/2;mu=mode[0]+alpha*(vi-mode[1])
            for i in range(6):
                vr=mu+np.sqrt(2)*conditional_sd*zz[i];v=np.array([vr,vi])
                if np.any(v<LO) or np.any(v>HI):continue
                weight=w[k]*(bb-aa)/2*ww[i]*np.sqrt(2)*conditional_sd*np.exp(zz[i]**2);nodes.append(v);weights.append(weight);regions.append('core');events.append([vr>1.3,center>.105]);phis.append(node(v,weight,'triangular_local_GH6_GL6',2,'local_core'))
    independent=moments(nodes,weights,phis,regions,events);independent.update(method='triangular_local_GH6_GL6',node_count=len(nodes),scope='localimagcore/realconditionalGauss quadrature; no tailcertificate')
    a,b=refs;gh12=old['independent_references'][-1];diff=dict(relative_evidence_refinement=float(abs(np.expm1(a['log_evidence']-b['log_evidence']))),relative_evidence_independent=float(abs(np.expm1(independent['log_evidence']-b['log_evidence']))),relative_evidence_legacy_GH12=float(abs(np.expm1(gh12['log_evidence']-b['log_evidence']))),amplitude_mean_refinement=float(np.max(abs(a['amplitude_mean']-b['amplitude_mean']))),amplitude_mean_independent=float(np.max(abs(independent['amplitude_mean']-b['amplitude_mean']))),event_refinement=float(np.max(abs(a['event_probability']-b['event_probability']))),event_independent=float(np.max(abs(independent['event_probability']-b['event_probability']))))
    tol=c['tolerances'];smooth=all(diff[k]<tol['relative_evidence'] for k in ('relative_evidence_refinement','relative_evidence_independent','relative_evidence_legacy_GH12')) and max(diff['amplitude_mean_refinement'],diff['amplitude_mean_independent'])<tol['amplitude_mean'];event=max(diff['event_refinement'],diff['event_independent'])<tol['event_probability']
    result=dict(parent=parent,status='NUMERICAL_FIXED_SHAPE_TASK_REFERENCE_PASS' if smooth and event else 'REFERENCE_UNQUALIFIED',smooth_evidence_and_mean_pass=bool(smooth),event_numerical_pass=bool(event),strict_certificate=False,fullsupport_certificate=False,V1_status=old['status'],revision_reason='peakpartition+exacthardeventdomain splits; legacyoutputsretained',profile_phi=shift,profile_mode=mode,axis_references=refs,independent_local_reference=independent,differences=diff,prior_density=1/.15,node_count=len(saved),material_mean_rule='host+fraction*E(delta),256ROIaveragingforROItask',GT_access=False)
    write(OUT/f'{parent}_conditional_revision2.json',result)
    with (OUT/f'{parent}_integration_nodes_revision2.csv').open('w') as f:
        writer=csv.writer(f);writer.writerow(['delta_real','delta_imag','integration_weight','phi','method','level','region']);writer.writerows(saved)
    np.savez_compressed(OUT/f'{parent}_conditional_revision2_arrays.npz',amplitude_mean=b['amplitude_mean'],amplitude_variance=b['amplitude_variance'],material_mean=HOST+fraction*(b['amplitude_mean'][0]+1j*b['amplitude_mean'][1]),event_probability=b['event_probability'])
    return result

def run2():
    c=json.loads(CONFIG.read_text())
    if c['planned_RHS256_max']>c['extra_RHS256_cap']:raise RuntimeError('PREFLIGHT_REJECTED_NODES_EXCEED_ALLOCATION')
    if (OUT/'RECEIPT_revision2.json').exists():raise RuntimeError('revision2 no automatic replay')
    book=RevisionBook(c);status='FAILED';error=None;rows=[]
    try:
        book.event('REVISION2_FROZEN_CONFIG',path=CONFIG,planned_RHS256=c['planned_RHS256_max'])
        for parent in c['parents']:
            shape=load_cached_shape(parent,book);rows.append(run_parent(shape,c));write(OUT/'REVISION2_PARTIAL.json',dict(rows=rows,RHS256=book.rhs))
        write(OUT/'RESULTS_revision2.json',dict(status='COMPLETE',parents=rows,RHS256=book.rhs,LU256=book.lu,RHS576=0,additional_RHS256=book.rhs-c['source_V1_RHS256'],qualification='numericalfixedshapeoriginalUniform2tasks; no rigorousfullsupporttailorroundingcertificate'))
        status='COMPLETE';print(json.dumps(plain(dict(status=status,RHS256=book.rhs,parents=[dict(parent=r['parent'],status=r['status'],differences=r['differences']) for r in rows]))),flush=True)
    except BaseException as exc:error=repr(exc);write(OUT/'FAILURE_revision2.json',dict(error=error));raise
    finally:book.close(status,error)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('command',choices=['freeze','run']);ns=a.parse_args();freeze2() if ns.command=='freeze' else run2()
