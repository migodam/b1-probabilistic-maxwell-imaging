"""Proper finite correlated-field prior for passive effective TM materials.

No topology labels, peak normalization, neural networks or solver imports.
The probability reference is dz_64 du_2, not material Lebesgue volume.
"""
import numpy as np
from scipy.ndimage import label, gaussian_filter


def saturation(s):
    """C2 quintic with exact host/pure-phase plateaus; not an air decoder."""
    t=np.clip((np.asarray(s)+1)/2,0,1)
    return t*t*t*(10+t*(-15+6*t))


def saturation_prime(s):
    t=np.clip((np.asarray(s)+1)/2,0,1)
    return 15*t*t*(1-t)*(1-t)


class MaterialPrior:
    dimension=66
    gaussian_dimension=64

    def __init__(self,config,points):
        self.config=config;self.points=np.asarray(points,float)
        if self.points.shape!=(576,2):raise ValueError('frozen 24x24 mesh required')
        self.host=complex(*config['host_chi']);self.mask=np.linalg.norm(self.points,axis=1)<=config['roi_radius']
        # Discrete DCT modes have unit full-grid material RMS before weighting.
        xx=self.points[:,0]+.5;yy=self.points[:,1]+.5
        modes=[(k,l) for k in range(8) for l in range(8) if k or l]
        ell=config['correlation_rolloff_length'];weights=np.array([(1+(np.pi*ell)**2*(k*k+l*l))**-1 for k,l in modes])
        columns=[(1 if k==0 else np.sqrt(2))*(1 if l==0 else np.sqrt(2))*np.cos(np.pi*k*xx)*np.cos(np.pi*l*yy) for k,l in modes]
        B=np.column_stack(columns)*weights[None,:]/np.linalg.norm(weights)
        self.A=np.c_[np.full(576,config['offset_SD']),config['field_SD']*B]
        self.temperature=float(config['phase_transition_halfwidth'])
        self.real_interval=np.array(config['inclusion_increment_real_interval']);self.imag_interval=np.array(config['inclusion_increment_imag_interval'])

    def legal(self,q):
        q=np.asarray(q)
        return bool(q.shape==(66,) and np.all(np.isfinite(q)) and np.all(q[64:]>=0) and np.all(q[64:]<=1))

    def draw(self,rng,n=1):
        return np.c_[rng.normal(size=(n,64)),rng.uniform(size=(n,2))]

    def logprior(self,q):
        if not self.legal(q):return -np.inf
        return float(-.5*np.dot(q[:64],q[:64])-32*np.log(2*np.pi))

    def decode(self,q):
        q=np.asarray(q,float)
        if not self.legal(q):raise ValueError('outside declared Gaussian x uniform latent prior')
        f=self.config['field_mean']+self.A@q[:64]
        fraction=saturation(f/self.temperature)*self.mask
        delta=complex(self.real_interval[0]+np.diff(self.real_interval)[0]*q[64],self.imag_interval[0]+np.diff(self.imag_interval)[0]*q[65])
        chi=self.host+delta*fraction
        return chi,dict(fraction=fraction,field=f,inclusion_increment=delta)

    def derivative(self,q):
        _,a=self.decode(q);df=saturation_prime(a['field']/self.temperature)/self.temperature*self.mask
        D=np.empty((576,66),complex);D[:,:64]=a['inclusion_increment']*df[:,None]*self.A
        D[:,64]=np.diff(self.real_interval)[0]*a['fraction'];D[:,65]=1j*np.diff(self.imag_interval)[0]*a['fraction']
        return D

    def material_rms(self,a,b):return float(np.linalg.norm(np.asarray(a)-np.asarray(b))/24)

    def local_prior_reversible(self,q,rng,beta=.15,property_sd=.04):
        """Product pCN and fixed diagonal reflected Gaussian; acceptance uses likelihood only."""
        if not self.legal(q) or not 0<beta<=1:raise ValueError('invalid starting point or beta')
        new=np.asarray(q).copy();new[:64]=np.sqrt(1-beta*beta)*q[:64]+beta*rng.normal(size=64)
        t=np.mod(q[64:]+property_sd*rng.normal(size=2),2);new[64:]=np.where(t<=1,t,2-t)
        return new

    def independent_prior_proposal(self,rng):return self.draw(rng,1)[0]


def topology(fraction,threshold=.5,min_cells=4):
    a=np.asarray(fraction).reshape(24,24)>=threshold
    fg,n=label(a,structure=np.array([[0,1,0],[1,1,1],[0,1,0]]));components=sum(np.sum(fg==k)>=min_cells for k in range(1,n+1))
    inv,m=label(~a,structure=np.ones((3,3)));edges=set(np.r_[inv[0],inv[-1],inv[:,0],inv[:,-1]])
    holes=sum(k not in edges and np.sum(inv==k)>=min_cells for k in range(1,m+1))
    if components==0:event='empty_at_fraction_threshold'
    elif components>=2:event='multiple_components'
    elif holes==0:event='one_component_no_closed_host_pocket'
    elif holes==1:event='one_component_one_closed_host_pocket'
    else:event='one_component_multiple_closed_host_pockets'
    return dict(components=int(components),closed_host_pockets=int(holes),event=event,threshold=threshold,min_cells=min_cells)


def structural_change(prior,q0,q1):
    x0,a=prior.decode(q0);x1,b=prior.decode(q1);df=b['fraction']-a['fraction'];changes=abs(df)>=.15
    cc,n=label(changes.reshape(24,24),structure=np.array([[0,1,0],[1,1,1],[0,1,0]]))
    largest=max([int(np.sum(cc==i)) for i in range(1,n+1)]+[0]);smooth=gaussian_filter(df.reshape(24,24),1)
    return dict(material_RMS=prior.material_rms(x0,x1),fraction_RMS=float(np.linalg.norm(df)/24),smoothed_fraction_RMS=float(np.linalg.norm(smooth)/24),largest_connected_changed_region=largest,task_structure_pass=bool(largest>=4 and np.linalg.norm(smooth)/24>=.01),topology_before=topology(a['fraction']),topology_after=topology(b['fraction']))

