"""P1 material risk actions. No solver calls and no posterior sampling claim.

Threshold outputs regress joint delta*event labels, not independent phase and
property factors. Inference accepts an explicit online tensor only.
"""
from pathlib import Path
import numpy as np
import torch
from torch import nn
from scipy.special import ndtr
from p3r2.residual_risk import load_residual_predictor, online_features, normalize
from p3r2.risk_models import observation_features
from p1_material.prior import MaterialPrior

KINDS=('r2_cnn','local_cnn','no_geometry_cnn','phase8','phase16','direct_moment','affine','sharpen','spatial')


def frozen_geometry(checkpoint):
    return torch.load(checkpoint,map_location='cpu',weights_only=False)['geometry']


def phase_labels(x, geometry, thresholds, s=None):
    """x is Gaussianized66 training truth; only this label builder sees it."""
    cfg=geometry['config'];roi=geometry['roi']
    if s is None:
        s=(cfg['field_mean']+np.asarray(x[:,:64])@np.asarray(geometry['A'])[roi].T)/cfg['phase_transition_halfwidth']
    if s.shape!=(len(x),len(roi)): raise ValueError('s must be standardized ROI field')
    u=ndtr(x[:,64:]);r=cfg['inclusion_increment_real_interval'];i=cfg['inclusion_increment_imag_interval']
    delta=np.stack((r[0]+np.diff(r)[0]*u[:,0],i[0]+np.diff(i)[0]*u[:,1]),1)
    return (delta[:,:,None,None]*(s[:,None,None,:]>np.asarray(thresholds)[None,None,:,None])).astype('float32')


def feature_image(y, geometry, anchor=None, normalization=None):
    """8 online channels: R2 image6 plus inferred u2, or raw-y compressed control.

    No label dictionary is accepted. Anchor is whitelisted inside R2's feature
    builder; normalization must have been fit on TRAIN or frozen R2 TRAIN.
    """
    if anchor is not None:
        allowed={k:anchor[k] for k in ('observation','residual','q0','chi0','J')}
        f=online_features(allowed,geometry)
        if normalization is None: raise ValueError('TRAIN feature normalization required')
        f=normalize(f,normalization)
        return np.concatenate((f['image'],np.broadcast_to(f['base'][:,-2:,None,None],(len(y),2,24,24))),1).astype('float32')
    # Fixed, label-free raw measurement compression. This is a weaker information
    # set than raw512 and must never be described as full-y conditional mean.
    a=np.asarray(y).reshape(len(y),8,64)
    return np.broadcast_to(a.mean(2)[:,:,None,None],(len(y),8,24,24)).copy().astype('float32')


class LocalRisk(nn.Module):
    def __init__(self,kind,geometry):
        super().__init__();self.kind=kind;self.cfg=geometry['config']
        self.register_buffer('roi',torch.as_tensor(geometry['roi'],dtype=torch.long))
        self.k=int(kind[5:]) if kind.startswith('phase') else 0
        self.local=nn.Sequential(nn.Conv2d(8,16,3,padding=1),nn.GELU(),nn.Conv2d(16,16,3,padding=1),nn.GELU(),nn.Conv2d(16,16,3,padding=1),nn.GELU())
        # K+1 nonnegative masses yield bounded monotone joint moment functions.
        self.head=nn.Conv2d(16,2*(self.k+1) if self.k else 2,1)
        nn.init.zeros_(self.head.weight);nn.init.zeros_(self.head.bias)
        if self.k:
            t,w=np.polynomial.legendre.leggauss(self.k);w=w*15*(1-t*t)**2/16
            self.register_buffer('thresholds',torch.tensor(t,dtype=torch.float32))
            self.register_buffer('weights',torch.tensor(w,dtype=torch.float32))
        bounds=[self.cfg['inclusion_increment_real_interval'][1],self.cfg['inclusion_increment_imag_interval'][1]]
        self.register_buffer('bounds',torch.tensor(bounds,dtype=torch.float32))
        self.register_buffer('host',torch.tensor(self.cfg['host_chi'],dtype=torch.float32))

    def forward(self,image,anchor=None,return_moments=False):
        if self.kind=='no_geometry_cnn':
            image=torch.cat((image[:,:2],torch.zeros_like(image[:,2:6]),image[:,6:]),1)
        raw=self.head(self.local(image)).flatten(2)[:,:,self.roi]
        if self.k:
            mass=torch.softmax(raw.reshape(len(image),2,self.k+1,-1),dim=2)
            # Positive masses sum to the maximum legal increment. Tail sums can
            # represent any bounded nonnegative decreasing moment sequence.
            moments=torch.flip(torch.cumsum(torch.flip(mass[:,:,1:],dims=[2]),dim=2),dims=[2])*self.bounds[None,:,None,None]
            mean=self.host[None,:,None]+torch.sum(moments*self.weights[None,None,:,None],dim=2)
            return (mean,moments) if return_moments else mean
        raw=raw.reshape(len(image),2,-1)
        return raw if anchor is None else anchor+raw


def load_r2_baseline(path,device='cpu'):
    """Exact historical5842 parameter predictor; no reimplementation of decoder."""
    return load_residual_predictor(Path(path),device)


def roi_components(material,roi):
    a=np.asarray(material)[:,roi]
    return np.stack((a.real,a.imag),1).astype('float32')


def full_material(components,geometry):
    out=np.full((len(components),576),complex(*geometry['config']['host_chi']),dtype=np.complex64)
    out[:,geometry['roi']]=components[:,0]+1j*components[:,1]
    return out


def classical(image,anchor,kind,coefficient=1.):
    """Frozen scalar calibration/sharpen/spatial controls, trained/selected offline."""
    if kind=='affine': return anchor*coefficient[0]+np.asarray(coefficient[1])[:,None]
    if kind=='sharpen': return np.clip((anchor-coefficient[0])*coefficient[1]+coefficient[0],coefficient[2],coefficient[3])
    if kind=='spatial':
        from scipy.ndimage import uniform_filter
        # Operates on full-grid image; caller gathers ROI afterwards.
        return uniform_filter(image[:,:2],size=(1,1,3,3),mode='nearest')
    raise ValueError(kind)


class RiskActionPredictor:
    """Saved neural risk action with strict online-only context inference."""
    def __init__(self,path,device='cpu'):
        self.p=torch.load(path,map_location='cpu',weights_only=False);self.device=device
        p=self.p;self.geometry=p['geometry'];self.kind=p['kind']
        if self.kind=='r2_cnn':
            from p3r2.residual_risk import ResidualPredictor
            self.baseline=ResidualPredictor(p['historical_payload'],device);self.model=None
        elif self.kind in ('affine','sharpen','spatial'):self.model=None
        else:
            self.model=LocalRisk(self.kind,self.geometry).to(device);self.model.load_state_dict(p['state_dict']);self.model.eval()

    def predict(self,y,anchor=None,batch_size=128):
        if self.kind=='r2_cnn':return self.baseline.predict_full(anchor,batch_size)
        p=self.p;g=self.geometry;rows=[]
        for start in range(0,len(y),batch_size):
            yy=y[start:start+batch_size]
            aa=None if anchor is None else {k:anchor[k][start:start+batch_size] for k in ('observation','residual','q0','chi0','J')}
            image=feature_image(yy,g,aa,p['feature_normalization'])
            if aa is None:
                image=(image-p['feature_mean'])/p['feature_scale'];base=np.zeros((len(yy),2,len(g['roi'])),dtype='float32')
            else:base=roi_components(aa['chi0'],g['roi'])
            if p.get('context_mode')=='decoder_only':image[:,2:6]=0
            if self.model is not None:
                with torch.inference_mode():value=self.model(torch.tensor(image,device=self.device),None if self.model.k else torch.tensor(base,device=self.device)).cpu().numpy()
            elif self.kind=='affine':
                coef=np.asarray(p['coefficient']);value=base*coef[None,:,0,None]+coef[None,:,1,None]
            elif self.kind=='spatial':
                from scipy.ndimage import uniform_filter
                chi=aa['chi0'];im=np.stack((chi.real,chi.imag),1).reshape(-1,2,24,24)
                value=uniform_filter(im,size=(1,1,3,3),mode='nearest').reshape(-1,2,576)[:,:,g['roi']]
            else:
                cfg=g['config'];host=np.array(cfg['host_chi'])[None,:,None];q=aa['q0'];r=cfg['inclusion_increment_real_interval'];i=cfg['inclusion_increment_imag_interval']
                delta=np.stack((r[0]+np.diff(r)[0]*q[:,64],i[0]+np.diff(i)[0]*q[:,65]),1)[:,:,None]
                fraction=np.clip((base-host)/np.maximum(delta,1e-8),0,1)
                value=host+delta*np.clip(.5+p['slope']*(fraction-.5),0,1)
            rows.append(full_material(value,g))
        return np.concatenate(rows)
