"""Final phase revision: joint exact-bin moments tilted from raw online anchor.

GT never enters the baseline. Anchor q0 is native(z,u); raw chi0 supplies phase.
The learned output remains a joint weighted-event moment risk action.
"""
import numpy as np
import torch
from .risk_bin import BinRisk
from .risk import feature_image,full_material

FORMAT='p3r4_anchorbin_risk_v2'
EPSILON=1e-4
FRACTION_TOLERANCE=2e-5


def anchor_mass_from_components(anchor,delta,host,bounds,k,epsilon=EPSILON,tolerance=FRACTION_TOLERANCE):
    """anchor[N,2,R],delta[N,2]. Zero-delta components are never divided.

    Fraction is common to Re/Im; both identifiable ratios must agree. For zero
    delta in both components, a host anchor has zero joint moments regardless
    of phase, so choose phase0. Fixed smoothing still adds eps*B/2 at update0.
    """
    if anchor.ndim!=3 or anchor.shape[1]!=2 or delta.shape!=anchor.shape[:2]:raise ValueError('raw component dimensions')
    if not all(torch.isfinite(a).all() for a in (anchor,delta,host,bounds)):raise ValueError('finite raw anchor and properties required')
    if torch.any(bounds<=0) or torch.any(delta<0) or torch.any(delta>bounds[None]+tolerance):raise ValueError('nonnegative legal delta and positive component bounds required')
    contrast=anchor-host[None,:,None];active=delta>1e-12
    ratio=contrast/torch.where(active,delta,torch.ones_like(delta))[:,:,None]
    for c in range(2):
        if torch.any(active[:,c,None]&((ratio[:,c]<-tolerance)|(ratio[:,c]>1+tolerance))):raise ValueError('anchor fraction outside original phase support')
        if torch.any((~active[:,c,None])&(torch.abs(contrast[:,c])>tolerance)):raise ValueError('nonhost material with zero increment')
    identifiable=active.all(1)
    if torch.any(identifiable[:,None]&(torch.abs(ratio[:,0]-ratio[:,1])>tolerance)):raise ValueError('raw anchor Re/Im phases disagree')
    denominator=(delta*delta).sum(1)[:,None]
    phase=(delta[:,:,None]*contrast).sum(1)/torch.where(denominator>0,denominator,torch.ones_like(denominator))
    phase=torch.clamp(phase,0,1)  # Only previously checked floating tolerance.
    ramp=torch.clamp(k*phase[:,None,:]-torch.arange(k,device=anchor.device,dtype=anchor.dtype)[None,:,None],0,1)
    moments=delta[:,:,None,None]*ramp[:,None]
    scaled=moments/bounds[None,:,None,None]
    masses=torch.cat((1-scaled[:,:,:1],scaled[:,:,:-1]-scaled[:,:,1:],scaled[:,:,-1:]),2)
    if torch.any(masses<-tolerance):raise ValueError('negative baseline mass')
    masses=(1-epsilon)*masses+epsilon/(k+1)
    return masses,moments,phase


class AnchorBinRisk(BinRisk):
    def __init__(self,kind,geometry):
        super().__init__(kind,geometry);self.epsilon=EPSILON

    def baseline(self,q0,chi0):
        complex_form=torch.is_complex(chi0) and chi0.shape==(len(q0),576)
        component_form=not torch.is_complex(chi0) and chi0.shape==(len(q0),2,576)
        if q0.ndim!=2 or q0.shape[1]!=66 or not (complex_form or component_form):raise ValueError('explicit native q0[N,66] and raw chi0 complex[N,576] or realified[N,2,576] required')
        if not torch.isfinite(q0).all() or torch.any((q0[:,64:]<0)|(q0[:,64:]>1)):raise ValueError('native online attributes required')
        re=self.cfg['inclusion_increment_real_interval'];im=self.cfg['inclusion_increment_imag_interval']
        delta=torch.stack((re[0]+(re[1]-re[0])*q0[:,64],im[0]+(im[1]-im[0])*q0[:,65]),1)
        anchor=torch.stack((chi0[:,self.roi].real,chi0[:,self.roi].imag),1) if complex_form else chi0[:,:,self.roi]
        return anchor_mass_from_components(anchor,delta,self.host,self.bounds,self.k,self.epsilon)

    def forward(self,image,q0=None,chi0=None,return_moments=False):
        if q0 is None or chi0 is None:raise ValueError('raw online q0/chi0 are mandatory; normalized features cannot substitute')
        base_mass,_,_=self.baseline(q0,chi0)
        correction=self.head(self.local(image)).flatten(2)[:,:,self.roi].reshape(len(image),2,self.k+1,-1)
        masses=torch.softmax(torch.log(base_mass)+correction,dim=2)
        moments=torch.flip(torch.cumsum(torch.flip(masses[:,:,1:],dims=[2]),dim=2),dims=[2])*self.bounds[None,:,None,None]
        mean=self.host[None,:,None]+torch.sum(moments*self.weights[None,None,:,None],dim=2)
        return (mean,moments) if return_moments else mean


class AnchorBinPredictor:
    def __init__(self,path,device='cpu'):
        self.p=torch.load(path,map_location='cpu',weights_only=False)
        if self.p.get('format')!=FORMAT or self.p.get('anchor_epsilon')!=EPSILON:raise ValueError('explicit revision2 format and frozen epsilon required')
        self.geometry=self.p['geometry'];self.kind=self.p['kind'];self.device=device
        self.model=AnchorBinRisk(self.kind,self.geometry).to(device);self.model.load_state_dict(self.p['state_dict']);self.model.eval()

    def predict(self,y,anchor=None,batch_size=128):
        if anchor is None:raise ValueError('explicit raw online anchor required')
        rows=[]
        for start in range(0,len(y),batch_size):
            online={k:anchor[k][start:start+batch_size] for k in ('observation','residual','q0','chi0','J')}
            image=feature_image(y[start:start+batch_size],self.geometry,online,self.p['feature_normalization'])
            if self.p.get('context_mode')=='decoder_only':image[:,2:6]=0
            with torch.inference_mode():
                value=self.model(torch.tensor(image,device=self.device),torch.tensor(online['q0'],dtype=torch.float32,device=self.device),torch.tensor(np.stack((online['chi0'].real,online['chi0'].imag),1),dtype=torch.float32,device=self.device)).cpu().numpy()
            rows.append(full_material(value,self.geometry))
        return np.concatenate(rows)


def load_risk_compatible(path,device='cpu'):
    """Format-aware route: old checkpoints keep their original inference path."""
    header=torch.load(path,map_location='cpu',weights_only=False)
    if header.get('format')==FORMAT:return AnchorBinPredictor(path,device)
    from .risk import RiskActionPredictor
    return RiskActionPredictor(path,device)
