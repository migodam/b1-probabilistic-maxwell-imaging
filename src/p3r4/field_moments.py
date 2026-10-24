"""Raw512 measurement-only material mean and pointwise joint moments.

These are risk/moment actions, never an ambient Gaussian posterior density or
joint spatial sampler. Every inference path consumes observation data only.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import torch
from torch import nn

FORMAT = 'p3r4_field_moments_v1'
REFERENCE_FORMAT = 'p3r4_field_moments_reference_v1'
KINDS = ('train_mean', 'mlp_small', 'mlp_medium', 'cnn_small', 'cnn_medium',
         'mlp_small_moments', 'mlp_medium_moments', 'cnn_small_moments', 'cnn_medium_moments',
         'r2_reference')


def packed_observation(y):
    y = np.asarray(y)
    if y.ndim == 3 and y.shape[1:] == (32, 8) and np.iscomplexobj(y):
        y = np.concatenate((y.real.transpose(0,2,1), y.imag.transpose(0,2,1)), axis=2).reshape(len(y),512)
    if y.ndim != 2 or y.shape[1] != 512 or np.iscomplexobj(y) or not np.isfinite(y).all():
        raise ValueError('finite real source-major y[B,512] or complex observation[B,32,8] required')
    return y.astype(np.float32,copy=False)


def unpack_observation(y):
    y = packed_observation(y).reshape(-1,8,64)
    return (y[:,:,:32]+1j*y[:,:,32:]).transpose(0,2,1)


def moment_labels(material, roi, second_moments=False):
    material = np.asarray(material)
    if material.ndim != 2 or material.shape[1] != 576 or not np.iscomplexobj(material) or not np.isfinite(material).all():
        raise ValueError('finite full576 complex TRAIN/validation materials required')
    c = material[:,roi]
    r, i = c.real.astype(np.float64), c.imag.astype(np.float64)
    fields = (r,i,r*r,i*i,r*i) if second_moments else (r,i)
    return np.stack(fields,axis=1)


def material_from_mean(mean_components, roi, host):
    means = np.asarray(mean_components)
    result = np.full((len(means),576), complex(*host), dtype=np.complex128)
    result[:,roi] = means[:,0]+1j*means[:,1]
    return result


def pointwise_covariance(mean, second):
    """[B,2,ROI],[B,3,ROI] -> symmetric per-pixel real2x2 covariance."""
    cov = np.empty((len(mean),mean.shape[-1],2,2),dtype=np.float64)
    cov[:,:,0,0] = second[:,0]-mean[:,0]**2
    cov[:,:,1,1] = second[:,1]-mean[:,1]**2
    cov[:,:,0,1] = cov[:,:,1,0] = second[:,2]-mean[:,0]*mean[:,1]
    return cov


def project_pointwise_psd(covariance):
    """PSD correction of each2x2 block; no spatial covariance is supplied."""
    if np.iscomplexobj(covariance):
        raise ValueError('real joint Re/Im covariance blocks required')
    covariance = np.asarray(covariance,dtype=np.float64)
    if covariance.shape[-2:] != (2,2) or not np.isfinite(covariance).all():
        raise ValueError('finite per-pixel2x2 covariance required')
    eigenvalue,eigenvector = np.linalg.eigh(covariance)
    projected = (eigenvector*np.maximum(eigenvalue,0)[...,None,:])@eigenvector.swapaxes(-1,-2)
    change = np.linalg.norm(projected-covariance,axis=(-2,-1))
    return projected, dict(minimum_raw_eigenvalue=float(eigenvalue.min()),
                           negative_eigenvalue_count=int((eigenvalue<0).sum()),
                           altered_pixel_count=int((change>1e-12).sum()),
                           maximum_frobenius_correction=float(change.max()))


class FieldMomentRegressor(nn.Module):
    """Two64/128-wide MLP layers, or measurement-CNN with raw512 skip.

    Both receive every measurement, with no anchor/context/geometry features.
    The dense ROI head jointly predicts all spatial moment actions.
    """
    def __init__(self,kind='mlp_small',roi_count=256,second_moments=False):
        super().__init__()
        kind=kind.replace('_moments','')
        if kind not in ('mlp_small','mlp_medium','cnn_small','cnn_medium'):
            raise ValueError('MLP/CNN small/medium architecture required')
        self.kind=kind
        self.roi_count=roi_count
        self.channels=5 if second_moments else 2
        width=64 if kind.endswith('small') else 128
        if kind.startswith('mlp'):
            self.encoder=nn.Sequential(nn.Linear(512,width),nn.SiLU(),nn.Linear(width,width),nn.SiLU())
        else:
            c=8 if width==64 else 16
            self.measurement=nn.Sequential(nn.Conv2d(2,c,3,padding=1),nn.SiLU(),
                nn.Conv2d(c,2*c,3,stride=(2,1),padding=1),nn.SiLU(),
                nn.Conv2d(2*c,2*c,3,stride=2,padding=1),nn.SiLU(),
                nn.AdaptiveAvgPool2d((4,2)),nn.Flatten(),nn.Linear(16*c,width),nn.SiLU())
            self.raw_skip=nn.Sequential(nn.Linear(512,32),nn.SiLU())
            self.merge=nn.Sequential(nn.Linear(width+32,width),nn.SiLU())
        self.head=nn.Linear(width,self.channels*roi_count)
        nn.init.zeros_(self.head.weight)
        nn.init.zeros_(self.head.bias)

    def forward(self,y):
        if y.ndim!=2 or y.shape[1]!=512:
            raise ValueError('normalized raw512 input required')
        if self.kind.startswith('mlp'):
            h=self.encoder(y)
        else:
            image=y.reshape(-1,8,2,32).permute(0,2,3,1)
            h=self.merge(torch.cat((self.measurement(image),self.raw_skip(y)),dim=1))
        return self.head(h).reshape(-1,self.channels,self.roi_count)


class FieldMomentPredictor:
    def __init__(self,payload,device='cpu'):
        if payload.get('format')!=FORMAT:
            raise ValueError('measurement-only field-moment checkpoint required')
        self.payload=self.p=payload
        self.device=torch.device(device)
        self.roi=np.asarray(payload['geometry']['roi'],dtype=int)
        self.geometry=payload['geometry']
        self.model=None
        if payload['architecture']!='train_mean':
            self.model=FieldMomentRegressor(payload['architecture'],len(self.roi),payload['second_moments']).to(self.device)
            self.model.load_state_dict(payload['state_dict']);self.model.eval()

    def components(self,y,batch_size=128):
        y=packed_observation(y)
        if batch_size<1:raise ValueError('positive batch size required')
        p=self.payload;parts=[]
        for start in range(0,len(y),batch_size):
            current=y[start:start+batch_size]
            if self.model is None:
                result=np.broadcast_to(p['target_mean'],(len(current),)+p['target_mean'].shape).copy()
            else:
                z=((current-p['y_mean'])/p['y_scale']).astype(np.float32)
                with torch.inference_mode():
                    residual=self.model(torch.from_numpy(z).to(self.device)).cpu().numpy()
                result=residual*p['target_scale']+p['target_mean']
            if not np.isfinite(result).all():raise FloatingPointError('nonfinite field-moment prediction')
            parts.append(result)
        channels=5 if p['second_moments'] else 2
        return np.concatenate(parts).astype(np.float64) if parts else np.empty((0,channels,len(self.roi)))

    def predict(self,y,batch_size=128):
        return material_from_mean(self.components(y,batch_size)[:,:2],self.roi,self.geometry['config']['host_chi'])

    predict_full=predict

    def predict_moments(self,y,batch_size=128):
        raw=self.components(y,batch_size)
        mean=material_from_mean(raw[:,:2],self.roi,self.geometry['config']['host_chi'])
        result=dict(material_mean=mean,second_moment_status='NOT_ESTIMATED_MEAN_ONLY')
        if not self.payload['second_moments']:return result
        cov=pointwise_covariance(raw[:,:2],raw[:,2:])
        if len(raw):
            psd,warning=project_pointwise_psd(cov)
        else:
            psd,warning=cov.copy(),{}
        full_cov=np.zeros((len(raw),576,2,2),np.float64)
        full_psd=np.zeros_like(full_cov)
        full_cov[:,self.roi]=cov;full_psd[:,self.roi]=psd
        host=complex(*self.geometry['config']['host_chi'])
        second=np.broadcast_to(np.array([host.real**2,host.imag**2,host.real*host.imag])[None,:,None],(len(raw),3,576)).copy()
        second[:,:,self.roi]=raw[:,2:]
        moment_status=('EXACT_UNCONDITIONAL_PRIOR_MOMENTS_NOT_POSTERIOR' if self.payload.get('moment_source')=='EXACT_ORIGINAL_P1_UNCONDITIONAL_PRIOR' else
                       'EMPIRICAL_TRAIN_POINTWISE_MOMENTS_NOT_POSTERIOR' if self.payload['architecture']=='train_mean' else
                       'LEARNED_POINTWISE_JOINT_RE_IM_MOMENTS_ONLY')
        result.update(second_moment_status=moment_status,
            raw_second_moments=second,raw_pointwise_covariance=full_cov,psd_pointwise_covariance=full_psd,
            raw_complex_variance=np.trace(full_cov,axis1=-2,axis2=-1),
            psd_complex_variance=np.trace(full_psd,axis1=-2,axis2=-1),psd_projection_warning=warning,
            raw_total_material_variance=np.trace(full_cov,axis1=-2,axis2=-1).mean(1),
            psd_total_material_variance=np.trace(full_psd,axis1=-2,axis2=-1).mean(1),
            raw_material_energy_second_moment=(second[:,0]+second[:,1]).mean(1),
            interpretation='pointwise2x2 covariance actions; no cross-pixel covariance, full posterior, density, or joint sampler')
        return result


class HistoricalR2Reference:
    """Exact historical measurement-only MLP, with a packed-y adapter."""
    def __init__(self,payload,device='cpu',identity=None):
        if payload.get('kind') not in ('mlp_small','mlp_medium'):
            raise ValueError('historical raw-measurement R2 MLP checkpoint required')
        try:
            from p3r2.risk_models import RiskPredictor
        except ModuleNotFoundError:
            from src.p3r2.risk_models import RiskPredictor
        self.reference=RiskPredictor(payload,device)
        self.payload=self.p=payload
        self.identity=identity
        self.model=self.reference.model

    def predict(self,y,batch_size=128):
        return self.reference.predict_full(unpack_observation(y),batch_size)

    predict_full=predict

    def predict_moments(self,y,batch_size=128):
        return dict(material_mean=self.predict(y,batch_size),second_moment_status='NOT_ESTIMATED_HISTORICAL_MEAN_ONLY')


def load_field_moments(path,device='cpu'):
    payload=torch.load(Path(path),map_location='cpu',weights_only=False)
    if payload.get('format')==FORMAT:return FieldMomentPredictor(payload,device)
    if payload.get('format')==REFERENCE_FORMAT:
        result=HistoricalR2Reference(payload['reference_payload'],device,payload.get('reference_identity'))
        result.geometry=payload['geometry']
        return result
    if payload.get('format')=='p3r2_material_risk_v1':
        return HistoricalR2Reference(payload,device)
    raise ValueError('unsupported field-moment or historical R2 measurement-only checkpoint')
