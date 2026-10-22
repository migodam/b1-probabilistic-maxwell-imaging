"""Exact channel-removal controls for the original5842-parameter local CNN.

Training delegates unchanged optimization/checkpoint selection to the original
train_residual under a temporary feature hook. Checkpoints require this loader;
their explicit new kind cannot be silently loaded as an unmasked original CNN.
"""
from contextlib import contextmanager
from pathlib import Path
import json
import numpy as np
import torch
from p3r2 import residual_risk as original
from p3r2.risk_training import material_targets

FAMILIES=('sensitivity_only','backprojection_only')
FORMAT='p3r2_masked_mechanism_risk_v1'
EXPECTED_PARAMETERS=5842


def mask_contract(family):
    if family not in FAMILIES: raise ValueError('unknown mechanism mask')
    keep=[0,1,2,3] if family=='sensitivity_only' else [0,1,4,5]
    return dict(version=1,family=family,feature='image[N,6,24,24]',keep_image_channels=keep,
        zero_image_channels=[k for k in range(6) if k not in keep],
        anchor_channels=[0,1],broadcast_u2='Original last two normalized base features; unchanged',
        timing='Zero raw feature channels before original train-only normalization; retained-channel statistics unchanged',
        objective='SUPERVISED_RISK',posterior_claim=False)


def masked_features(features,family):
    result=dict(features);result['image']=np.array(features['image'],copy=True)
    result['image'][:,mask_contract(family)['zero_image_channels']]=0
    return result


@contextmanager
def feature_hook(family):
    saved=original.online_features
    def hook(data,geometry): return masked_features(saved(data,geometry),family)
    original.online_features=hook
    try: yield
    finally: original.online_features=saved


class MaskedPredictor:
    def __init__(self,payload,device='cpu'):
        self.p=payload;self.device=torch.device(device)
        family=payload.get('mechanism_family')
        if payload.get('format')!=FORMAT or payload.get('kind')!='mechanism_'+str(family) or payload.get('mask')!=mask_contract(family):
            raise ValueError('explicit mechanism checkpoint/mask version required')
        if payload.get('base_kind')!='local_geometry_cnn': raise ValueError('architecture changes prohibited')
        self.model=original.ResidualNet('local_geometry_cnn',payload['dimensions'],payload['geometry'],payload['residual_scale'],payload['width']).to(self.device)
        self.model.load_state_dict(payload['state_dict']);self.model.eval()
        if original.parameter_count(self.model)!=EXPECTED_PARAMETERS: raise ValueError('parameter count differs from5842')
    def predict(self,online,batch_size=128):
        fields={k:online[k] for k in ('observation','residual','q0','chi0','J')};roi=self.p['geometry']['roi'];predictions=[]
        with torch.inference_mode():
            for begin in range(0,len(fields['q0']),batch_size):
                batch={k:v[begin:begin+batch_size] for k,v in fields.items()};anchor=material_targets(batch['chi0'],roi)
                f=masked_features(original.online_features(batch,self.p['geometry']),self.p['mechanism_family'])
                f=original.normalize(f,self.p['normalization'])
                # Enforce removed channels at inference as well, including any
                # inconsistent normalization in a malformed external payload.
                if np.any(f['image'][:,self.p['mask']['zero_image_channels']]!=0): raise ValueError('removed channels must remain exactly zero')
                args=[torch.from_numpy(f[k]).to(self.device) for k in ('base','geometry','image')]
                args += [torch.as_tensor(batch['q0'],dtype=torch.float32,device=self.device),torch.from_numpy(anchor).to(self.device)]
                y=self.model(*args).cpu().numpy();c=len(roi);predictions.append(y[:,:c]+1j*y[:,c:])
        return np.concatenate(predictions)
    def predict_full(self,online,batch_size=128):
        out=np.asarray(online['chi0']).copy();out[:,self.p['geometry']['roi']]=self.predict(online,batch_size);return out


def load_mechanism_predictor(path,device='cpu'):
    return MaskedPredictor(torch.load(path,map_location='cpu',weights_only=False),device)


def train_masked_residual(family,seed,train,validation,geometry,output,**kwargs):
    output=Path(output)
    with feature_hook(family):
        result=original.train_residual('local_geometry_cnn',seed,train,validation,geometry,output,**kwargs)
    payload=torch.load(output/'BEST.pt',map_location='cpu',weights_only=False)
    payload.update(format=FORMAT,kind='mechanism_'+family,base_kind='local_geometry_cnn',mechanism_family=family,
                   mask=mask_contract(family),objective='SUPERVISED_RISK',parameter_count=EXPECTED_PARAMETERS,
                   loader='p3r2_mechanism.masked_risk.load_mechanism_predictor')
    if result['parameter_count']!=EXPECTED_PARAMETERS: raise ValueError('training architecture parameter count changed')
    torch.save(payload,output/'BEST.pt')
    result.update(kind=family,base_kind='local_geometry_cnn',mechanism_family=family,mask=mask_contract(family),
                  objective='SUPERVISED_RISK',calibrated_density='NOT_CLAIMED',parameter_count=EXPECTED_PARAMETERS,
                  loader='p3r2_mechanism.masked_risk.load_mechanism_predictor',
                  input_contract='Same anchor real/imag and broadcast u2; only declared sensitivity/backprojection channels removed')
    original.write_json(output/'RESULT.json',result)
    # Constructor checks mask metadata, state layout and5842 parameters. No
    # extra training, inference cohort or labels are opened by this check.
    load_mechanism_predictor(output/'BEST.pt')
    return result
