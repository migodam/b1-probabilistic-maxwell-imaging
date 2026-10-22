"""Fixed decoder-envelope sensitivity control; no train-J choice of frame.

Original optimization and CNN stay unchanged. A temporary feature hook removes
both residual backprojections before train-only normalization. The distinct kind
and format require this metadata-aware loader, never the original loader.
"""
from contextlib import contextmanager
from pathlib import Path
import copy
import numpy as np
import torch
from p3r2 import residual_risk as original
from p3r2.risk_training import material_targets

FAMILY='decoder_random_sensitivity'
FORMAT='p3r2_random_sensitivity_risk_v1'
KIND='mechanism_decoder_random_sensitivity'
FRAME_SEED=2026102303
EXPECTED_PARAMETERS=5842
LOADER='p3r2_mechanism.random_sensitivity.load_random_sensitivity_predictor'


def frame():
    q,r=np.linalg.qr(np.random.default_rng(FRAME_SEED).normal(size=(64,64)))
    q*=np.where(np.diag(r)<0,-1.,1.)[None,:]
    return q.astype('float32')


def feature_contract():
    return dict(version=1,family=FAMILY,feature='image[N,6,24,24]',
        keep_image_channels=[0,1,2,3],zero_image_channels=[4,5],anchor_channels=[0,1],
        broadcast_u2='Original last two normalized base features; unchanged',
        frame_seed=FRAME_SEED,frame_shape=[64,64],frame_method='NumPy default_rng Gaussian QR; diagonal R signs canonicalized',
        first_frame_columns=[0,16],second_frame_columns=[16,28],
        sensitivity='abs(delta)*saturation_prime(field/temperature)/temperature*norm(A_i V); original ROI mask',
        frame_fit='Single fixed random orthonormal64 frame; independent of J and all GT; shared across training seeds',
        timing='Zero raw residual-backprojection channels before original train-only normalization',
        objective='SUPERVISED_RISK',posterior_claim=False)


def random_geometry(geometry):
    result=copy.deepcopy(geometry);q=frame()
    if result['weak'].shape!=(64,16) or result['strong'].shape!=(64,12):raise ValueError('fixed16+12 ranks required')
    result.update(weak=q[:,:16].copy(),strong=q[:,16:28].copy(),frame_fit=feature_contract()['frame_fit'])
    return result


def check_geometry(geometry):
    q=frame()
    if not np.array_equal(geometry['weak'],q[:,:16]) or not np.array_equal(geometry['strong'],q[:,16:28]):
        raise ValueError('random sensitivity frame differs from the fixed declared QR frame')
    if geometry.get('frame_fit')!=feature_contract()['frame_fit']:raise ValueError('frame provenance mismatch')


def random_features(features):
    result=dict(features);result['image']=np.array(features['image'],copy=True)
    result['image'][:,4:6]=0
    return result


@contextmanager
def feature_hook():
    saved=original.online_features
    def hook(data,geometry):return random_features(saved(data,geometry))
    original.online_features=hook
    try:yield
    finally:original.online_features=saved


def validate_payload(payload):
    if (payload.get('format')!=FORMAT or payload.get('kind')!=KIND or
        payload.get('mechanism_family')!=FAMILY or payload.get('feature_contract')!=feature_contract() or
        payload.get('base_kind')!='local_geometry_cnn' or payload.get('loader')!=LOADER or
        payload.get('objective')!='SUPERVISED_RISK' or payload.get('parameter_count')!=EXPECTED_PARAMETERS):
        raise ValueError('explicit random-sensitivity checkpoint/version/loader required')
    check_geometry(payload['geometry'])
    state=payload['state_dict']
    for key in ('weak','strong'):
        if not np.array_equal(state[key].cpu().numpy(),payload['geometry'][key]):raise ValueError('state/frame metadata mismatch')
    stats=payload['normalization']
    if np.any(stats['image_mean'][:,4:6]!=0) or np.any(stats['image_scale'][:,4:6]!=1):
        raise ValueError('removed-channel normalization must remain exactly zero/unit')


class RandomSensitivityPredictor:
    def __init__(self,payload,device='cpu'):
        validate_payload(payload);self.p=payload;self.device=torch.device(device)
        self.model=original.ResidualNet('local_geometry_cnn',payload['dimensions'],payload['geometry'],payload['residual_scale'],payload['width']).to(self.device)
        self.model.load_state_dict(payload['state_dict']);self.model.eval()
        if original.parameter_count(self.model)!=EXPECTED_PARAMETERS:raise ValueError('parameter count differs from5842')
    def predict(self,online,batch_size=128):
        fields={k:online[k] for k in ('observation','residual','q0','chi0','J')};roi=self.p['geometry']['roi'];predictions=[]
        with torch.inference_mode():
            for begin in range(0,len(fields['q0']),batch_size):
                batch={k:v[begin:begin+batch_size] for k,v in fields.items()};anchor=material_targets(batch['chi0'],roi)
                f=random_features(original.online_features(batch,self.p['geometry']))
                f=original.normalize(f,self.p['normalization'])
                if np.any(f['image'][:,4:6]!=0):raise ValueError('backprojection channels must remain exactly zero')
                args=[torch.from_numpy(f[k]).to(self.device) for k in ('base','geometry','image')]
                args += [torch.as_tensor(batch['q0'],dtype=torch.float32,device=self.device),torch.from_numpy(anchor).to(self.device)]
                y=self.model(*args).cpu().numpy();c=len(roi);predictions.append(y[:,:c]+1j*y[:,c:])
        return np.concatenate(predictions)
    def predict_full(self,online,batch_size=128):
        out=np.asarray(online['chi0']).copy();out[:,self.p['geometry']['roi']]=self.predict(online,batch_size);return out


def load_random_sensitivity_predictor(path,device='cpu'):
    return RandomSensitivityPredictor(torch.load(path,map_location='cpu',weights_only=False),device)


def train_random_sensitivity(seed,train,validation,geometry,output,**kwargs):
    output=Path(output);geometry=random_geometry(geometry);check_geometry(geometry)
    with feature_hook():
        result=original.train_residual('local_geometry_cnn',seed,train,validation,geometry,output,**kwargs)
    payload=torch.load(output/'BEST.pt',map_location='cpu',weights_only=False)
    payload.update(format=FORMAT,kind=KIND,base_kind='local_geometry_cnn',mechanism_family=FAMILY,
        feature_contract=feature_contract(),objective='SUPERVISED_RISK',parameter_count=EXPECTED_PARAMETERS,loader=LOADER)
    if result['parameter_count']!=EXPECTED_PARAMETERS:raise ValueError('training architecture parameter count changed')
    torch.save(payload,output/'BEST.pt')
    result.update(kind=FAMILY,base_kind='local_geometry_cnn',mechanism_family=FAMILY,feature_contract=feature_contract(),
        objective='SUPERVISED_RISK',calibrated_density='NOT_CLAIMED',parameter_count=EXPECTED_PARAMETERS,loader=LOADER,
        input_contract='Same anchor Re/Im and broadcast u2; fixed random decoder-envelope sensitivity; zero backprojection')
    original.write_json(output/'RESULT.json',result)
    load_random_sensitivity_predictor(output/'BEST.pt')
    # Interface rejection check only: the new kind must never load as an original CNN.
    try:original.load_residual_predictor(output/'BEST.pt')
    except (RuntimeError,ValueError,KeyError):pass
    else:raise ValueError('unmasked original loader incorrectly accepted random sensitivity checkpoint')
    return result
