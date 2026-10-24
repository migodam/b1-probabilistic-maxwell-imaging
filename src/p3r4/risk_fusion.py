"""Frozen spline model moments as CNN risk context, never posterior assertion."""
from pathlib import Path
import numpy as np
import torch
from torch import nn
from .risk import feature_image,roi_components,full_material

FORMAT='p3r4_spline_moment_context_risk_v1';K=32;PERMUTATION_SEED=20261025122
PROVIDERS={'spline':32,'diag_gaussian':32,'field_moments':0}

def moment_contract(value):
    provider=value.get('provider',value.get('moment_provider','spline'))
    count=value.get('K',value.get('moment_K',32))
    if provider not in PROVIDERS or count!=PROVIDERS[provider]:raise ValueError('explicit provider/K contract: spline32,diag_gaussian32,field_moments0 only')
    return provider,int(count)


def matching_moment_gates(train,validation):
    provider,count=moment_contract(train);other=moment_contract(validation)
    identity=train.get('checkpoint_identity')
    if any(g.get('status')!='MOMENTS_COMPLETE' or g.get('truth_read') is not False for g in (train,validation)) or other!=(provider,count) or identity is None or validation.get('checkpoint_identity')!=identity or train.get('sample_seed')!=validation.get('sample_seed'):
        raise ValueError('matching sealed TRAIN/validation provider/checkpoint/K required')
    return provider,count,identity


def parent_seed(parent,seed):return int(np.random.SeedSequence([seed,*str(parent).encode('utf8')]).generate_state(1)[0])


def read_moments(directory,parents,expected_checkpoint=None,expected_provider=None,expected_K=None):
    import json
    from .evaluation import checkpoint_identity
    gate=json.loads((Path(directory)/'MOMENTS_COMPLETE.json').read_text())
    if gate.get('status')!='MOMENTS_COMPLETE' or gate.get('truth_read') is not False:raise ValueError('closed observation-only moment gate required')
    provider,count=moment_contract(gate)
    if expected_provider is not None and provider!=expected_provider:raise ValueError('moment provider mismatch')
    if expected_K is not None and count!=expected_K:raise ValueError('moment K mismatch')
    if expected_checkpoint is not None and gate.get('checkpoint_identity',expected_checkpoint)!=expected_checkpoint:raise ValueError('moment gate checkpoint mismatch')
    sealed={r['path']:r for r in gate['artifacts']}
    rows=[]
    for parent in parents:
        path=Path(directory)/f'{parent}.npz'
        identity_current=checkpoint_identity(path)
        if sealed.get(identity_current['path'])!=identity_current:raise ValueError('moment artifact changed after generation')
        with np.load(Path(directory)/f'{parent}.npz',allow_pickle=False) as z:
            file_provider=str(z['provider'].item()) if 'provider' in z else 'spline'
            if str(z['parent'].item())!=str(parent) or int(z['K'])!=count or file_provider!=provider:raise ValueError('matching parent/provider/K required')
            if provider=='field_moments' and 'x_samples' in z:raise ValueError('deterministic field moments must not invent latent samples')
            image=z['moment_image'];identity=str(z['checkpoint_identity'].item())
        if image.shape!=(4,24,24) or not np.isfinite(image).all() or np.any(image[2:]<-1e-12):raise ValueError('finite material meanRe/Im and varianceRe/Im required')
        if expected_checkpoint is not None:
            import json
            if json.loads(identity)!=expected_checkpoint:raise ValueError('moment source frozen provider identity mismatch')
        rows.append(image)
    return np.asarray(rows,dtype='float32')


def transform_moments(raw,mean,scale,parents,control):
    a=((raw-mean)/scale).astype('float32')
    if control=='permuted':
        for i,parent in enumerate(parents):
            permutation=np.random.default_rng(parent_seed(parent,PERMUTATION_SEED)).permutation(576)
            a[i]=a[i].reshape(4,576)[:,permutation].reshape(4,24,24)
    elif control!='moments':raise ValueError('moments/permuted only')
    return a


class FusionRisk(nn.Module):
    def __init__(self,geometry):
        super().__init__();self.register_buffer('roi',torch.as_tensor(geometry['roi'],dtype=torch.long))
        self.local=nn.Sequential(nn.Conv2d(12,16,3,padding=1),nn.GELU(),nn.Conv2d(16,16,3,padding=1),nn.GELU(),nn.Conv2d(16,16,3,padding=1),nn.GELU())
        self.head=nn.Conv2d(16,2,1);nn.init.zeros_(self.head.weight);nn.init.zeros_(self.head.bias)
    def forward(self,image,moments,anchor):
        if image.shape[1:]!=(8,24,24) or moments.shape!=(len(image),4,24,24):raise ValueError('original8 plus fixed4 moment channels required')
        return anchor+self.head(self.local(torch.cat((image,moments),1))).flatten(2)[:,:,self.roi]


class FusionPredictor:
    def __init__(self,path,device='cpu'):
        self.p=torch.load(path,map_location='cpu',weights_only=False)
        if self.p.get('format')!=FORMAT:raise ValueError('fusion risk format required')
        self.geometry=self.p['geometry'];self.device=device
        self.moment_provider,self.moment_K=moment_contract(self.p)
        self.model=FusionRisk(self.geometry).to(device);self.model.load_state_dict(self.p['state_dict']);self.model.eval()
    def predict(self,y,anchor=None,moments=None,parents=None,batch_size=128,moment_provider=None,moment_K=None):
        if anchor is None or moments is None or parents is None or len(parents)!=len(y):raise ValueError('online momentarrays/parent identities plus original allowed anchor required')
        if self.p.get('moment_contract_version',1)>=2 and (moment_provider is None or moment_K is None):raise ValueError('new checkpoints require explicit inference provider/K metadata')
        if moment_provider is not None and moment_provider!=self.moment_provider:raise ValueError('predictor provider mismatch')
        if moment_K is not None and moment_K!=self.moment_K:raise ValueError('predictor K mismatch')
        normalized=transform_moments(moments,self.p['moment_mean'],self.p['moment_scale'],parents,self.p['control']);rows=[]
        for start in range(0,len(y),batch_size):
            sl=slice(start,start+batch_size);online={k:anchor[k][sl] for k in ('observation','residual','q0','chi0','J')}
            image=feature_image(y[sl],self.geometry,online,self.p['feature_normalization']);base=roi_components(online['chi0'],self.geometry['roi'])
            with torch.inference_mode():result=self.model(torch.tensor(image,device=self.device),torch.tensor(normalized[sl],device=self.device),torch.tensor(base,device=self.device)).cpu().numpy()
            rows.append(full_material(result,self.geometry))
        return np.concatenate(rows)
