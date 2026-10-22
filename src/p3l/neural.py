"""Normalized conditional-density pilot modules; no data or training code.

Inputs are explicit scalar/image contexts; labels are supplied in declared
coordinates (default eight fixed Gaussian weak coordinates plus two probit
properties). Densities are with respect to these input label coordinates.
Parent-owned standardization and property backtransforms require their exact
Jacobians externally. No clipping, target likelihood, truth lookup, posterior
qualification, teacher distillation, or physics operation occurs here.
"""
from math import prod

import torch
from torch import nn
from torch.nn import functional as F


MODEL_TYPES = ('mlp_diag','mlp_full','mlp_mdn','mlp_flow','cnn_full')
MIN_CHOLESKY_DIAGONAL = .02


def parameter_count(model):
    """Number of trainable scalar parameters, including all context branches."""
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


class ContextEncoder(nn.Module):
    """CNN+MLP versus parameter-matched flattened-image/context MLP.

    At defaults both have128 output features, with only208 parameters
    difference. MLP-only still receives all image pixels, so comparisons
    change inductive bias rather than silently remove input information.
    """
    def __init__(self,scalar_dim=1154,image_shape=(4,24,24),*,cnn=False):
        super().__init__()
        self.scalar_dim=int(scalar_dim)
        self.image_shape=tuple(image_shape)
        self.cnn=bool(cnn)
        if self.scalar_dim<1 or self.image_shape!=(4,24,24):
            raise ValueError('positive scalar_dim and image_shape=(4,24,24) required')
        if cnn:
            self.scalar_net=nn.Sequential(nn.Linear(self.scalar_dim,256),nn.SiLU(),
                                          nn.Linear(256,128),nn.SiLU())
            self.image_net=nn.Sequential(nn.Conv2d(4,16,3,stride=2,padding=1),nn.SiLU(),
                nn.Conv2d(16,32,3,stride=2,padding=1),nn.SiLU(),
                nn.Conv2d(32,64,3,stride=2,padding=1),nn.SiLU(),nn.Flatten(),
                nn.Linear(64*3*3,128),nn.SiLU())
            self.fusion=nn.Sequential(nn.Linear(256,128),nn.SiLU())
        else:
            self.fusion=nn.Sequential(nn.Linear(self.scalar_dim+prod(self.image_shape),128),
                                      nn.SiLU(),nn.Linear(128,128),nn.SiLU())

    def forward(self,scalar,image):
        if scalar.ndim!=2 or scalar.shape[1]!=self.scalar_dim:
            raise ValueError('scalar context must have shape (batch,scalar_dim)')
        if image.ndim!=4 or tuple(image.shape[1:])!=self.image_shape or image.shape[0]!=scalar.shape[0]:
            raise ValueError('image context must have shape (batch,4,24,24) matching scalar batch')
        if self.cnn:
            return self.fusion(torch.cat((self.scalar_net(scalar),self.image_net(image)),dim=-1))
        return self.fusion(torch.cat((scalar,image.flatten(start_dim=1)),dim=-1))


class ConditionalDensity(nn.Module):
    def __init__(self,scalar_dim=1154,image_shape=(4,24,24),output_dim=10,*,cnn=False):
        super().__init__()
        if not isinstance(output_dim,int) or output_dim<1:
            raise ValueError('positive integer output_dim required')
        self.output_dim=output_dim
        self.encoder=ContextEncoder(scalar_dim,image_shape,cnn=cnn)

    def parameter_count(self):
        return parameter_count(self)

    def _validate_values(self,values,context):
        if values.ndim not in (2,3) or values.shape[0]!=context.shape[0] or values.shape[-1]!=self.output_dim:
            raise ValueError('values must have shape (batch,d) or (batch,count,d)')
        return values.ndim==3

    @staticmethod
    def _count(count):
        if not isinstance(count,int) or isinstance(count,bool) or count<1:
            raise ValueError('positive integer sample count required')


def _cholesky(raw,dimension):
    """Packed lower triangle with strictly positive softplus+.02 diagonal."""
    rows,cols=torch.tril_indices(dimension,dimension,device=raw.device)
    matrix=raw.new_zeros(*raw.shape[:-1],dimension,dimension)
    matrix[...,rows,cols]=raw
    diagonal=torch.arange(dimension,device=raw.device)
    matrix[...,diagonal,diagonal]=F.softplus(matrix[...,diagonal,diagonal])+MIN_CHOLESKY_DIAGONAL
    return matrix


class ConditionalGaussian(ConditionalDensity):
    def __init__(self,*,full=True,scalar_dim=1154,image_shape=(4,24,24),output_dim=10,cnn=False):
        super().__init__(scalar_dim,image_shape,output_dim,cnn=cnn)
        self.full=bool(full)
        packed=output_dim*(output_dim+1)//2 if full else output_dim
        self.head=nn.Linear(128,output_dim+packed)

    def distribution_parameters(self,scalar,image):
        raw=self.head(self.encoder(scalar,image))
        mean=raw[...,:self.output_dim]
        scale=raw[...,self.output_dim:]
        chol=_cholesky(scale,self.output_dim) if self.full else torch.diag_embed(F.softplus(scale)+MIN_CHOLESKY_DIAGONAL)
        return mean,chol

    def log_prob(self,values,scalar,image):
        mean,chol=self.distribution_parameters(scalar,image)
        extra=self._validate_values(values,mean)
        if extra:mean,chol=mean[:,None,:],chol[:,None,:,:]
        return torch.distributions.MultivariateNormal(mean,scale_tril=chol).log_prob(values)

    def sample(self,count,scalar,image):
        self._count(count)
        mean,chol=self.distribution_parameters(scalar,image)
        return torch.distributions.MultivariateNormal(mean,scale_tril=chol).sample((count,)).transpose(0,1)


class ConditionalMDN(ConditionalDensity):
    """Four normalized full-covariance Gaussian components; no mode semantics."""
    components=4

    def __init__(self,*,scalar_dim=1154,image_shape=(4,24,24),output_dim=10):
        super().__init__(scalar_dim,image_shape,output_dim)
        self.packed=output_dim*(output_dim+1)//2
        self.head=nn.Linear(128,self.components*(1+output_dim+self.packed))

    def distribution_parameters(self,scalar,image):
        raw=self.head(self.encoder(scalar,image)).reshape(-1,self.components,1+self.output_dim+self.packed)
        logits=raw[...,0]
        means=raw[...,1:1+self.output_dim]
        chol=_cholesky(raw[...,1+self.output_dim:],self.output_dim)
        return logits,means,chol

    def log_prob(self,values,scalar,image):
        logits,means,chol=self.distribution_parameters(scalar,image)
        extra=self._validate_values(values,means)
        if extra:
            logits,means,chol=logits[:,None,:],means[:,None,:,:],chol[:,None,:,:,:]
        component=torch.distributions.MultivariateNormal(means,scale_tril=chol).log_prob(values.unsqueeze(-2))
        return torch.logsumexp(F.log_softmax(logits,dim=-1)+component,dim=-1)

    def sample(self,count,scalar,image):
        self._count(count)
        logits,means,chol=self.distribution_parameters(scalar,image)
        # Independent component choices per context and draw; no averaging of
        # means/covariances masquerades as mixture sampling.
        with torch.no_grad():
            component=torch.multinomial(F.softmax(logits,dim=-1),count,replacement=True)
            batch=torch.arange(means.shape[0],device=means.device)[:,None]
            chosen_mean=means[batch,component]
            chosen_chol=chol[batch,component]
            noise=torch.randn_like(chosen_mean)
            return chosen_mean+torch.matmul(chosen_chol,noise.unsqueeze(-1)).squeeze(-1)


class AffineCoupling(nn.Module):
    """Conditional RealNVP bijection, scales bounded by exp(±2)."""
    def __init__(self,dimension,layer_index):
        super().__init__()
        mask=((torch.arange(dimension)+layer_index)%2==0).float()
        self.register_buffer('mask',mask)
        self.net=nn.Sequential(nn.Linear(dimension+128,128),nn.SiLU(),
                               nn.Linear(128,128),nn.SiLU(),nn.Linear(128,2*dimension))

    def _affine_parameters(self,values,context):
        if values.ndim==3:context=context[:,None,:].expand(-1,values.shape[1],-1)
        masked=values*self.mask
        shift,raw_scale=self.net(torch.cat((masked,context),dim=-1)).chunk(2,dim=-1)
        active=1-self.mask
        return shift*active,2*torch.tanh(raw_scale)*active

    def forward_transform(self,values,context):
        shift,log_scale=self._affine_parameters(values,context)
        return values*self.mask+(1-self.mask)*(values*torch.exp(log_scale)+shift),log_scale.sum(dim=-1)

    def inverse_transform(self,values,context):
        shift,log_scale=self._affine_parameters(values,context)
        return values*self.mask+(1-self.mask)*(values-shift)*torch.exp(-log_scale),-log_scale.sum(dim=-1)


class ConditionalRealNVP(ConditionalDensity):
    """Four affine-coupling bijections from a normalized standard Gaussian."""
    def __init__(self,*,scalar_dim=1154,image_shape=(4,24,24),output_dim=10):
        super().__init__(scalar_dim,image_shape,output_dim)
        if output_dim<2:raise ValueError('RealNVP requires output_dim>=2')
        self.layers=nn.ModuleList([AffineCoupling(output_dim,k) for k in range(4)])

    def forward_from_base(self,base,scalar,image):
        context=self.encoder(scalar,image)
        self._validate_values(base,context)
        transformed=base
        logdet=base.new_zeros(base.shape[:-1])
        for layer in self.layers:
            transformed,change=layer.forward_transform(transformed,context)
            logdet=logdet+change
        return transformed,logdet

    def inverse_to_base(self,values,scalar,image):
        context=self.encoder(scalar,image)
        self._validate_values(values,context)
        transformed=values
        logdet=values.new_zeros(values.shape[:-1])
        for layer in reversed(self.layers):
            transformed,change=layer.inverse_transform(transformed,context)
            logdet=logdet+change
        return transformed,logdet

    def log_prob(self,values,scalar,image):
        base,inverse_logdet=self.inverse_to_base(values,scalar,image)
        return -.5*(base.square()+torch.log(base.new_tensor(2*torch.pi))).sum(dim=-1)+inverse_logdet

    def sample(self,count,scalar,image):
        self._count(count)
        base=torch.randn(scalar.shape[0],count,self.output_dim,device=scalar.device,dtype=scalar.dtype)
        with torch.no_grad():
            return self.forward_from_base(base,scalar,image)[0]


def build_model(model_type,*,scalar_dim=1154,image_shape=(4,24,24),output_dim=10):
    """Construct an untrained normalized pilot; no configuration/data lookup."""
    options=dict(scalar_dim=scalar_dim,image_shape=image_shape,output_dim=output_dim)
    if model_type=='mlp_diag':return ConditionalGaussian(full=False,**options)
    if model_type=='mlp_full':return ConditionalGaussian(full=True,**options)
    if model_type=='mlp_mdn':return ConditionalMDN(**options)
    if model_type=='mlp_flow':return ConditionalRealNVP(**options)
    if model_type=='cnn_full':return ConditionalGaussian(full=True,cnn=True,**options)
    raise ValueError(f'unknown model_type; expected one of {MODEL_TYPES}')
