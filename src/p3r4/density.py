"""Full-dimensional conditional probability models; no Maxwell calls.

All densities in this module use real Gaussianized x=(z64, eta2) coordinates.
Explicit flows have algebraic inverse/Jacobian. CFM log density is deliberately
opt-in and labelled a numerical continuous-flow approximation, never exact MH.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import torch
from torch import nn
from torch.nn import functional as F

LOG2PI = math.log(2 * math.pi)
KINDS = ('gaussian', 'diag_gaussian', 'spline', 'triangular', 'cfm')


def normal_log_prob(x):
    return -.5 * (x.square() + LOG2PI).sum(-1)


def randn(shape, like, generator=None):
    # CPU generators work on both CPU and MPS, and avoid MPS generator gaps.
    if generator is None:
        return torch.randn(shape, device=like.device, dtype=like.dtype)
    return torch.randn(shape, generator=generator, device=generator.device,
                       dtype=like.dtype).to(like.device)


def rational_quadratic_spline(inputs, parameters, inverse=False, bound=5.,
                              min_bin=1e-3, min_derivative=1e-3):
    """Monotone RQ spline on [-bound,bound], identity C1 tails.

    parameters[..., 3*K-1] contains K widths, K heights, K-1 interior
    derivatives. Returns elementwise output and log absolute determinant.
    Boundary derivatives are exactly one. Knots and roots use only PyTorch
    operations supported on CPU/MPS; no searchsorted or external flow library.
    """
    bins = (parameters.shape[-1] + 1) // 3
    if parameters.shape[:-1] != inputs.shape or bins < 2 or parameters.shape[-1] != 3*bins-1:
        raise ValueError('parameters must have shape inputs.shape + (3*K-1,)')
    if bound <= 0 or not 0 < bins*min_bin < 1 or not 0 < min_derivative < 1:
        raise ValueError('invalid spline constraints')
    widths = min_bin + (1-bins*min_bin)*F.softmax(parameters[..., :bins], -1)
    heights = min_bin + (1-bins*min_bin)*F.softmax(parameters[..., bins:2*bins], -1)
    cw = F.pad(torch.cumsum(widths, -1), (1, 0))*2*bound-bound
    ch = F.pad(torch.cumsum(heights, -1), (1, 0))*2*bound-bound
    # Pin exact endpoints; cumsum roundoff must not perturb the tail seam.
    cw = torch.cat((torch.full_like(cw[..., :1], -bound), cw[..., 1:-1],
                    torch.full_like(cw[..., :1], bound)), -1)
    ch = torch.cat((torch.full_like(ch[..., :1], -bound), ch[..., 1:-1],
                    torch.full_like(ch[..., :1], bound)), -1)
    widths, heights = cw[..., 1:]-cw[..., :-1], ch[..., 1:]-ch[..., :-1]
    deriv = min_derivative+F.softplus(parameters[..., 2*bins:])
    deriv = F.pad(deriv, (1, 1), value=1.)
    inside = (inputs >= -bound) & (inputs <= bound)
    safe_input = torch.where(inside, inputs, torch.zeros_like(inputs))
    knots = ch if inverse else cw
    idx = (safe_input[..., None] >= knots[..., 1:]).sum(-1).clamp(max=bins-1)
    def take(t):
        return t.gather(-1, idx[..., None]).squeeze(-1)
    xk, yk, wk, hk = map(take, (cw, ch, widths, heights))
    dk, dk1 = take(deriv), take(deriv[..., 1:])
    delta = hk/wk
    common = dk+dk1-2*delta
    if inverse:
        dy = safe_input-yk
        a = dy*common+hk*(delta-dk)
        b = hk*dk-dy*common
        c = -delta*dy
        discriminant = (b.square()-4*a*c).clamp_min(0)
        # Stable root in [0,1], including a=0 (the identity/linear case).
        theta = (2*c)/(-b-torch.sqrt(discriminant))
        output = xk+theta*wk
    else:
        theta = (safe_input-xk)/wk
        numerator = hk*(delta*theta.square()+dk*theta*(1-theta))
        denominator = delta+common*theta*(1-theta)
        output = yk+numerator/denominator
    denominator = delta+common*theta*(1-theta)
    numerator_derivative = delta.square()*(dk1*theta.square()+
                             2*delta*theta*(1-theta)+dk*(1-theta).square())
    logdet = torch.log(numerator_derivative)-2*torch.log(denominator)
    if inverse:
        logdet = -logdet
    return torch.where(inside, output, inputs), torch.where(inside, logdet, 0.)


class ContextEncoder(nn.Sequential):
    def __init__(self, input_dim=512, width=128):
        super().__init__(nn.Linear(input_dim, width), nn.SiLU(),
                         nn.Linear(width, width), nn.SiLU())


class SplineCoupling(nn.Module):
    def __init__(self, dimension, context_dim, layer, width=128, bins=10,
                 output_rank=32, bound=5.):
        super().__init__()
        if dimension < 2:
            raise ValueError('coupling needs at least two coordinates')
        self.dimension, self.bins, self.bound = dimension, bins, bound
        # Alternate masks and deterministic permutations; every coordinate is
        # transformed repeatedly, no complementary independent-prior block.
        perm = torch.roll(torch.arange(dimension), shifts=2*(layer//2))
        self.register_buffer('keep', perm[(torch.arange(dimension)+layer)%2 == 0])
        self.register_buffer('change', perm[(torch.arange(dimension)+layer)%2 == 1])
        nout = len(self.change)*(3*bins-1)
        self.net = nn.Sequential(nn.Linear(len(self.keep)+context_dim, width), nn.SiLU(),
                                 nn.Linear(width, width), nn.SiLU(),
                                 nn.Linear(width, output_rank), nn.SiLU(),
                                 nn.Linear(output_rank, nout))
        nn.init.zeros_(self.net[-1].weight)
        nn.init.zeros_(self.net[-1].bias)
        with torch.no_grad():
            self.net[-1].bias.reshape(len(self.change), 3*bins-1)[:, 2*bins:] = math.log(math.expm1(1-1e-3))

    def transform(self, x, context, inverse=False):
        params = self.net(torch.cat((x[:, self.keep], context), -1))
        params = params.reshape(len(x), len(self.change), 3*self.bins-1)
        changed, logdet = rational_quadratic_spline(x[:, self.change], params, inverse, self.bound)
        out = x.clone()
        out[:, self.change] = changed
        return out, logdet.sum(-1)


class SplineCore(nn.Module):
    def __init__(self, dimension, context_dim, layers=8, width=128, bins=10,
                 output_rank=32, bound=5.):
        super().__init__()
        self.dimension = dimension
        self.layers = nn.ModuleList([SplineCoupling(dimension, context_dim, i, width,
                                      bins, output_rank, bound) for i in range(layers)])
        self.affine = nn.Linear(context_dim, 2*dimension)
        nn.init.zeros_(self.affine.weight)
        nn.init.zeros_(self.affine.bias)

    def transform(self, x, context, inverse=False):
        mean, logscale = self.affine(context).chunk(2, -1)
        logscale = logscale.clamp(-5, 3)
        logdet = x.new_zeros(len(x))
        if inverse:
            x = (x-mean)*torch.exp(-logscale)
            logdet = logdet-logscale.sum(-1)
        for layer in (reversed(self.layers) if inverse else self.layers):
            x, ld = layer.transform(x, context, inverse)
            logdet = logdet+ld
        if not inverse:
            x = mean+torch.exp(logscale)*x
            logdet = logdet+logscale.sum(-1)
        return x, logdet


class ConditionalDensity(nn.Module):
    """Normalized full-rank Gaussian, full spline, or triangular spline.

    x: [B,66], y: [B,512] (training-standardized context). Triangular is
    q(a|y)q(b|a,y), and samples all b coordinates with full-rank innovation.
    sample returns [B,K,66]. log_prob is with respect to dx.
    """
    def __init__(self, kind='spline', dimension=66, context_dim=512, width=128,
                 layers=8, bins=10, output_rank=32, weak_dim=16, bound=5., diag_hidden=655):
        super().__init__()
        if kind not in ('gaussian', 'diag_gaussian', 'spline', 'triangular') or dimension < 4:
            raise ValueError('use gaussian, diag_gaussian, spline, triangular; dimension >=4')
        self.kind, self.dimension, self.context_dim = kind, dimension, context_dim
        self.encoder = ContextEncoder(context_dim, width)
        if kind == 'diag_gaussian':
            # Same raw-y encoder as the spline; head capacity is explicit.
            # width128, diag_hidden655 => 682943 parameters (spline682860).
            self.head = nn.Sequential(nn.Linear(width, diag_hidden), nn.SiLU(),
                                      nn.Linear(diag_hidden, diag_hidden), nn.SiLU(),
                                      nn.Linear(diag_hidden, 2*dimension))
            nn.init.zeros_(self.head[-1].weight)
            nn.init.zeros_(self.head[-1].bias)
        elif kind == 'gaussian':
            i, j = torch.tril_indices(dimension, dimension)
            self.register_buffer('tril_i', i)
            self.register_buffer('tril_j', j)
            self.head = nn.Linear(width, dimension+len(i))
            nn.init.zeros_(self.head.weight)
            nn.init.zeros_(self.head.bias)
        elif kind == 'spline':
            self.flow = SplineCore(dimension, width, layers, width, bins, output_rank, bound)
        else:
            if not 2 <= weak_dim <= dimension-2 or layers < 4:
                raise ValueError('triangular needs 2<=weak_dim<=D-2 and layers>=4')
            self.weak_dim = weak_dim
            opts = dict(width=width, bins=bins, output_rank=output_rank, bound=bound)
            self.weak = SplineCore(weak_dim, width, layers//2, **opts)
            self.strong = SplineCore(dimension-weak_dim, width+weak_dim, layers-layers//2, **opts)

    def gaussian_parameters(self, context):
        raw = self.head(context)
        mean, packed = raw[:, :self.dimension], raw[:, self.dimension:]
        diagonal = self.tril_i == self.tril_j
        values = .5*torch.tanh(packed)
        values = torch.where(diagonal[None, :], torch.exp(packed.clamp(-5, 3)), values)
        chol = raw.new_zeros(len(raw), self.dimension, self.dimension)
        chol[:, self.tril_i, self.tril_j] = values
        return mean, chol

    def transform(self, x, y, inverse=False):
        if x.shape != (len(y), self.dimension) or y.shape[1:] != (self.context_dim,):
            raise ValueError('aligned x[B,D], y[B,context_dim] required')
        context = self.encoder(y)
        if self.kind == 'diag_gaussian':
            mean, logscale = self.head(context).chunk(2, -1)
            logscale = logscale.clamp(-5, 3)
            ld = logscale.sum(-1)
            if inverse:
                return (x-mean)*torch.exp(-logscale), -ld
            return mean+torch.exp(logscale)*x, ld
        if self.kind == 'gaussian':
            mean, chol = self.gaussian_parameters(context)
            ld = torch.log(torch.diagonal(chol, dim1=-2, dim2=-1)).sum(-1)
            if inverse:
                return torch.linalg.solve_triangular(chol, (x-mean)[..., None], upper=False).squeeze(-1), -ld
            return mean+(chol@x[..., None]).squeeze(-1), ld
        if self.kind == 'spline':
            return self.flow.transform(x, context, inverse)
        k = self.weak_dim
        if inverse:
            a, b = x[:, :k], x[:, k:]
            zb, ldb = self.strong.transform(b, torch.cat((context, a), -1), True)
            za, lda = self.weak.transform(a, context, True)
            return torch.cat((za, zb), -1), lda+ldb
        a, lda = self.weak.transform(x[:, :k], context)
        b, ldb = self.strong.transform(x[:, k:], torch.cat((context, a), -1))
        return torch.cat((a, b), -1), lda+ldb

    def log_prob(self, x, y):
        z, ld = self.transform(x, y, True)
        return normal_log_prob(z)+ld

    def sample(self, y, n_samples=1, generator=None):
        if n_samples < 1:
            raise ValueError('positive n_samples required')
        c = y.repeat_interleave(n_samples, 0)
        z = randn((len(c), self.dimension), c, generator)
        return self.transform(z, c)[0].reshape(len(y), n_samples, self.dimension)


class VelocityBlock(nn.Module):
    def __init__(self, width, hidden=None):
        super().__init__()
        hidden = width if hidden is None else hidden
        self.net = nn.Sequential(nn.Linear(width, hidden), nn.SiLU(), nn.Linear(hidden, width))
    def forward(self, h):
        return h+self.net(h)/2


class ConditionalFlowMatching(nn.Module):
    """Independent conditional CFM with non-degenerate Gaussian base.

    base_mean/base_chol optionally supply y-only Gaussian parameters. Without
    them buffers base_location/base_scale define a diagonal, full-rank base.
    The base must be frozen and independent of hidden endpoint x1 given y.
    No endpoint smoothing, OT rematching, projection, or physical refresh.
    """
    kind = 'cfm'
    def __init__(self, dimension=66, context_dim=512, width=128, blocks=4,
                 time_features=32, velocity_hidden=None, **unused):
        super().__init__()
        if time_features < 2 or time_features % 2:
            raise ValueError('even positive time_features required')
        self.dimension, self.context_dim = dimension, context_dim
        self.encoder = ContextEncoder(context_dim, width)
        self.register_buffer('frequencies', 2**torch.arange(time_features//2, dtype=torch.float32)*math.pi)
        self.register_buffer('base_location', torch.zeros(dimension))
        self.register_buffer('base_scale', torch.ones(dimension))
        self.input = nn.Linear(dimension+width+time_features, width)
        self.blocks = nn.ModuleList([VelocityBlock(width, velocity_hidden) for _ in range(blocks)])
        self.output = nn.Linear(width, dimension)
        nn.init.zeros_(self.output.weight)
        nn.init.zeros_(self.output.bias)

    def velocity(self, x, t, y):
        if not torch.is_tensor(t):
            t = x.new_full((len(x),), t)
        t = t.reshape(-1, 1).expand(len(x), 1)
        phase = t*self.frequencies[None, :]
        h = F.silu(self.input(torch.cat((x, self.encoder(y), phase.sin(), phase.cos()), -1)))
        for block in self.blocks:
            h = block(h)
        return self.output(F.silu(h))

    def base_sample(self, y, generator=None, base_mean=None, base_chol=None):
        if (base_mean is None) != (base_chol is None):
            raise ValueError('both conditional base_mean and base_chol required')
        eps = randn((len(y), self.dimension), y, generator)
        if base_mean is None:
            return self.base_location+self.base_scale*eps
        self._check_base(base_mean, base_chol, len(y))
        return base_mean+(base_chol@eps[..., None]).squeeze(-1)

    def _check_base(self, mean, chol, batch):
        if mean.shape != (batch, self.dimension) or chol.shape != (batch, self.dimension, self.dimension):
            raise ValueError('full Gaussian base shapes [B,D], [B,D,D] required')
        if not bool(torch.isfinite(chol).all()) or not bool(torch.isfinite(mean).all()):
            raise ValueError('finite base parameters required')
        if not bool((torch.diagonal(chol, dim1=-2, dim2=-1)>0).all()) or not torch.allclose(chol, torch.tril(chol)):
            raise ValueError('positive-diagonal lower-triangular full-rank base required')

    def base_log_prob(self, x, base_mean=None, base_chol=None):
        if base_mean is None and base_chol is None:
            return normal_log_prob((x-self.base_location)/self.base_scale)-self.base_scale.log().sum()
        self._check_base(base_mean, base_chol, len(x))
        eps = torch.linalg.solve_triangular(base_chol, (x-base_mean)[..., None], upper=False).squeeze(-1)
        return normal_log_prob(eps)-torch.diagonal(base_chol, dim1=-2, dim2=-1).log().sum(-1)

    def loss(self, x1, y, generator=None, base_mean=None, base_chol=None):
        x0 = self.base_sample(y, generator, base_mean, base_chol)
        if generator is None:
            t = torch.rand((len(y), 1), device=y.device, dtype=y.dtype)
        else:
            t = torch.rand((len(y), 1), generator=generator, device=generator.device,
                           dtype=y.dtype).to(y.device)
        xt = (1-t)*x0+t*x1
        return (self.velocity(xt, t, y)-(x1-x0)).square().mean()

    def integrate(self, x, y, nfe=32, reverse=False, divergence=False):
        """Fixed RK4, exactly nfe velocity calls; full trace if requested.

        Continuous-flow density is numerical. No stochastic trace is used.
        Full-trace audit requires gradients and is intentionally expensive.
        """
        if nfe < 4 or nfe % 4:
            raise ValueError('RK4 nfe must be positive multiple of four')
        dt = (-1. if reverse else 1.)/(nfe//4)
        t = 1. if reverse else 0.
        accumulated = x.new_zeros(len(x))
        def rhs(state, time):
            if not divergence:
                return self.velocity(state, time, y), state.new_zeros(len(state))
            with torch.enable_grad():
                state = state.detach().requires_grad_(True)
                v = self.velocity(state, time, y)
                tr = state.new_zeros(len(state))
                for d in range(self.dimension):
                    tr = tr+torch.autograd.grad(v[:, d].sum(), state, retain_graph=d<self.dimension-1)[0][:, d]
            return v.detach(), -tr.detach()
        for _ in range(nfe//4):
            k1, d1 = rhs(x, t)
            k2, d2 = rhs(x+.5*dt*k1, t+.5*dt)
            k3, d3 = rhs(x+.5*dt*k2, t+.5*dt)
            k4, d4 = rhs(x+dt*k3, t+dt)
            x = x+dt*(k1+2*k2+2*k3+k4)/6
            accumulated = accumulated+dt*(d1+2*d2+2*d3+d4)/6
            t += dt
        return x, accumulated

    def sample(self, y, n_samples=1, generator=None, nfe=32, base_mean=None, base_chol=None):
        if n_samples < 1:
            raise ValueError('positive n_samples required')
        c = y.repeat_interleave(n_samples, 0)
        mean = base_mean.repeat_interleave(n_samples, 0) if base_mean is not None else None
        chol = base_chol.repeat_interleave(n_samples, 0) if base_chol is not None else None
        x0 = self.base_sample(c, generator, mean, chol)
        with torch.no_grad():
            result, _ = self.integrate(x0, c, nfe=nfe)
        return result.reshape(len(y), n_samples, self.dimension)

    def log_prob(self, x, y, *, allow_approximate=False, nfe=64, base_mean=None, base_chol=None):
        if not allow_approximate:
            raise ValueError('CFM density is numerical continuous-flow approximation; explicitly set allow_approximate=True; unsuitable for exact MH')
        with torch.no_grad():
            x0, delta = self.integrate(x, y, nfe=nfe, reverse=True, divergence=True)
        return self.base_log_prob(x0, base_mean, base_chol)-delta


def make_model(kind, **config):
    return ConditionalFlowMatching(**config) if kind == 'cfm' else ConditionalDensity(kind, **config)


class DensityPredictor:
    """Saved-model inference on raw Gaussianized x and source-major real y.

    Training-only normalizers and a frozen orthogonal frame are restored.
    Bounded chunks are used for all inference. Sampling never reads truth.
    """
    def __init__(self, payload, device='cpu'):
        self.payload = payload
        self.device = torch.device(device)
        self.model = make_model(payload['kind'], **payload['model_config']).to(self.device)
        self.model.load_state_dict(payload['state_dict'])
        self.model.eval()
        self.stats = payload['normalization']

    def _y(self, y):
        y = np.asarray(y)
        if y.ndim != 2 or y.shape[1] != self.model.context_dim or np.iscomplexobj(y) or not np.isfinite(y).all():
            raise ValueError('finite real source-major y[B,512] required')
        return np.asarray((y-self.stats['y_mean'])/self.stats['y_scale'], dtype=np.float32)

    def _base(self, base_mean, base_chol, count):
        required = self.payload.get('contract', {}).get('cfm_base') == 'anchor'
        if base_mean is None and base_chol is None:
            if required:
                raise ValueError('checkpoint requires its y-only anchor base_mean and fullrank base_chol')
            return None, None
        if self.model.kind != 'cfm' or base_mean is None or base_chol is None:
            raise ValueError('conditional Gaussian base requires CFM and both mean/chol')
        mean, chol = np.asarray(base_mean), np.asarray(base_chol)
        d = self.model.dimension
        if mean.shape != (count, d) or chol.shape != (count, d, d):
            raise ValueError('raw conditional base shapes [B,D], [B,D,D] required')
        if not np.isfinite(mean).all() or not np.isfinite(chol).all() or np.iscomplexobj(mean) or np.iscomplexobj(chol):
            raise ValueError('finite real conditional base required')
        if np.any(np.diagonal(chol, axis1=1, axis2=2) <= 0) or not np.allclose(chol, np.tril(chol)):
            raise ValueError('positive-diagonal lower-triangular conditional base required')
        mean = ((mean@self.stats['frame']-self.stats['x_mean'])/self.stats['x_scale']).astype(np.float32)
        transform = self.stats['frame'].T/self.stats['x_scale'][:, None]
        converted = np.empty_like(chol, dtype=np.float32)
        for begin in range(0, count, 128):
            a = transform@chol[begin:begin+128]
            converted[begin:begin+128] = np.linalg.cholesky(a@a.transpose(0, 2, 1))
        return mean, converted

    def log_prob(self, x, y, batch_size=128, base_mean=None, base_chol=None, **kwargs):
        x, cy = np.asarray(x), self._y(y)
        if x.shape != (len(cy), self.model.dimension) or np.iscomplexobj(x) or not np.isfinite(x).all():
            raise ValueError('finite aligned Gaussianized x[B,66] required')
        if batch_size < 1:
            raise ValueError('positive batch_size required')
        frame = self.stats['frame']
        sx = ((x@frame-self.stats['x_mean'])/self.stats['x_scale']).astype(np.float32)
        mean, chol = self._base(base_mean, base_chol, len(x))
        values = []
        for start in range(0, len(x), batch_size):
            tx = torch.from_numpy(sx[start:start+batch_size]).to(self.device)
            ty = torch.from_numpy(cy[start:start+batch_size]).to(self.device)
            base_kwargs = {} if mean is None else dict(
                base_mean=torch.from_numpy(mean[start:start+batch_size]).to(self.device),
                base_chol=torch.from_numpy(chol[start:start+batch_size]).to(self.device))
            with torch.no_grad():
                lp = self.model.log_prob(tx, ty, **base_kwargs, **kwargs)
            values.append(lp.detach().cpu().numpy())
        return (np.concatenate(values).astype(np.float64)-np.log(self.stats['x_scale']).sum()
                +self.stats.get('frame_logabsdet', 0.)) if values else np.empty(0)

    def sample(self, y, n_samples=1, seed=0, batch_size=128, base_mean=None, base_chol=None, **kwargs):
        cy = self._y(y)
        if batch_size < 1 or n_samples < 1:
            raise ValueError('positive batch_size and n_samples required')
        generator = torch.Generator(device='cpu').manual_seed(seed)
        mean, chol = self._base(base_mean, base_chol, len(cy))
        # Limit effective network batch even when K exceeds batch_size.
        chunks = []
        for start in range(0, len(cy), max(1, batch_size//n_samples)):
            rows = cy[start:start+max(1, batch_size//n_samples)]
            base_kwargs = {} if mean is None else dict(
                base_mean=torch.from_numpy(mean[start:start+len(rows)]).to(self.device),
                base_chol=torch.from_numpy(chol[start:start+len(rows)]).to(self.device))
            # If one observation requests many samples, split sample count too.
            parts = []
            for offset in range(0, n_samples, batch_size):
                k = min(batch_size, n_samples-offset)
                with torch.no_grad():
                    value = self.model.sample(torch.from_numpy(rows).to(self.device), k, generator, **base_kwargs, **kwargs)
                parts.append(value.detach().cpu().numpy())
            chunks.append(np.concatenate(parts, axis=1))
        if not chunks:
            return np.empty((0, n_samples, self.model.dimension))
        sx = np.concatenate(chunks, axis=0)
        return (sx*self.stats['x_scale']+self.stats['x_mean'])@self.stats.get('frame_inverse', self.stats['frame'].T)


def load_density(path, device='cpu'):
    return DensityPredictor(torch.load(Path(path), map_location='cpu', weights_only=False), device)
