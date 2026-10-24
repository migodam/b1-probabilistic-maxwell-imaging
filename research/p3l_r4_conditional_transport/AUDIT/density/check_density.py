"""Zero-physics correctness audit. No campaign data and no optimizer updates."""
import json
from pathlib import Path
import resource
import sys
import time

import numpy as np
import torch

from src.p3r4.density import (ConditionalDensity, ConditionalFlowMatching,
    rational_quadratic_spline, make_model, normal_log_prob, DensityPredictor)
from src.p3r4.train_density import fit_normalization

wall, cpu = time.monotonic(), time.process_time()
torch.set_num_threads(1)
torch.manual_seed(139)
results = {}

# Random RQ splines, identity tails, exact boundary continuity and gradients.
p = torch.randn(29, 5, 29, dtype=torch.float64)*1.1
x = torch.linspace(-7, 7, 145, dtype=torch.float64).reshape(29, 5)
y, ld = rational_quadratic_spline(x, p)
z, ild = rational_quadratic_spline(y, p, inverse=True)
results['rqs_roundtrip_max'] = float((x-z).abs().max())
results['rqs_logdet_cancel_max'] = float((ld+ild).abs().max())
assert results['rqs_roundtrip_max'] < 1e-9
assert results['rqs_logdet_cancel_max'] < 1e-9
inside = x.abs() <= 5
assert torch.equal(y[~inside], x[~inside])
xi = x.clone().requires_grad_(True)
yi, ldi = rational_quadratic_spline(xi, p)
real_ld = torch.autograd.grad(yi.sum(), xi)[0].log()
results['rqs_autodiff_logdet_max'] = float((real_ld-ldi).abs().max())
assert results['rqs_autodiff_logdet_max'] < 1e-9
# The inverse also differentiates correctly, including exactly initialized flow.
identity = torch.zeros_like(p)
identity[..., 20:] = np.log(np.expm1(1-1e-3))
zi, zld = rational_quadratic_spline(x, identity, inverse=True)
assert torch.allclose(zi, x, atol=1e-12)
for edge in (-5., 5.):
    ei = torch.full((29, 5), edge, dtype=torch.float64, requires_grad=True)
    eo, el = rational_quadratic_spline(ei, p)
    assert torch.allclose(eo, ei, atol=1e-12)
    assert torch.allclose(torch.autograd.grad(eo.sum(), ei)[0], torch.ones_like(ei), atol=1e-10)
    assert torch.allclose(el, torch.zeros_like(el), atol=1e-10)

# Numerically integrate a 1D transformed standard normal over almost all mass.
grid = torch.linspace(-10, 10, 40001, dtype=torch.float64)
param = torch.randn(1, 29, dtype=torch.float64).expand(len(grid), 29)
base, grid_ild = rational_quadratic_spline(grid, param, inverse=True)
q = torch.exp(-.5*(base.square()+np.log(2*np.pi))+grid_ild)
mass = torch.trapezoid(q, grid)
results['rqs_integrated_probability_mass'] = float(mass)
assert abs(float(mass)-1) < 2e-5

# Full66 inverse, Jacobian, all-dimensional change and NLL gradient tests.
for kind in ('gaussian', 'diag_gaussian', 'spline', 'triangular'):
    model = ConditionalDensity(kind, width=32, layers=8, output_rank=12).double()
    c = torch.randn(2, 512, dtype=torch.float64)
    z = torch.randn(2, 66, dtype=torch.float64)
    # Initialization is a normalized unit Gaussian in all66 coordinates.
    assert torch.allclose(model.log_prob(z, c), normal_log_prob(z), atol=1e-10)
    # Perturb every last-layer map so tests are not identity-only.
    with torch.no_grad():
        if kind == 'diag_gaussian':
            model.head[-1].weight.normal_(0, .02)
            model.head[-1].bias.normal_(0, .02)
        elif kind == 'gaussian':
            model.head.weight.normal_(0, .02)
            model.head.bias.normal_(0, .02)
        else:
            cores = [model.flow] if kind == 'spline' else [model.weak, model.strong]
            for core in cores:
                core.affine.weight.normal_(0, .02)
                for layer in core.layers:
                    layer.net[-1].weight.normal_(0, .02)
                    layer.net[-1].bias.add_(torch.randn_like(layer.net[-1].bias)*.07)
            assert all(all(any(int(i) in layer.change.tolist() for layer in core.layers)
                       for i in range(core.dimension)) for core in cores)
    mapped, forward_ld = model.transform(z, c)
    recovered, inverse_ld = model.transform(mapped, c, inverse=True)
    err = float((recovered-z).abs().max())
    cancel = float((forward_ld+inverse_ld).abs().max())
    assert err < 1e-9 and cancel < 1e-9, (kind, err, cancel)
    ci, zi = c[:1], z[0].detach().requires_grad_(True)
    jac = torch.autograd.functional.jacobian(lambda v: model.transform(v[None], ci)[0][0], zi)
    actual_ld = torch.linalg.slogdet(jac)[1]
    analytic_ld = model.transform(zi[None], ci)[1][0]
    jacerr = abs(float(actual_ld-analytic_ld))
    assert jacerr < 1e-8, (kind, jacerr)
    (-model.log_prob(mapped.detach(), c).mean()).backward()
    assert all(torch.isfinite(t.grad).all() for t in model.parameters() if t.grad is not None)
    samples = model.sample(c, 7, torch.Generator().manual_seed(1))
    assert samples.shape == (2, 7, 66) and torch.isfinite(samples).all()
    results[kind] = dict(inverse_max=err, logdet_cancellation_max=cancel,
                         full66_autodiff_logdet_error=jacerr, gradients_finite=True)
    if kind == 'gaussian':
        loc, chol = model.gaussian_parameters(model.encoder(c))
        reference = torch.distributions.MultivariateNormal(loc, scale_tril=chol).log_prob(mapped)
        assert torch.allclose(reference, model.log_prob(mapped, c), atol=1e-10)
    if kind == 'triangular':
        # Strong conditional must actually depend on sampled a; its Jacobian
        # includes this shear while determinant still factors block-triangularly.
        z2 = z.clone()
        z2[:, :model.weak_dim] += .5
        assert float((model.transform(z2, c)[0][:, model.weak_dim:]-mapped[:, model.weak_dim:]).abs().max()) > 1e-5

# Train-only statistics include correct frame and x Jacobian in raw dx66.
rng = np.random.default_rng(12)
frame = np.linalg.qr(rng.normal(size=(66, 66)))[0]
rawx = rng.normal(size=(40, 66))*np.linspace(.4, 2, 66)+.3
rawy = rng.normal(size=(40, 512))*3
stats = fit_normalization(rawx, rawy, frame)
m = ConditionalDensity('gaussian', width=16)
payload = dict(kind='gaussian', model_config=dict(width=16), state_dict=m.state_dict(), normalization=stats)
predictor = DensityPredictor(payload)
rawlp = predictor.log_prob(rawx, rawy, batch_size=7)
stdx = (rawx@frame-stats['x_mean'])/stats['x_scale']
reference = -.5*(stdx**2+np.log(2*np.pi)).sum(1)-np.log(stats['x_scale']).sum()
assert np.max(abs(rawlp-reference)) < 2e-5
samples = predictor.sample(rawy[:3], n_samples=145, batch_size=32)
assert samples.shape == (3, 145, 66) and np.isfinite(samples).all()
results['raw_coordinate_logdet_max'] = float(np.max(abs(rawlp-reference)))

# CFM has a known fullrank Gaussian source and no endpoint smoothing. Audit an
# analytic linear velocity v_i=rate_i*x_i + shift_i, whose exact density is known.
cfm = ConditionalFlowMatching(width=16).double()
rate = torch.linspace(-.3, .3, 66, dtype=torch.float64)
shift = torch.linspace(-.2, .2, 66, dtype=torch.float64)
cfm.velocity = lambda x, t, y: rate*x+shift
c = torch.zeros(2, 512, dtype=torch.float64)
z = torch.randn(2, 66, dtype=torch.float64)
exact = z*rate.exp()+shift*torch.expm1(rate)/rate
for nfe in (16, 32, 64):
    mapped, delta = cfm.integrate(z, c, nfe=nfe, divergence=True)
    inverse, _ = cfm.integrate(mapped, c, nfe=nfe, reverse=True)
    lp = cfm.log_prob(mapped, c, allow_approximate=True, nfe=nfe)
    target_lp = normal_log_prob(z)-rate.sum()
    error = float((lp-target_lp).abs().max())
    results[f'cfm_nfe{nfe}'] = dict(sample_max=float((mapped-exact).abs().max()),
                                    inverse_max=float((inverse-z).abs().max()),
                                    log_density_max=error)
    assert error < 1e-6 and float((mapped-exact).abs().max()) < 2e-6
try:
    cfm.log_prob(z, c)
    raise AssertionError('CFM density silently presented as exact')
except ValueError:
    pass
mean = torch.randn(2, 66, dtype=torch.float64)
chol = torch.eye(66, dtype=torch.float64)[None].repeat(2, 1, 1)*.3
base = cfm.base_sample(c, torch.Generator().manual_seed(3), mean, chol)
ref = torch.distributions.MultivariateNormal(mean, scale_tril=chol).log_prob(base)
assert torch.allclose(cfm.base_log_prob(base, mean, chol), ref, atol=1e-10)

# CPU/MPS real model gradients and inversion; no optimizer or physics calls.
for device in ('cpu', 'mps'):
    if device == 'mps' and not torch.backends.mps.is_available():
        results['mps'] = 'NOT_AVAILABLE'
        continue
    checks = {}
    for kind in ('gaussian', 'diag_gaussian', 'spline', 'triangular', 'cfm'):
        model = make_model(kind, width=32, layers=8, output_rank=12).to(device)
        x, y = torch.randn(4, 66, device=device), torch.randn(4, 512, device=device)
        if kind == 'cfm':
            loss = model.loss(x, y, torch.Generator().manual_seed(4))
        else:
            loss = -model.log_prob(x, y).mean()
        loss.backward()
        assert torch.isfinite(loss) and all(torch.isfinite(p.grad).all() for p in model.parameters() if p.grad is not None)
        samples = model.sample(y, 2, torch.Generator().manual_seed(4))
        assert torch.isfinite(samples).all()
        checks[kind] = dict(loss=float(loss.detach().cpu()), gradients_finite=True)
        if kind != 'cfm':
            with torch.no_grad():
                for name, parameter in model.named_parameters():
                    if name.endswith(('head.weight', 'head.4.weight', 'net.6.weight', 'affine.weight')):
                        parameter.add_(torch.randn_like(parameter)*.02)
                mapped, ld = model.transform(x, y)
                back, ild = model.transform(mapped, y, inverse=True)
                err = float((back-x).abs().max().cpu())
                cancel = float((ld+ild).abs().max().cpu())
                assert err < 2e-4 and cancel < 2e-4, (device, kind, err, cancel)
                checks[kind].update(perturbed_inverse_max=err, perturbed_logdet_cancellation_max=cancel)
        elif device == 'mps':
            lp = model.log_prob(x[:1], y[:1], allow_approximate=True, nfe=4)
            assert torch.isfinite(lp).all()
            checks[kind]['full_trace_mps_finite'] = True
    results[device] = checks

for kind in ('gaussian', 'diag_gaussian', 'spline', 'triangular', 'cfm'):
    m = make_model(kind)
    results.setdefault('default_parameter_counts', {})[kind] = sum(t.numel() for t in m.parameters())
results.update(status='PASS', physics_calls=0, optimizer_updates=0,
               wall_seconds=time.monotonic()-wall, cpu_seconds=time.process_time()-cpu,
               peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),
               torch_version=torch.__version__)
output = Path(__file__).parent/'NUMERICAL_CORRECTNESS.json'
output.write_text(json.dumps(results, indent=2)+'\n')
print(json.dumps(results, indent=2))
