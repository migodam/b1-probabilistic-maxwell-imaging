# R5 lossless wavelet coordinates

Implemented 2026-11-08 in
`research/delegated/z_alg_r2_publication/repository/research/research_gn_nn_feedback_wavelet_r5/src/wavelets.py`
(path relative to Inv_SLAM).

## Interface and dependency

`WaveletPyramid2D(kind='haar', size=24, levels=2)` is a parameter-free PyTorch
module. `analyze(x)` / `forward(x)` return a dictionary; `synthesize(bands)`
reconstructs the field. The final two axes must be 24 by 24; arbitrary leading
batch/channel axes are preserved. `coefficient_vector(x)` concatenates the
spatially flattened bands in the order below, with final dimension 576.
`energy_by_band(error)` returns squared coefficient sums over spatial axes,
preserving leading axes. The returned values are energies, not fractions or
means; summing over bands gives the spatial squared error.

Only Python standard library and existing PyTorch are required. PyWavelets is
not installed and is not required. Constants are constructed in the input's
device/dtype and cached. The derived cache is not a checkpoint dependency;
`state_dict()` requires no wavelet weights. CPU and MPS real float32 paths were
verified. Other real floating dtypes follow native PyTorch device support;
MPS float64 is not claimed.

## Exact coordinate contract

| Order | Name | Spatial shape | Coefficients per field |
|---|---|---|---:|
| 1 | ll2 | 6 by 6 | 36 |
| 2 | lh2 | 6 by 6 | 36 |
| 3 | hl2 | 6 by 6 | 36 |
| 4 | hh2 | 6 by 6 | 36 |
| 5 | lh1 | 12 by 12 | 144 |
| 6 | hl1 | 12 by 12 | 144 |
| 7 | hh1 | 12 by 12 | 144 |
| Total | | | **576** |

The first letter specifies the row-axis filter and the second the column-axis
filter. Thus `lh` means row-low / column-high. These are explicit axis labels,
not an ambiguous horizontal/vertical edge convention. Only LL is recursively
decomposed. All details and the coarsest LL are retained; there is no coefficient
truncation, padding, interpolation, resampling or reduced-dimensional bottleneck.
`levels` may also be another positive integer if `size` is divisible by
`2**levels`; the band names then follow the same coarse-to-fine rule.

## Filters, periodic boundary, and transpose synthesis

For an even one-dimensional period N, the implementation defines correlation
rows

    L[k,j] = sum_t h[t] 1{j = (2k+t) mod N}
    H[k,j] = sum_t g[t] 1{j = (2k+t) mod N}
    A = [L; H],  k = 0,...,N/2-1.

Haar uses `h = [1,1]/sqrt(2)`. Daubechies-2 uses

    h = [1+sqrt(3), 3+sqrt(3), 3-sqrt(3), 1-sqrt(3)] / (4 sqrt(2)).

For both filters of length F,

    g[t] = (-1)**t h[F-1-t].

In particular Haar high pass is `[+1,-1]/sqrt(2)`. Each h has sum sqrt(2)
and unit squared norm, and each g has zero sum and unit squared norm.
These are the usual Haar/db2 orthogonal filter banks with an explicit
correlation/phase choice. The low filter corresponds to the usual
reconstruction-low ordering; coefficient values are not promised to match the
default phase/sign convention of another wavelet package. The official
[PyWavelets filter-bank documentation](https://pywavelets.readthedocs.io/en/latest/ref/wavelets.html)
distinguishes decomposition and reconstruction filter orderings.

At each scale, analysis computes `Y = A X A.T`, then splits Y into its four
quadrants. Synthesis reassembles the quadrants and computes `X = A.T Y A`.
The periodic even-length filter bank is orthogonal, so synthesis is the exact
algebraic transpose/inverse. Floating arithmetic introduces normal rounding
error; “lossless” means complete invertible coordinates, not bitwise identity.

Both axes wrap periodically at each scale (24 and 12 for the default pyramid).
For db2, a filter near an edge therefore reads pixels from the opposite edge.
This is a coordinate boundary assumption, not a change to the Maxwell physical
boundary condition and not evidence that the material is physically periodic.
There is no extension beyond the period and no hidden coefficient cropping.

## One decisive implementation check

Executed once using
`/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python`
with PyTorch **2.14.0**, seed **271828**, float32 random input and independent
probe of shape **(2, 3, 24, 24)**. The same sample was used on CPU and MPS for
both banks. No training, observations or GT were read.

The check measured reconstruction error, coefficient count, total Parseval
energy, autograd of coefficient energy versus `2*x`, and autograd through
analysis -> seven trainable scalar band gains initialized to one -> synthesis
-> dot product with the independent probe. The resulting input gradient was
compared with the probe; every band-gain gradient was finite and nonzero.
All relative checks passed the preselected **3e-6** tolerance.

A focused cache-lifecycle check also passed: first populate db2 constants inside
`torch.inference_mode()`, then use the same module in a normal autograd pass.
The constants are deliberately created outside inference mode so validation
before training cannot leave unusable inference tensors in the cache.

| Device / bank | Roundtrip max abs | Roundtrip relative L2 | Parseval relative error | Energy-gradient relative L2 | End-to-end gradient relative L2 | Smallest absolute band-gain gradient |
|---|---:|---:|---:|---:|---:|---:|
| CPU / Haar | 4.77e-7 | 1.162e-7 | 7.14e-8 | 1.101e-7 | 1.135e-7 | 2.738 |
| CPU / db2 | 5.96e-7 | 1.072e-7 | 7.14e-8 | 1.062e-7 | 1.045e-7 | 6.091 |
| MPS / Haar | 4.77e-7 | 1.145e-7 | 7.14e-8 | 1.148e-7 | 1.158e-7 | 2.738 |
| MPS / db2 | 5.96e-7 | 1.031e-7 | 1.43e-7 | 1.087e-7 | 1.041e-7 | 6.091 |

Each case produced exactly **576 coefficients per field**. These checks establish
the intended transform and gradient wiring at the specified configuration;
they do not establish CNN performance, fixed-point convergence or improved
material recovery. No frozen R4 file was modified and no dependency was added.
