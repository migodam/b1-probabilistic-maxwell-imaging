# R5 local directed-interval certification

## Delivered API

`certification.py` provides:

```python
certify(G, E, S, target, m, box=(0.5, 30.0), timeout=30.0)
```

The returned dictionary is JSON serializable and has status
`CERTIFIED_LOCAL_REAL_ROOT` or `UNKNOWN`. It reports the successful `X`,
interval `F(x)`, interval `J(X)`, `K`, contraction matrix, pivot margins,
inclusion margins, radius, precision, all attempt failures, CPU, and wall time.
The radius schedule is fixed to `[1e-7, 1e-6, 1e-5] * max(1, norm_inf(m))`;
the precision schedule is fixed to 64, 128, and 256 bits.

## Directed proof path

For a real box `X`, the implementation evaluates

```text
A(X) = I - diag(X) G
j(X) = A(X)^-1 diag(X) E
F(X) = S j(X) - target
Etotal(X) = E + G j(X)
J_complex(X)_i,s = (S A(X)^-1)_i Etotal(X)_i,s
```

The complex Jacobian is serialized as real equations in source-major,
real-then-imaginary order. Point values of `G`, `E`, `S`, and `target` are
first interpreted as binary64 and then converted exactly to `mpmath.iv`
intervals. `C` is the binary64 inverse of the point real Jacobian; it is
inserted into the interval calculation as an exact binary64 constant.

`A(X)^-1` is enclosed by interval Gauss-Jordan elimination. Every selected
pivot must exclude zero using endpoint bounds; no float pivot test is used.
The proof accepts only when both

```text
K = x - C F(x) + (I - C J(X))(X - x)
K is strictly inside X
each row sum of abs(I - C J(X)) is strictly less than 1
```

Interval endpoints remain `mpmath.iv` values through all proof checks. They are
converted to binary64 and widened one ULP outward only when producing JSON.

## Regression tests

`test_certification.py` imports the historical R2 physics implementation and
reads, without modifying:

- `experiments/high4_public/E2_HIGH4_00.npz` and both archived roots in
  `results/high4_E2_HIGH4_00_FullJ8.json`;
- `experiments/general6_public/GENERAL6_00.npz` and both archived roots in
  `results/general6_GENERAL6_00_general_quotient20.json`.

All four historical roots certify at the first precision/radius attempt:

| fixture | n | s | radius | max contraction row sum | minimum inclusion margin |
|---|---:|---:|---:|---:|---:|
| E2_HIGH4_00 root 0 | 4 | 2 | 2.6475566059220584e-06 | 0.022562142457780873 | 2.5878220002006204e-06 |
| E2_HIGH4_00 root 1 | 4 | 2 | 2.3722458430409908e-06 | 0.017149948757322304 | 2.3315619025049412e-06 |
| GENERAL6_00 root 0 | 6 | 3 | 2.6113298058117306e-06 | 0.18449574358477103 | 2.129549692346348e-06 |
| GENERAL6_00 root 1 | 6 | 3 | 2.520478733649232e-06 | 0.38438136235341874 | 1.5516525220873876e-06 |

Each fixture also runs one fresh R2 forward and material-adjoint Jacobian.
The fresh response residual is between `4.10e-16` and `1.94e-14`; every R2
point Jacobian entry is contained in the reported interval `J(X)`.

The shifted-target negative control returns `UNKNOWN` and retains all nine
precision/radius failures. The complete test run used `0.327445` CPU seconds
and `0.331998` wall seconds on the recorded local run.

## Cost accounting

For every certified target, the certifier records:

- one interval point inversion and one interval box inversion;
- two complex LU-equivalent binary64 point solves, with `s + 1` complex RHS;
- one real LU-equivalent binary64 point inverse, with `n` real RHS.

Fixture accounting separately records Green setup, receiver SVD, one dense
R2 forward (`1` LU, `s` complex RHS), and its one-receiver material Jacobian
(`1` additional complex RHS). Exact values are in `tests.json`.

## Scope and exclusions

This is a local algebraic result for the scalar-TM map at exact binary64
`G,E,S,target`. It does not certify global root coverage or completeness,
roots outside the reported `X`, complex material roots, Hankel-integral or
geometry accuracy, PDE or discretization error, or physical material
interpretation beyond this algebraic model. A failed attempt is `UNKNOWN`,
not evidence that a root is absent.
