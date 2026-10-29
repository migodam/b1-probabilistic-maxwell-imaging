# 实际运行算法与信息边界

## 1. 物理世界和任务

```
freeze(public geometry, two macro families, two equal-mean textures,
       eta={1/8,1/12,1/16}, train and held-out target parameters,
       static epsinf+i ell/f, nested sources/frequencies/noise)
for each prescribed world and target:
    build public interface-fitted mesh (not material/target dependent)
    assemble LW, bW; cut-integrate VX, bX
    reserve FE; factor LW-VX
    reserve source RHS; solve bW+bX, check each backward residual
    save complex mean and selected field
score pair mean distances and own-worst independent refinements
```

## 2. Measurement-only 有限宏观模型

```
online inputs = y, six public macro wall candidates, 12 public target candidates
joint_weight[g,x] = prior[g,x] * exp(-.5 * ||F(g,x)-y||_Sigma^2)
normalize jointly once
produce target marginal / coherent joint MAP / wall-marginal-selected conditional
wall-first control selects g using target-absent model, then fits x
save all outputs, seal
scorer then opens target/wall labels
known-wall oracle runs only in scorer/explicit diagnostic block
```

Counting prior and finite-grid resolution are explicit. No hidden noise/channel renormalization and no data reuse as independent prior. Continuous FWI uses the same scalar full objective and per-source weighted adjoint; all starts are selected finally by objective, not truth.

## 3. 有限材料–目标补偿

```
fix reference world and dx (offline mechanism, not online reconstruction)
form source-shared actual Maxwell material J using receiver adjoints
for each bounded continuation step:
    solve box-constrained linear correction with numerical damping
    test full nonlinear objective at alpha=1, then at .5 if necessary
    accept only residual decrease
retain two kernel alternatives and all rejected trials
evaluate same endpoint at A0, A1, A2
reoptimize for added data separately; validate selected endpoint independently
```

Numerical damping is not posterior precision. No local minimizer is called a global profile optimum. Later two A2 updates and third mesh were bounded, recorded development refinements within original phase cap.

## 4. 返回传播，完整自散射

```
for each declared wall, common public target ROI P:
    factor L0 and LW; reserve all source and receiver-adjoint actions
    H = P* (GW-G0) P  [operator actions, not full2257 columns]
    Q = qr(H Omega) with fixed seed
    collect H conj(Q)
    core = Q* H conj(Q)
    Hr = Qr core_r Qr^T
for each target:
    exact isolated-target factor (L0-VX)
    q = bX+VX uW
    j0=(I-VX G0)^-1 q; KXQr=(I-VX G0)^-1 VX Qr
    j=j0+KXQr core_r (I-Qr^T KXQr core_r)^-1 Qr^T j0
    output dW+BW j
    paid true-feedback validation; compare direct full on same mesh
```

The initial Hermitian-right projection remains saved. Corrected transpose is reciprocal complex-symmetric compression, not a new Maxwell equation. Q/cores use coefficient Euclidean coordinates; field L2 uses the actual ROI mass matrix.

## 5. E3 consensus logic

```
online gets prior, y, query(atom_id), never table/truth/outerTX
query all anchors plus fixed prior-stratified wall curves:790 unique atoms/case
same nodes -> ordinary response approximation / contrast / sampled consensus
all unqueried cells retain default [0,1] bounds
seal every case before full reference/scoring
measure conditional diameter and missing bounds
if reliable cheap bound missing: no certified early stop, no cold Bayesian run
```

Later micro-group consensus is fully paid, posthoc diagnostic; it does not manufacture saved online queries.
