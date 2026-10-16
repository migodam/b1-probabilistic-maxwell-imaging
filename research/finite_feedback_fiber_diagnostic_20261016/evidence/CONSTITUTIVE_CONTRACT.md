# Read-only constitutive and packing evidence

No model imports, physics runs, tests, new measurements, or historical changes were performed. Paths below are relative to `/Volumes/migodam's-external-brain/Research/Inv_SLAM`.

## A30 exact convention (decisive)

Let `P=Gaussian/A30/a30-maxwell-material-manifold/src/a30/physics.py`, `V=Gaussian/A30/a30-maxwell-material-manifold/vendor/frozen_a29/src/a24/physics.py`, and `C=Gaussian/A30/a30-maxwell-material-manifold/vendor/frozen_a29/src/a24/contracts.py`.

- A30 explicitly imports the vendored A24 adapter (P:59–62), using scalar TM, 24×24 cells, volume 1/576, 8 sources, 32 receivers, `scalar_chi`, and `sqrt_cell_volume` normalization (P:145–152, 175–177, 239–257).
- **Here X=diag(chi), exactly.** V:119–124 returns `a=chi` and derivative identically one. V:191–197 constructs `L=I-diag(chi)G`, solves `L current=diag(chi) incident`, sets `field=incident+G current`, and `data=S current`. P:381–386 independently audits these same identities, including relative `w=(chi-chi_b)field`.
- V:3–4 documents mass normalized current `c=sqrt(v)j` and exciting field `e=sqrt(v)E+Gc`. Thus the user’s symbol `j` must mean the adapter’s normalized `current` for direct reuse; `Etot` must mean its normalized `field`. Do not insert an extra volume or 3-vector repeat.
- Native arrays are complex128 `current,field:(576,8)`, `data:(32,8)`, `chi:(576,)` (P:151,196–198). Material real mass coordinates are `sqrt(v)[all Re chi; all Im chi]` (C:115–123). Data real coordinates are **per source** `[Re all 32 receivers; Im all 32 receivers]`, then next source (C:101–113), not global all-data-real then all-data-imaginary.
- Domain check allows real chi in [-.5,4] and imaginary chi in [0,.5], tolerance 1e-10 (V:107–117). An increment can be signed; the endpoint material must satisfy these bounds.
- Existing exact Jacobian is `S L0^-1 diag(field0) dchi`, independently for each source (P:431–440). Cached transfer `T=S L0^-1` costs q adjoint RHS on first construction (P:442–451). Material pullback and mass-coordinate factor `1/sqrt(v)` are explicit at P:467–475 and V:210–220. The saved background real mass J construction in P:269–276 uses the same source packing.
- `L_response` is the **known-background affine response** `yb+Db((chi-chi_b)eb)` (P:406–411); it is not the current-material full operator named L in `L=I-XG`.

## Exact finite increment identity to reuse

For an admitted A30 anchor, put `X0=diag(chi0)`, `DeltaX=diag(chi1-chi0)`, `R0=L0^-1`, and `e0=incident+G current0`. Subtracting the two source equations gives

`(I-R0 DeltaX G) DeltaCurrent = R0 DeltaX e0`,

hence

`DeltaCurrent=(I-R0 DeltaX G)^-1 R0 DeltaX e0`, and `DeltaData=S DeltaCurrent`.

This is an algebraic inference from the cited equations, not a newly measured physics result. It requires invertibility of the endpoint operator; no Neumann smallness is needed for the exact inverse, while a Neumann truncation would need its own justification. `R0 DeltaX G` ordering is essential. A30 `solve_anchor` already applies R0 with a residual audit (P:413–429); explicit inverse construction is unnecessary. `anchor_jvp` is the first-order term only, not the finite feedback sum.

## A17 closure_final: related machinery, different material law

Let `F=Gaussian/A17/CLOSURE_R1/remote_snapshots/closure_final_closed_v1/source`, `E=Gaussian/A17/EXCHANGE_R1`, `D=E/vendor/Gaussian/A9/code/a9_engine.py`, `O=E/vendor/Gaussian/A16/CURRENT_SELECTION_R1/code/operators.py`, `M=E/vendor/Gaussian/A16/CURRENT_SELECTION_R1/code/material_gauge_persistence.py`.

- F/research/A17_BASELINE_LOCK.json:4 identifies E as baseline; F/code/closure_runtime.py:31–46 loads that baseline’s declared modules. F/code/closure_state_v2.py:40–68 uses `DenseDDA`, saved material gauge, 3D volume `(edge/n)^3`, and `model.state(anchor chi)`.
- A17 uses **CM plus radiation reaction polarizability**, `a(chi)=3v chi/(chi+3-i3 c_k v chi)` and `da=9v/(chi+3-i3 c_k v chi)^2`, `c_k=k^3/(6pi)` (D:17–47; actual implementation D:309–329). Its `X=diag(repeat(a(chi),3))`, not diag chi. Exact finite DeltaX must use `a(chi1)-a(chi0)`, never `da(chi0) DeltaChi` except as a tangent approximation.
- A17 current is normalized dipole `c=p/sqrt(v)` with array `(P,3N)`; exciting physical field is `incident+Goff p`, and `state.field` means **receiver field** `(P,m)` (D:925–943). This naming differs critically from A30 `state.field` which is internal exciting field.
- The A17 tangent injection is `repeat(da*dchi,3)*exciting/sqrt(v)`, solved by existing L factorization, observed through GS (D:1008–1044). O:294–311 uses the same injection in B, includes `/scale` in S, and sums over illumination/components in B adjoint. D:27–31 defines L with physical Goff; using normalized c changes RHS by `1/sqrt(v)`.
- A17 real data packing is per illumination `[Re m; Im m]` (E/vendor/Gaussian/A14/src/pairtsom_adapter.py:66–88). Its physical real material coordinates have volume metric `v I` (D:21). M:32–61 restores exact saved Q, verifies same tangent space and basis identity, and M:70–76 rejects mismatched ell/step coordinates. F/code/closure_state_v2.py:63–65 passes `allow_legacy=False`. Reusing a tangent span does not make coefficients portable across a newly generated Q.

## A26 N3 inference boundary

Let `A=Gaussian/A26/a26-om-feedback-prior-completion`.

- A/A26_RESEARCH_REPORT_ZH.md:21–25 describes the same scalar relative current `w=diag(h)(eb+Gb w)`, the L approximation `w=diag(h)eb`, R closure, shared boxed background-adjoint initialization, N3 ordinary BP features, and absence of a decoder-after physics update. These are documented historical definitions, not a new claim of accuracy.
- A/src/a26/encoder.py:410–415 states N1 and N3 share exact L physical base, N3 does not run R or compute beta. Lines 425–439 implement common projected initial image, one L gradient update to `base`, and N3 BP feature from the **pre-box initializer increment**. Lines 451–457 define 12 channels and reshape to image grid.
- A/src/a26/decoder.py:89–97 adds the CNN residual in mass coordinates, scales by `sqrt(volume)`, projects to material bounds, and unpacks to chi. N3 output can be used as an anchor material only with its corresponding admitted paid exact state; the CNN prediction alone is not an exact internal field/current/LU cache.

## Reuse pitfalls

1. A17 polarizability and 3-vector injection cannot be copied into scalar A30.
2. A17 `field` is receiver data; A30 `field` is internal total field.
3. A30 `L_response` versus matrix L are different objects.
4. Material real packing is global allRe/allIm; receiver packing is per source Re/Im.
5. Finite DeltaX and differential dX coincide for A30 raw chi, but not for A17 polarizability or a nonlinear material chart.
6. A30 caches are whitelist-bound and audited; missing cache has no fallback (P:348–403). Existing anchor solve/JVP/receiver transfer operations consume their documented RHS reservations; read-only source reasoning does not authorize physics execution.
