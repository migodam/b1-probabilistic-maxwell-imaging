# Linear known-reference physical-existence material chart

`src/wr/native_world.py` introduces `LinearWorld(ModulatedWorld)`. This is an
**ORACLE physical-existence search representation around a known reference**,
not an online inverse algorithm, a new material/noise prior, a prior probability,
or a demonstrated physical witness. No existing running source was changed.

The chart has 9 x-nodes from -1.6 to 1.6 and 6 y-nodes from 0.4 to 1.4, hence
54 bilinear tensor hats and 108 real parameters. Ordering is y-major/x-minor.
Points outside the rectangle clamp to end nodes. In one dimension each row has
at most two adjacent nonnegative weights summing to one; the tensor product is
therefore nonnegative and sums to one everywhere, including outside the node
rectangle. Thus parameters in [-1,1] imply both Bq fields are in [-1,1].

The decoder is exactly affine:

`epsilon = 1 + c0*(1 + 0.6 B q_c) + i*sigma0*(1 + 0.75 B q_s)/f`.

No tanh or chart saturation appears. Since original c0 <= 4 and sigma0 <= 0.036,
wall contrast <= 6.4 and wall conductivity <= 0.063 globally. Both multipliers
remain strictly positive (at least 0.4 and 0.25), so the wall material support
is identical to the supplied reference support, including the inherited,
explicitly declared cavity overrides. No added pixel holes or new features
are introduced by the chart. Original material bounds permit contrast <= 8
and conductivity <= 0.15, so this representation stays within them.

The inherited fixed probe override is retained exactly, including zero
parameter derivatives inside the probe. **A fixed probe of contrast 7 exceeds
the wall-chart bound 6.4 while remaining within the original contrast bound 8.**
Accordingly 6.4 is a global wall-modulation bound, not a whole-world bound when
such a probe is present. Probe/cavity geometry is unchanged. With neither
explicit cavity edits nor a probe, support is precisely the original world.

The analytic derivative is constant in q:

`[0.6*c0*B, i*0.75*sigma0*B/f]`, zeroed at fixed probe points.

`coefficient`, `components`, `derivative`, inherited `jacobian`, and
`material_metadata` follow the existing parent interface. `lower` and `upper`
expose the required [-1,1] box. Out-of-box, complex, nonfinite, or wrong-sized
parameters are rejected instead of silently clipping. Frequency must be
finite and positive. Search drivers must use this box rather than the old
DCT chart's [-2,2] box.

## Mechanism verification

Source AST parsing passed. `MECHANISM_CHECK.json` records one small zero-PDE
mechanism check using 1633 points including all nodes, exterior points, wall
and air points. Results:

- Nonnegative basis; partition error 2.22e-16; nodal identity error zero.
- Directional finite-difference decoder derivative error 1.27e-9 at step 1e-6.
- Affine identity error 8.88e-16.
- Parameter-box extremes preserve support and satisfy wall material bounds.
- The zero-parameter reference equals the parent decoder exactly, including
  explicit cavity and fixed contrast-7 probe overrides; their derivatives zero.
- Maxwell/PDE calls, new RHS and new factorizations: all zero.

The NumPy runtime lacks `skfem`, so importing the parent module normally was
blocked. Verification executed the actual class/function ASTs in an isolated
NumPy namespace, excluding backend imports; it did not replace the numerical
methods with stubs. Full backend import/integration and physical search remain
NOT_RUN. No package was installed and no physical ledger entry was created.
The finite sample check supports the decoder mechanism; the partition argument
above supplies its global algebraic bound. Neither proves PDE or witness status.
