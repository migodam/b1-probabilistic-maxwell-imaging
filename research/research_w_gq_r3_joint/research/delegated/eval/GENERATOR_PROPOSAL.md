# Sol-Eval independent generator proposal (pre-physics)

No physical data have been generated at this stage. Root approval/freeze is required.

- One exposed development world and three independent blind worlds: normal nonuniform wall + circle; cavity wall + ellipse; strong nonlinear target. All target coordinates, radii, real contrast, loss and allowed shape draw use separate independent random draws. Shape strata are explicit conditioning, not a correlated discrete target label library.
- Public model family must disclose support alternatives and all bounds before generation. Proposed rectangular occupied wall x∈[-.75,.75], y∈[.5,1.25], 2×2 macro values, plus one declared cavity-support branch. This is limited macro-support inference, not arbitrary unknown support.
- Proposed continuous target bounds x∈[-.8,.8], y∈[1.8,2.4], major semiaxis [.085,.155], contrast [1.5,6.5] with strong stratum [6.5,8], loss [.005,.20]. Circle or rotated ellipse as allowed by public model.
- Proposed truth wall has independent macro eps/loss draws and predeclared unresolved zero-mean ±.3 texture. If the timed core cannot represent this model error, retain it as explicit failure boundary rather than provide online true means.
- A2 measurements have frequencies [.85,1,1.15], four sources, 16 receivers. Independent N(0,.01²) noise on each real and imaginary channel. Online NPZ contains `y` only; `config.json` contains public bounds, geometry and noise contract only.
- Generator seeds are frozen public integers in protocol; sampled parameters, clean responses and material truth remain in inputs/private. Generator/scorer CLI separation prevents accidental truth use. This is data-flow isolation, not an operating-system inaccessible sandbox.
- Endpoint outputs must be sealed before scoring: endpoint NPZ contains z, predicted data, target material image, wall image and public pixel axes; JSON carries objective, stop reason, all branch summaries and cost. Scoring reads private truth only after ONLINE_SEALED.json.
- Report all branch errors for diagnosis after seal; selected branch must minimize the predeclared measurement objective. No GT-selected endpoint.
- One independent finer-mesh full-wave evaluation is reserved after selected blind result is fixed. Report prediction discrepancy and noise-unit norm; no automatic continuum qualification.
- All physical jobs acquire R3/.fem.lock; one BLAS thread; hard 4GiB RSS and common deadline. Failures, rejected candidates and no-run items remain visible.
