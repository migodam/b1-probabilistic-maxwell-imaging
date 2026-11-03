# Fixed finite-eta port-profile comparison

Both arms are sealed development predictions. No GT was opened by these programs. Lambda0 is the known state-profile/WRI control; lambda1 is the session B port augmentation, with no originality or guaranteed-recovery claim.

Shared eta=0.0042844779316861901, computed once as trace(Cref Cref*)/48 at the public host0.1. All data/state/port terms divide observed||Y||² and preserve the same historical prior. Both arms start from the host, run250L-BFGS iterations, select evaluated materials only by original physicalPhi including the host, then call the same original scaled corrector once for1000iterations.

|Arm|Aux endpoint objective|Selected physical residual before corrector|Final original residual|Final original Phi|Final PGinf|LU|Complex RHS|Wall(s)|
|---|---:|---:|---:|---:|---:|---:|---:|---:|
|lambda0|0.0129772451|0.485690525|0.3781299153|0.07150441947|1.55124e-06|1290|46464|18.5911|
|lambda1|0.00852455093|0.168120228|0.09325484465|0.004368925422|6.65299e-07|1297|47124|14.7917|

All four optimizations (two auxiliary, two correctors) reached their iteration limit. The final projected gradients exceed1e-8; no first-order convergence or basin certificate is claimed. Lower original residual is not a material-recovery conclusion.

The one new4×4 lambda1 directional wiring check had absolute error5.621e-13 and relative error2.213e-10; it used only a public constant0.55 observation. Its5LU/106RHS are separate. The shared576-cell host reference costs1LU/60RHS. Combined total: 2593LU/93754complexRHS; numerical process timer sum33.429461s wall/32.101670s CPU. No failed run, eta sweep, lambda sweep, or derivative-check repetition occurred.

For each profile gradient the original Maxwell LU is retained for48receiver-transpose,12profile-current and12adjoint RHS. Both arms use the same derivative path, including a charged adjoint whose lambda0 term vanishes. The shared reference C is reused at the initial host. Small48×48 Gram/Cholesky and12small RHS are counted separately from physical solves.

The two result directories preserve every evaluated material, auxiliary/physical objective histories, the auxiliary optimizer endpoint, the original-Phi-selected seed, the common corrector trace and sealed prediction. Material scoring is delegated separately; online code has no truth access.
