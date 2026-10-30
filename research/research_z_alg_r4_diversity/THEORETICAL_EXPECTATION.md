# Why lower contrast does not automatically favor GlobalDual

This is a mathematical expectation declared before new truth is disclosed. It is not a result of the 72-observation experiment. The numerical results must decide whether there is a useful advantage.

## 1. What changes when materials are weak

For fixed G, S and incident fields E0, let m=epsilon a, with independent real material coordinates and a bounded fixed a. In a region where the Neumann expansion converges,

$$J=\epsilon D_aE_0+\epsilon^2D_aGD_aE_0+O(\epsilon^3),$$

$$Q=E_0+\epsilon GD_aE_0+O(\epsilon^2),\qquad H=Q^*Q=H_0+O(\epsilon),\quad H_0=E_0^*E_0.$$

Thus multiple scattering weakens continuously. It is not removed from any executed forward or inverse calculation. The generator records the Born/full response difference from its already-paid full forward, solely as an offline regime descriptor.

If H0 is positive definite and the field perturbation is small relative to its smallest eigenvalue, the inverse Gram is close to a fixed matrix. Explicitly, if delta=||H-H0||2 < lambda_min(H0),

$$\|H^{-1}-H_0^{-1}\|_2\le\frac{\delta}{\lambda_{\min}(H_0)[\lambda_{\min}(H_0)-\delta]}.$$

Consequently, with R=F(m)-Y,

$$\left|\operatorname{tr}(RH^{-1}R^*)-\operatorname{tr}(RH_0^{-1}R^*)\right|
\le \|R\|_F^2\,\|H^{-1}-H_0^{-1}\|_2.$$

Low nominal amplitude alone is not a numerical certificate of this condition: nearly dependent source fields can make the smallest Gram eigenvalue small. The experiment does not assume that every low-amplitude case is a Born problem or that H is exactly constant.

The actual paired experiment holds a nonzero background b fixed while varying stroke amplitude A. Its precise small-A expansion is about that background, not about vacuum. Set Mb=I-Db G, Jb=Mb^-1 Db E0, Qb=E0+GJb and m=b+A a. Where the reference and updated resolvents are invertible,

$$F(b+Aa)-F(b)=A S M_b^{-1}D_aQ_b+A^2 S M_b^{-1}D_aG M_b^{-1}D_aQ_b+O(A^3),$$

$$Q(b+Aa)=Q_b+A G M_b^{-1}D_aQ_b+O(A^2),\qquad H(b+Aa)=Q_b^*Q_b+O(A).$$

The same inverse-Gram bound applies with Hb=Qb*Qb replacing H0. This statement does not grant the solver the hidden true background b; it describes the controlled material family. Strong propagation about b or a nearly singular Hb can invalidate a naive amplitude-only argument.

## 2. A fixed Gram changes weighting, not information

The leading weak-material GlobalDual loss is a fixed weighted least-squares loss, ||R H0^(-1/2)||F squared. An invertible fixed source weighting leaves the exact zero-residual solution set unchanged. It adds no frequencies, receivers, known support or material labels.

It can improve or worsen conditioning of an optimization method. In noisy data it also changes which source combinations are emphasized; that need not be the statistically appropriate noise weighting. Even with identical generic regularizer coefficients, weighted and unweighted penalized minima need not be identical.

The experiment therefore retains the original full-data objective for final online selection and correction. It compares GlobalDual directly with a matched data-adjoint L-BFGS procedure, as well as native and reciprocal mature TRF. Faster convergence than Full-J alone would not isolate the benefit of the Gram objective, because first-order source-adjoint methods also avoid the receiver Full-J work.

## 3. What may reasonably improve

At a zero-residual interior root the derivative of the residual weight contributes no first-order gradient term. Away from a root, the exact GlobalDual gradient includes a material-dependent field/metric term that is quadratic in R before fixed normalization. This term and the source weighting can steer an iteration differently, which was useful for an old R3 target. They do not imply a globally better basin or better material distance.

Without regularization, at an exact root m*, the material Hessian of the raw GlobalDual objective is the weighted data Gauss–Newton form:

$$D^2\psi(m_*)[h,h]=2\operatorname{tr}\{DF(m_*)[h]H(m_*)^{-1}DF(m_*)[h]^*\}.$$

Terms involving the derivative of H^-1 vanish there because R=0. Positive definite H preserves the exact Jacobian nullspace; it can change curvature and conditioning but not create additional local information. Shared regularization, finite residuals and active bounds require their own analysis and do not inherit an unregularized-root claim automatically.

For weak targets ordinary solvers may already reach the available accuracy floor. A tie, or an added Gram cost without material improvement, is then a scientifically meaningful answer. For intermediate contrast or noisy observations there may be a finite-budget advantage, but it must be shown rather than presumed.

## 4. Claims allowed by this experiment

- Accuracy: posthoc RMSE divided by true material increment, with Otsu segmentation based only on the reconstruction. Raw RMSE is retained so the normalization is explicit.
- Search: online top1 uses measurements and the common regularizer. Phase-best candidates are finite optimizer records, not independent solution branches or posterior draws.
- Cost: full LU, source/adjoint/receiver RHS, Gram and measured wall/CPU; hardware results stay separate. The CUDA path accelerates physical kernels while SciPy optimization remains on the i9 host.
- Scope: twelve shapes and one paired noise direction per shape. There are 72 observations but only twelve independent shape constructions; contrast/noise variants are paired, not72 independent objects.

No universal GlobalDual superiority, low-contrast uniqueness, posterior guarantee or hardware speedup is claimed in advance.

## 5. A useful exact source-space property, not additional measurements

For any invertible source mixing B, replace incident fields and measured columns by E0 B and Y B. Maxwell linearity in illumination gives Q'=QB and R'=RB. Therefore

$$\operatorname{tr}\{RB[B^*HB]^{-1}B^*R^*\}=\operatorname{tr}(RH^{-1}R^*).$$

The raw GlobalDual objective is invariant under this invertible source-basis change. This explains why its natural object is an illuminated field subspace rather than twelve arbitrarily scaled source coordinates. Truncating source columns is not invertible and does not inherit this identity.

A fixed source prewhitener derived from the public-reference Gram could improve conditioning without adding information. It was not implemented in this experiment. In addition, the executed objective has a fixed reference normalization and a material regularizer: recomputing that scalar after changing source units can alter the relative regularization. The exact raw identity must not be promoted to an invariance theorem for every implementation or for an unweighted noise model.

All executed GlobalDual derivatives retain the complete material-dependent Gram derivative. An approximation of objective values by a fixed Gram does not, by itself, justify discarding that derivative.

## 6. Small-noise local accuracy has a further limitation

This is the standard linear unbiased least-squares argument, derived here to delimit the imaging claim. It is not a new Maxwell theorem. Suppose an interior true root has a full-column-rank real data Jacobian A, no regularizer is used, and real-packed measurement noise has covariance sigma² I. At first order about that root, any fixed positive definite data weight W gives

$$\widehat h_W=L_W\eta,\qquad L_W=(A^TWA)^{-1}A^TW,\qquad L_W A=I.$$

The unweighted estimator is L0=(A^T A)^-1 A^T. Write D=L_W-L0. Since DA=0, L0 D^T=(A^T A)^-1 A^T D^T=0. Consequently,

$$\operatorname{Cov}(\widehat h_W)-\operatorname{Cov}(\widehat h_0)=\sigma^2 DD^T\succeq0.$$

GlobalDual supplies the real-packed weight induced by H(m*)^-1. Its material-dependent weight derivative contributes only at higher noise order in this root-local expansion. Therefore a source-field Gram weight has no general small-noise material-variance advantage over the correctly noise-weighted ordinary data loss. For nonuniform known noise the corresponding efficient comparator is noise-whitened least squares.

This conclusion concerns the local linearized estimator, not finite-budget nonlinear optimization. Basin changes, active bounds, regularization, discretization errors and incomplete convergence can dominate the experiment. A GlobalDual candidate can still be useful far from the root, but that benefit must appear in the actual paired reconstructions; it cannot be inferred from the existence of a Gram normalization alone.
