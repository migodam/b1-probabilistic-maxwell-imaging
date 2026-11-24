# Two-sided port test for a supported material class

Let \(H\) be the declared active support and suppose the material is zero outside \(H\). The finite forward model gives \(J_H=M_H(B_H+G_{HH}J_H)\), \(F(M)=S_HJ_H\), and \(0\preceq M_H\preceq uI\), with a common real diagonal \(M_H\) for all illuminations. Define

\[
A_H=\tfrac1u I-\operatorname{Re}_H G_{HH}+t\operatorname{Im}_H G_{HH},\qquad z=(1+it)/2.
\]

Assume \(A_H\succeq\gamma I\) with \(\gamma>0\). For any source coefficient vector \(a\), put \(j=J_Ha\), \(b=B_Ha\) and \(e=b+G_{HH}j\). Componentwise \(j_i=m_ie_i\) implies \(u^{-1}\|j\|^2-\operatorname{Re}(j^*e)\le0\) and \(\operatorname{Im}(j^*e)=0\). A real multiple \(t\) of the optical equality gives \(j^* A_H j\le2\operatorname{Re}(j^*zb)\). Cauchy–Schwarz yields

\[
\|J_Ha\|_2\le C_H\|B_Ha\|_2,\qquad C_H=\frac{\sqrt{1+t^2}}{\gamma}.
\]

Equivalently \(J_H^*J_H\preceq C_H^2 B_H^*B_H\). This gives \(\|J_HW\|_F\le C_H\|B_HW\|_F\) for any matrix \(W\) of coherent illuminations. The material remains fixed when source columns are mixed. The same inequalities hold for passive complex \(m_i\) with \(\operatorname{Im}m_i\ge0\), \(|m_i|^2\le u\operatorname{Re}m_i\), and \(t\ge0\), in the stated convention \(\operatorname{Im}(j^*e)\le0\).

For source matrix \(W\) and receiver matrix \(V\), both nonzero, every admissible \(M\) obeys

\[
\|Y-F(M)\|_F\ \ge\ 
\frac{\big(\|V^*YW\|_F-\|V^*S_H\|_F C_H\|B_HW\|_F\big)_+}
{\|V\|_2\|W\|_2}.
\]

Indeed \(\|V^*S_HJ_HW\|_F\le\|V^*S_H\|_F\|J_HW\|_F\) and \(\|V^*(Y-F(M))W\|_F\le\|V\|_2\|W\|_2\|Y-F(M)\|_F\). The corresponding operator norm bound for \(V^*S_H\) can tighten the estimate. The reference implementation uses its Frobenius norm. Binary64 eigenspaces propose prefixes of the weakest receiver directions of \(S_HS_H^*\) and weakest incident directions of \(B_H^*B_H\). It checks the chosen rank pair with outward arithmetic, including both basis norms. The inequality permits these choices to depend on \(Y\). Setting one matrix to the identity gives the matched one-sided controls. [Implementation](fast_support.py) · [Ablation replay](ablation.py)

### Exact tile reuse

For a shared finite convolution kernel, \(G_{ij}\) depends on the difference of native grid indices. All translated copies of the same lattice tile have identical principal Green blocks under the same coordinate ordering. A tile may meet a physical domain edge when it still consists of exactly the same in-domain indices, up to translation. We use identical 24×24 tiles of the native 96² finite convolution. One outward Cholesky residual check certifies \(A_H\succeq\gamma I\) for all 16 locations, all data vectors and both apertures of that kernel. This reuse requires unchanged lattice, frequency, kernel and material cap. Receiver and incident matrices are cut at each tile, so their directional norms change with location and aperture. The paid shared work is one exact inverse DFT for the stored Green symbol and one 576-dimensional Cholesky with an interval residual certificate. Each subsequent pair of matrices \(V,W\) costs two small numerical eigendecompositions and verified products. [Cold and repeated replay](fast_development.json)

The exact active-state fallback uses \(K=S_H^*S_H+\alpha A_H\) at a different fixed \(t,\alpha\), with a verified positive gap, verified negative subspace and source eigenvalue tail. The soft-background route uses a variable cap \(U_i>0\), weights \(\mu_i=\alpha/U_i\), \(\lambda_i=\alpha t\), and one-sided Fourier compression of the native finite Green operator. Its diagonal is \(1/U_i+c\) and the retained terms are those with \(-\operatorname{Re}g+t\operatorname{Im}g<c\). Zero rows are eliminated exactly. A high count of negative directions in the lower envelope can exhaust the source tail. In that event the case stays unresolved. [Exact active fallback](restricted.py) · [Positive-background proof](native_box.py)

### Perturbed finite operators and measurement gain

For actual Green operator \(G+\Delta G\), the change in \(A_H\) has norm at most \(\sqrt{1+t^2}\|\Delta G_{HH}\|_2\). Thus the same result uses \(\gamma_*=\gamma-\sqrt{1+t^2}\|\Delta G_{HH}\|_2>0\) and \(C_*=\sqrt{1+t^2}/\gamma_*\). Native \(S+\Delta S\) and \(B+\Delta B\) can be covered by increasing the two weak factors respectively by \(\|V^*\Delta S_H\|_F\) and \(\|\Delta B_HW\|_F\). The Green perturbation bound with \(t=32\), \(\gamma=.009\) requires \(\|\Delta G_{HH}\|_2<.009/\sqrt{1025}\), approximately \(2.81\times10^{-4}\), when only the global norm bound is available. Measured data from a different physical model need their own verified tolerances.

Alternatively use a single additive model-error bound \(b\) over the whole class. Suppose \(Y=\mathcal G(F(M)+\Delta)+N\), with \(\|\Delta\|_F\le b\), \(\|N\|_F\le\delta\) and \(\|\mathcal G-I\|_{F\to F}\le\eta<1\). Then \(\|F(M)+\Delta\|_F\le(\|Y\|_F+\delta)/(1-\eta)\) and

\[
\|Y-F(M)\|_F\le b+\frac{\delta+\eta\|Y\|_F}{1-\eta}.
\]

Use this radius, or a tighter calibrated bound, for strict exclusion. Uncertain or instrument-specific port masks need actual common receivers within each source group. Disjoint observation groups allow squared bounds to add. For the present chamber studies, calibration and model perturbation scales give a sensitivity analysis. [Existing chamber error contract](../../z_alg_optical_inertia/fresnel/REPORT_ZH.md)

### Attribution and scope

Energy and power bounds on electromagnetic response are established tools. The conservation-law step follows that tradition and the finite-model forms already proved in Paper 1. Orthogonal weak transmitting and receiving directions also have an established singular-mode interpretation. The implemented combination tests a specific supported material class directly against finite noisy multi-source data, using translation to amortize one native coercivity proof and retaining the measured weak-port signal. It establishes sufficient exclusion for a declared class. Support recovery, completeness of the test, and solver speedups across broader classes require separate evidence. [Molesky et al., physical power bounds](https://www.princeton.edu/~arod/papers/Molesky20_PRR.pdf) · [Miller, orthogonal communication modes](https://doi.org/10.1364/AO.39.001681) · [Paper 1 finite-model forms](../../z_alg_optical_inertia/ip_submission_v2/manuscript_ip.pdf)

### Exposed follow-up using the same energy inequality

For full receive/source bases, let \(a_i=\|S_H^*v_i\|_2\), \(b_j=\|B_Hw_j\|_2\). The energy inequality gives \(|v_i^*F(M)w_j|\le C_Ha_i b_j\) for every pair. Triangle inequalities on the complex entries, followed by summing their squares, give
\[
 d(Y)\ \ge\ \frac{\sqrt{\sum_{ij} (|v_i^*Yw_j|-C_H a_i b_j)_+^2}}
 {\|V\|_2\|W\|_2}.
\]
For orthonormal bases this includes every rectangular lower bound obtained with the same basis ordering, because \(\|(x-y)_+\|_2\ge\|x\|_2-\|y\|_2\) for nonnegative real vectors. In the finite implementation, all row norms and basis scaling are verified. Norms of \(S_H^*v_i\) and \(B_Hw_j\) can be prepared once per known location and aperture. [Exposed comparison](entrywise.py)

The coherent energy inequality also gives a transfer representation \(J_H=\mathcal T B_H\) with \(\|\mathcal T\|_2\le C_H\). In finite dimensions define \(\mathcal T(B_Ha)=J_Ha\), which is well defined by the energy inequality, then extend it by zero on the orthogonal complement of the range. For a positive scalar \(\rho\), use the Hermitian block
\[
\mathcal K_\rho(Y)=\begin{pmatrix}
  \rho C_H S_HS_H^*&Y\\ Y^*&C_H B_H^*B_H/\rho
\end{pmatrix}.
\]
At exactly supported \(F=S_H\mathcal T B_H\) the block is positive semidefinite by congruence from \(\begin{psmallmatrix}I&\mathcal T/C_H\\\mathcal T^*/C_H&I\end{psmallmatrix}\). The off-diagonal perturbation has Frobenius norm \(\sqrt2\|Y-F\|_F\), giving a lower bound \(\|\mathcal K_\rho(Y)_-\|_F/\sqrt2\) after outward verification of the negative eigenvalues. This is a distance to the positive semidefinite cone. The two exposed checks use the energy inequality and the known port factors. The finite Maxwell class also imposes the shared constitutive law. [Douglas, contraction factorization](https://doi.org/10.1090/S0002-9939-1966-0203464-1) · [Higham, Frobenius PSD nearness](https://www.sciencedirect.com/science/article/pii/0024379588902236)

On the 256 previously exposed site checks, the entrywise bound increased the strict count from 191 to 196. The block test added one more. The remaining 59 state-level checks excluded 43 wrong sites and retained all 16 generating sites. The two-aperture cache, entrywise tests, block tests and active-state calls used 5.307, .149, .571 and 132.538 process CPU seconds, respectively. An initial stopped pass spent 19.07 seconds before correcting the binary64 export of an upper endpoint. [Full replay](entrywise_exposed.json) · [Stopped pass](entrywise_stopped_rounding.txt)

### What the bounded transfer relaxation retains

For a reciprocal Green model with \(G_{HH}^T=G_{HH}\), the supported state transfer of a diagonal material has the identity
\[
\mathcal T_M=(I-M_HG_{HH})^{-1}M_H=M_H(I-G_{HH}M_H)^{-1}.
\]
It follows that \(\mathcal T_M^T=\mathcal T_M\). The same energy proof bounds its norm, since it applies to every incident vector on the tile. The positive gap gives uniqueness of the state at zero incident field, which in the finite model also gives existence for arbitrary incident vectors. A reciprocal-control witness therefore fits \(Y\) with \(S_H\mathcal T B_H\), enforces exact ordinary-transpose symmetry of \(\mathcal T\), and checks its norm and data discrepancy with intervals. The symmetric witnesses use \(\mathcal T=C_HQXQ^T\) with stored entrywise-symmetric \(X\). This algebraic control is scoped to the class of bounded symmetric maps.

On the 49 wrong sites first excluded by the finite-state fallback, arbitrary bounded transfer fits were verified for 39, and bounded symmetric transfer fits for 38. All 16 correct sites had fits in both relaxations. The arbitrary and symmetric searches cost 1.819 and 3.733 seconds of process CPU. The unresolved 10 and 11 wrong sites include finite proposal failures, so these counts certify lower bounds on the number of relaxed fits. The 38 strict symmetric fits demonstrate that finite-model exclusion in those cases uses constraints beyond a scalar current gain and reciprocity. [Arbitrary-transfer witnesses](contraction_control.json) · [Reciprocal witnesses and failures](reciprocal_control.json)
