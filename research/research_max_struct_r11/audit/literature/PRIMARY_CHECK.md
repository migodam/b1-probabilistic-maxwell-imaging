# Primary-text check and extraction correction

Both primary PDFs were obtained through Undermind and read locally with pdftotext. Temporary signed download URLs and full copyrighted PDFs are not repository artifacts. Stable bibliography is in verified_metadata.json.

**Budko & Samokhin (2005):** the automatic full-text extraction got the sign of Eq. (5.9) wrong. Direct inspection of page9 gives `F(theta)=I/3−theta thetaᵀ`, not `theta thetaᵀ−I/3`. The original extraction JSON is retained as a failed extraction record; this correction governs the research text. Eq. (5.1) has the local `I+chi/3` term and subtracts the Green integrals. Eq. (5.10) is `Phi=I+P_L chi`. These signs agree with the nondecaying longitudinal part used in the R11 background discussion. The assumptions immediately following Eq. (5.7) include Hölder-continuous constitutive tensors on the whole space, across the boundary. Do not apply that spectrum verbatim to sharp-interface voxel DDA.

**Costabel, Darrigrand & Sakly (2012):** direct primary text confirms Theorem3.2's gradient/divergence-free/harmonic decomposition, Corollary3.3's essential electric spectrum `{0,1/2,1}` for smooth interfaces, and Theorem5.1 for the electric/magnetic contrast operator. Its sharp homogeneous-body assumptions differ from Budko's whole-space smooth coefficients. Neither result certifies any finite rank-120 truncation or our imaging endpoints.

Scope: source equation verification, not a numerical rerun, global novelty certificate, or complete literature review.
