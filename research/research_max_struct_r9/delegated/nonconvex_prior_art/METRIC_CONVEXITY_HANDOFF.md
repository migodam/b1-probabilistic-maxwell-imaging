# van Leeuwen–Yang v1 §2.7: mathematical handoff

Read source: https://arxiv.org/html/2403.15292v1#S2.SS7
Canonical later journal identity: An analysis of constraint-relaxation in PDE-based inverse problems, DOI10.1088/1361-6420/adab89 (2025). This handoff reads the differently titled v1, not the final journal text.

## Original mathematical annotations (not a repaired implementation)

$$\{\mathcal{P}_{i}\}_{i=1}^{n}$$

$$\mathsf{U}=\mathsf{P}_{n}=\text{span}\{p_{i}\}_{i=1}^{n}$$

$$\mathsf{U}=\text{span}\{p_{i}\}_{i=1}^{n}$$

$$J_{0}(c)={\textstyle\frac{1}{2}}\text{trace}\left[\left(M-A(c)M^{-1}D\right)^{*}M^{-1}\left(M-A(c)M^{-1}D\right)\right],$$

$$M_{ij}=\langle p_{i},p_{j}\rangle_{\mathsf{U}}$$

$$A(c)_{ij}=\overline{\mathcal{A}_{c}(p_{j},p_{i})}$$

$$u_{j}=\sum_{k=1}^{n}U_{kj}\,p_{k},$$

$$\sum_{k=1}^{n}\mathcal{A}_{c}(p_{k},p_{i})U_{kj}=\langle p_{j},p_{i}\rangle_{\mathsf{U}}\quad\text{for $i=1,2,\ldots,n$}.$$

$$\mathcal{A}_{c}(p_{k},p_{i})$$

$$\sum_{k=1}^{n}{W_{kj}}p_{k}$$

$$\sum_{k=1}^{n}\overline{\mathcal{A}_{c}(p_{i},p_{k})}W_{kj}=\langle p_{j},p_{i}\rangle_{\mathsf{U}}\quad\text{for $i=1,2,\ldots,n$}.$$

$$e_{ij}=d_{ij}-\left\langle p_{i},\sum_{k=1}^{n}U_{kj}p_{k}\right\rangle_{\mathsf{U}}=d_{ij}-\sum_{k=1}^{n}\langle p_{i},p_{k}\rangle_{\mathsf{U}}\overline{U_{kj}},$$

$$g_{ij}=\left\langle\sum_{k=1}^{n}W_{ki}p_{k}\,,\,\sum_{\ell=1}^{n}W_{\ell j}p_{\ell}\right\rangle_{\mathsf{U}}=\sum_{k,\ell}W_{ki}\overline{W_{\ell j}}\langle p_{k},p_{\ell}\rangle_{\mathsf{U}}.$$

$$m_{ij}=\langle p_{i},p_{j}\rangle_{\mathsf{U}}$$

$$a(c)_{ij}=\overline{\mathcal{A}_{c}(p_{j},p_{i})}$$

$$E(c)=D-M\overline{U(c)}=D-MA(c)^{-1}M,$$

$$G(c)=\overline{W(c)}^{*}M\overline{W(c)}=MA(c)^{-1}MA(c)^{-*}M\,.$$

$$\langle\cdot,\cdot\rangle_{\mathsf{U}}$$

$$A(c_{1})=A(c_{2})\Rightarrow c_{1}=c_{2}\,,\forall\,c_{1},c_{2}\in\mathsf{C}$$

$$\langle\cdot,\cdot\rangle_{\mathsf{U}}$$

$$\langle p_{i},p_{j}\rangle_{\mathsf{U}}=\mathcal{A}_{c}(p_{i},p_{j})$$

## Applicability notes (reviewer paraphrase)

U is chosen exactly as span of the real, sufficiently regular, linearly independent Riesz representations of the source/measurement operators. The source/measurement pairing gives square D. M is the basis inner-product Gram, not the state/adjoint Gram G(c). Noiseless data and well-posedness assumptions precede this section. If the PDE form A_c is affine in coefficient c and M is c-independent, Eq13 is quadratic in c. Remark2 explicitly retains identifiability and physical-approximation caveats. Remark3 needs invertible D for A(c)=M D^{-1} M.

In general §2.6, the PDE-residual projection space W_n depends on the current c. The preceding reduced metric is (I+rho^{-1}G(c))^{-1}; the special finite solution-space choice is what yields Eq13. No general Maxwell convexity follows.

The 1D Helmholtz Gram construction in §3.3.2 Eq21 needs frequency derivatives d'(k), boundary source responses b_i(k)=u_i(1;k), their derivatives, and c(1). Those are additional acquisition requirements.

v1 transcription issue: §2.4/Eq10 includes the inverse metric, whereas §3.3/Eq23–24 prints I+rho^{-1}G without the inverse. This document preserves the discrepancy; parent must verify intended formula before using it.