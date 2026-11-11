# Prior art audit — verified scope and remaining gaps

No global originality claim is made. Shared-slope total least squares, Schur projection and Douglas–Rachford are existing tools. R12 investigates a particular finite-illumination Maxwell implementation and its measured search behavior.

## Directly inspected primary text

Levinson and Markel, *Solution of the nonlinear inverse scattering problem by T-matrix completion. I. Theory*, Physical Review E94,043317(2016), [accepted manuscript](https://link.aps.org/accepted/10.1103/PhysRevE.94.043317), [DOI](https://doi.org/10.1103/PhysRevE.94.043317).

The abstract and introduction explicitly allow nonlocal interaction potentials at intermediate iterates and seek a data-compatible T-matrix corresponding to an approximately diagonal potential. Therefore allowing nonphysical intermediate scattering states is not a new general principle. Their full T/V representation differs from our actual-illumination J/Q state representation; this is an implementation/representation distinction, not proof of superior search or first invention. Dense background projection still requires quadratic storage in R12.

Evidence tier A for the inspected primary abstract/introduction; full update-by-update comparison remains ongoing. Do not infer convergence from the absence of an explicitly minimized scalar objective.

## Retrieval limitations

Semantic Scholar complementary search on2026-11-11 returned429 for paper and snippet endpoints; raw errors retained. Direct DOI access to PADMM-CIE10.1109/TMTT.2023.3298203 and EM-ADMM10.1017/S175907872000015X failed through the web reader. Their identities and descriptions currently remain inherited Pro references, not newly inspected full-text evidence. Specific overlap with their proximal updates remains unresolved.

## Current permissible contribution wording

An implementation and controlled experimental study of shared real-material field/current projection and reflected feasibility search, including actual bad-endpoint interventions, data-feasible selection limits and cost accounting. Independent24 scoring is now complete and does not support a universal DR replacement for G19. A toy ring escape alone does not establish an algorithm superiority claim.


## Additional primary-source inspection (2026-11-11)

**WRI:** van Leeuwen and Herrmann, *Mitigating local minima in full-waveform inversion by expanding the search space*, GJI195,661–667(2013), [DOI](https://doi.org/10.1093/gji/ggt258); [author laboratory formulation](https://slim.gatech.edu/research/full-waveform-inversion). The inspected original abstract and author formulation minimize data misfit plus a wave-equation penalty in an enlarged wavefield/material space. This directly precedes the broad claim that relaxing physical intermediate states can mitigate local minima. R12 instead uses a hard propagation/data ball and a joint shared-material field/current projection with reflected memory; that specific construction needs comparison, rather than a claim to invent lifted search. No WRI implementation was run here, so no empirical superiority over WRI is established.

**Electromagnetic ADMM:** Jian Liu, Huilin Zhou, Liangbing Chen and Tao Ouyang, IJMWT12,790–796(2020), [publisher abstract](https://www.cambridge.org/core/journals/international-journal-of-microwave-and-wireless-technologies/article/abs/alternating-direction-method-of-multiplier-for-solving-electromagnetic-inverse-scattering-problems/C4CDA44E47B372F5B5B936025DCB2EA0), DOI10.1017/S175907872000015X. The publisher abstract is now directly accessible: augmented-Lagrangian splitting updates contrast sources, material contrast and dual variables, with numerical and experimental examples. This confirms relevant EM splitting prior art. Full equations remain unavailable; the abstract's global-minimization wording is not taken as a proved nonconvex guarantee. R12 does not claim to outperform this unimplemented method.

**PADMM-CIE:** DOI10.1109/TMTT.2023.3298203 remains inaccessible through the primary reader. A secondary discovery record identifies a 3D CIE/Fourier-basis and regularization method, but this is not enough to adjudicate equation-by-equation novelty. Retain it as an unresolved closest comparison, not a fabricated verified baseline.

|Claim|Permissible status|
|---|---|
|First use of nonphysical intermediate states to escape inverse-scattering minima|Rejected: WRI and DCTMC already provide this principle|
|First splitting/dual-state approach for electromagnetic inversion|Rejected: CSI and EM-ADMM precede it|
|New general DR convergence theorem|Not claimed; nonconvex constitutive projection remains|
|Specific exact shared-real-material projection and reusable finite-data-ball implementation|Implemented; mathematical building blocks known; application-level distinction only|
|Controlled rescue of selected actual mature-solver bad endpoints|Empirical in the named exposed R11 cases, with AP and paid-budget controls|
|Independent broad material-recovery advantage|Confirm24 scored: no overallDR advantage overG19; see the paired, bounded-scope results|
|Continuum-accurate high-index3D or real-data superiority|Not established by present finite-DDA/measurement evidence|

The old retrieval failures remain in the audit trail; later access to the EM-ADMM abstract supersedes only that abstract-access limitation, not the missing full-text comparison.


**CIE primary text revisited:** the [2019 article](https://pmc.ncbi.nlm.nih.gov/articles/PMC8320915/) explicitly combines vector CIE-I with FFT-TSOM and relaxes the current/material contraction parameters separately. This confirms that dual-beta and vector-current subspace ideas are inherited. Our additional common physical corrector/prior and small-DDA adaptations are stated as repository comparison contracts, not the original paper's literal experimental implementation.

**Projection precedent:** [Elser's primary preprint abstract](https://arxiv.org/abs/math/0111080) already distinguishes the difference-map fixed point from the recovered object and discusses avoidance of stagnation. R12 does not rename this as newly discovered reflected memory. A new full-text fetch of CSI DOI10.1088/0266-5611/13/6/013 was blocked; verified metadata and inherited implementation lineage remain separate evidence levels.
