# E0/E1 Full66 decoder-only audit

This artifact reports decoder arithmetic and geometry only. Scientific interpretation and final theory decisions remain with the root Codex reviewer.

Independent frozen parents: 256; seed: 2026102807; deterministic edge controls: 15. New Maxwell RHS / states / LU / NN training = 0 / 0 / 0 / 0.

The original P1 MaterialPrior is imported from `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r4-conditional-transport/src/p1_material/prior.py` through sibling R4 `p3r4.core`. No Book or Target is constructed. All 576 material cells, including 256 ROI cells, enter [Re Dg; Im Dg]/24. x=(z,xi) is Gaussian66; native properties are ndtr(xi), and both property derivatives include the standard-normal density.

Observed rank66 counts are 88/256, 163/256, 186/256 and 191/256 at relative cutoffs 1e-6,1e-8,1e-10,1e-12 respectively. The machine-scale transition-A calculation has full column rank64 in 191/256 parents; the remaining 65 parents have plateau-inactive spatial directions. Four of those also have numerically witnessed curved radial same-material paths. These outcomes do not support imposing a universal rescaling quotient or a fixed large Full66 gauge dimension across this prior.

## Observed ranks

{
  "1e-06": {
    "11": 1,
    "15": 1,
    "16": 1,
    "17": 1,
    "20": 1,
    "22": 2,
    "25": 1,
    "28": 1,
    "32": 1,
    "34": 1,
    "37": 1,
    "38": 1,
    "39": 1,
    "41": 3,
    "42": 1,
    "44": 2,
    "46": 1,
    "48": 2,
    "49": 3,
    "50": 3,
    "51": 6,
    "52": 2,
    "53": 3,
    "54": 3,
    "55": 3,
    "56": 3,
    "57": 2,
    "58": 6,
    "59": 7,
    "60": 8,
    "61": 13,
    "62": 11,
    "63": 11,
    "64": 23,
    "65": 38,
    "66": 88
  },
  "1e-08": {
    "11": 1,
    "15": 1,
    "16": 1,
    "17": 1,
    "20": 1,
    "22": 2,
    "26": 1,
    "28": 1,
    "33": 1,
    "34": 1,
    "38": 1,
    "39": 2,
    "41": 1,
    "42": 2,
    "43": 1,
    "44": 1,
    "45": 1,
    "47": 1,
    "48": 1,
    "49": 2,
    "51": 1,
    "52": 4,
    "53": 2,
    "54": 4,
    "55": 3,
    "56": 3,
    "57": 2,
    "58": 4,
    "59": 3,
    "60": 1,
    "61": 3,
    "62": 5,
    "63": 5,
    "64": 9,
    "65": 20,
    "66": 163
  },
  "1e-10": {
    "11": 1,
    "15": 1,
    "16": 1,
    "17": 1,
    "20": 1,
    "22": 2,
    "26": 1,
    "28": 1,
    "33": 1,
    "34": 1,
    "38": 1,
    "39": 2,
    "41": 1,
    "42": 2,
    "43": 1,
    "44": 1,
    "45": 1,
    "48": 1,
    "49": 3,
    "52": 4,
    "53": 2,
    "54": 3,
    "55": 4,
    "56": 3,
    "57": 3,
    "58": 2,
    "59": 3,
    "60": 2,
    "61": 1,
    "62": 3,
    "63": 4,
    "64": 6,
    "65": 6,
    "66": 186
  },
  "1e-12": {
    "11": 1,
    "15": 1,
    "16": 1,
    "17": 1,
    "20": 1,
    "22": 2,
    "26": 1,
    "28": 1,
    "33": 1,
    "34": 1,
    "38": 1,
    "39": 2,
    "41": 1,
    "42": 2,
    "43": 1,
    "44": 1,
    "45": 1,
    "48": 1,
    "49": 3,
    "52": 4,
    "53": 2,
    "54": 3,
    "55": 4,
    "56": 2,
    "57": 3,
    "58": 3,
    "59": 2,
    "60": 3,
    "61": 1,
    "62": 2,
    "63": 5,
    "64": 3,
    "65": 4,
    "66": 191
  }
}

Transition-cell range: [10, 164]; pure-cell range: [0, 242]. Saturation-inactive subspace dimensions: {'0': 191, '1': 4, '2': 3, '3': 5, '4': 2, '5': 1, '6': 3, '7': 2, '8': 3, '9': 3, '10': 2, '11': 4, '12': 3, '13': 2, '14': 4, '17': 3, '18': 1, '21': 1, '22': 1, '23': 1, '24': 3, '27': 2, '28': 1, '32': 1, '33': 1, '38': 1, '40': 1, '44': 2, '45': 1, '49': 1, '50': 2, '54': 1}. Rank-formula mismatches at relative 1e-10: 1.

For B=D_z fraction, native interior properties give rank(Dg)=rank(B)+1+1{fraction not in range(B)} when the fraction is nonzero; the all-host case has rank zero. The property chain is invertible at finite xi, but near-tail conditioning changes numerical ranks. Formula checks use the same 1e-10 relative cutoff; near-seam/tail discrepancies are numerical qualifications, not new exact rank claims.

## Exact inactive directions versus numerical null directions

A transition cell has -temperature < f_i < temperature. The explicitly constructed subspace consists of dz in ker(A_transition), with both property coordinates fixed. Every transition field is then constant along the straight path; all other ROI cells remain exactly on their existing host/pure plateau until a plateau crossing. This supplies a local exact same-material subspace at strict strata, independently of a derivative cutoff. At seams, the crossing interval can collapse to zero. Floating-point SVD kernels are recorded separately and are not automatically called exact gauges.

Finite straight paths at t=0.001, 0.01, 0.1, 0.5 are saved for explicit inactive, SVD-null, radial-property, material-unit strongest, and material-unit weakest nonzero directions. Plateau changes are retained. Material-unit weak directions can require very large latent steps and leave the local stratum; this is explicitly recorded rather than treated as a linear prediction failure of Maxwell physics.

The inactive subspace is computed with a machine-epsilon rank threshold on A_transition, independently of the four material-Jacobian rank cutoffs. Its transition A-residual and a symmetric plateau-crossing interval for its first unit direction are saved. These are floating-point basis computations and finite material-equality witnesses, rather than symbolic certificates for every dimension.

## Global rescaling audit

The tested radial property tangent requests d(delta)=c delta and compensates using B dz=-c fraction. A pure-phase cell has fraction=1 and B-row=0, so it directly obstructs nonzero radial rescaling at that point. With no pure cells, representability of fraction in range(B) must still be checked; uniform-field reasoning alone does not establish it. The full material tangent residual after least-squares compensation is recorded for every parent.

A separate finite candidate alpha*delta with the spatial field divided by alpha is tested at alpha=0.8, 0.95, 1.05, 1.2 wherever native property intervals permit. The quintic fraction law is not homogeneous. These empirical finite residuals test that specific proposed global field/property rescaling; they do not prove the absence of every nonlinear same-material path. Saturation-inactive fibers and low-transition nonlinear amplitude fibers are distinct objects.

A stronger nonlinear test divides each transition fraction by alpha, inverts the original quintic, and fits the resulting fields with the original A. No-pure parents 092, 201, 233 and 251 (19,40,10,14 transition cells) admit corrected finite paths at all four prescribed tangent steps with material RMS residuals below 2.5e-16 and no host crossing. Their uncorrected straight radial paths have nonzero curvature residuals, reaching respectively about 1.77e-5, 8.87e-8, 3.97e-4 and 5.96e-6 at t=0.5. This is finite numerical evidence for particular curved same-material fibers in the original prior, beyond the straight inactive subspaces. The other no-pure parent 136 has 88 transition cells, numerical rank66 at 1e-10 and 1e-12, and a nonzero nonlinear field-fit/material residual; it is not promoted to an exact gauge. The uniform transition control also admits the nonlinear correction but is kept separate from the independent-parent results. Raw corrected x/material vectors and all residuals are saved.

## Numerical and cost qualifications

Three full 66-column central differences were evaluated at steps [0.0001, 1e-05, 1e-06] for every prior and edge point. Max prior relative Frobenius discrepancy: 1.52992e-08. Absolute singular ranks at 1e-8,1e-10,1e-12,1e-14, operator-norm FD errors, and FD-informed rank intervals are retained. Max projector idempotence discrepancy: 1.40822e-14.

Deterministic coordinate/material/J replay: True / True / True. Final decoder audit wall 7.180s; CPU 7.060s; peak RSS 490471424 bytes. Cost excludes plotting/report writing after the E0 timer receipt. Earlier implementation sweeps are recorded separately in the worker execution note; all sweeps consumed zero Maxwell RHS.

Raw complex128 material and Jacobian arrays plus float64 coordinates, fractions, fields, full singular spectra and FD errors: `raw/decoder/FULL66_DECODER_ARRAYS.npz`. Raw vectors and finite comparisons: `raw/decoder/DIRECTIONS_AND_PROJECTORS.json` and `raw/decoder/FINITE_STEP_COMPARISONS.json`. No old history or ledger is modified.
