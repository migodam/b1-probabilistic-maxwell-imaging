# Frozen A30 backend dependency addendum

A30_FROZEN_SOURCE.zip intentionally contains src/a30 and root license/notice files. Its physics adapter imports the unchanged vendored A29 backend through vendor/frozen_a29/src (a24.physics/contracts and related canonical packages). Therefore source.zip alone is insufficient for a new physics run.

Also extract A30_FROZEN_VENDOR.zip into the same NEW directory. Both archives use prefix A30_FROZEN/ and the same frozen A30 commit717d3d26a290248bd2b9aeb55f3d03694b407b84. Retain embedded source/provenance/licenses. Then B1_P3L_A30_ROOT points to A30_FROZEN. Use a fresh Python process to avoid mixing other versions of a20/a24/a25/a27.

This corrects the earlier incomplete dependency export note; original source and first receipt are preserved. No solver implementation was replaced, no arrays interpreted, no geometry/LU/RHS/physics/training was run. Saved-material network replay does not import A30 and is unaffected. Any new physical campaign still requires a separate output/root and authorization; do not reopen the sealed R2 Book.
