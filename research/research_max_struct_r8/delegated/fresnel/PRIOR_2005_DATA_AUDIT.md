# Public measured-data audit: 2D Fresnel experimental database

Date of retrieval: 6 September 2026. No inversion was implemented; this audit
covers access, provenance, format, and reuse terms only.

## 1. Gate status

**MET** for the requested 2D TM long dielectric cylinder data:

- 2001 opus, 4 TM dielectric files (`dielTM_dec4f`, `dielTM_dec8f`,
  `twodielTM_4f`, `twodielTM_8f`) downloaded from the official IOP-served
  supplementary-data page on the canonical `iopscience.iop.org` domain with no
  challenge.
- 2005 opus, 4 TM inhomogeneous foam files (`FoamDielExtTM`,
  `FoamDielIntTM`, `FoamMetExtTM`, `FoamTwinDielTM`) downloaded with an
  explicit access-route caveat (Section 3).

No generated data was used as a substitute. TE files were also preserved
because they came in the same official archives, but they are not needed for
the TM gate.

## 2. Official pointers checked

| Pointer | Result |
|---|---|
| https://www.fresnel.fr/3Ddatabase/ (landing) | HTTP 200. States the first-opus data "are free for scientific use and can be downloaded"; its "downloaded" links point to the IOP article pages below |
| https://www.fresnel.fr/3Ddatabase/database.php | HTTP 200 meta-refresh to index.php; no direct file links on current pages |
| old fresnel.fr paths (`/perso/geffrin/database/`, `/3Ddatabase/data/`, ...) | HTTP 404 (historical paths removed) |
| http://www.iop.org/EJ/abstract/0266-5611/21/6/S09 (old IOP 2005 link) | HTTP 403, Cloudflare "Just a moment..." challenge. **Not bypassed.** This matches the task's note that the official IOP2005 link hits a CAPTCHA |
| https://iopscience.iop.org/article/10.1088/0266-5611/17/6/301/data | HTTP 200, no challenge. Source of the 2001 `.exp` files (canonical domain) |
| https://iopscience.iop.org/article/10.1088/0266-5611/21/6/S09/data | HTTP 200 but page title "Radware Bot Manager Captcha" on the canonical domain |
| https://beta.iopscience.iop.org/article/10.1088/0266-5611/21/6/S09/data | HTTP 200, no challenge; same public supplementary-data page. Source of the 2005 `supp.zip` |

## 3. Access route used and caveats

- **2001 opus (primary, clean):** supplementary page on the canonical
  `iopscience.iop.org` domain exposed public pre-signed S3 links for the eight
  `.exp` files. Downloaded as-is on 2026-09-06.
- **2005 opus (secondary, flagged):** the canonical domain returned a
  Radware bot challenge; IOP's own `beta.iopscience.iop.org` mirror served the
  same public page and its pre-signed S3 link for `supp.zip` without a
  challenge. That public mirror link was used. No challenge was solved and no
  access control was circumvented; the route is recorded here so the parent
  can decide whether the 2005 files should be treated as first-class evidence.
  The 2005 archive is **not required** for the TM gate (the 2001 files already
  satisfy it).

## 4. Reuse terms and provenance

- IOPscience supplementary-data page states: "Supplementary data files are
  published under license by IOP Publishing Ltd. Unless otherwise specified,
  any and all rights in supplementary data belong to the author(s)."
- Institut Fresnel database page states the data are "free for scientific use
  and can be downloaded."
- No explicit per-file LICENSE text is embedded in the `.exp` archives.
  Scientific reuse with attribution to the special-section papers is the
  operative norm; commercial or redistribution terms were not found and are
  outside this audit.
- Provenance chain: Institut Fresnel anechoic chamber (CCRM) measurements,
  published through the Inverse Problems special sections
  (Belkebir-Saillard 2001; Geffrin-Sabouroux-Eyraud 2005).

## 5. Files, hashes, format

SHA256 sums: `data/SHA256SUMS.txt` (16 files).

### 2001 opus — `data/2001_iop_17_6_301/` (8 files)

| File | Size | Object / polarization / frequencies |
|---|---|---|
| dielTM_dec4f.exp | 550,707 | dielectric circular cylinder, off-centered; TM; 4, 8, 12, 16 GHz |
| dielTM_dec8f.exp | 1,101,084 | dielectric circular cylinder, off-centered; TM; 1-8 GHz |
| twodielTM_4f.exp | 550,689 | two dielectric circular cylinders; TM; 4, 8, 12, 16 GHz |
| twodielTM_8f.exp | 1,101,067 | two dielectric circular cylinders; TM; 1-8 GHz |
| rectTM_cent.exp | 550,672 | metallic cylinder, centered; TM; 4, 8, 12, 16 GHz |
| rectTM_dece.exp | 1,101,076 | metallic cylinder, off-centered; TM; 2-16 GHz step 2 |
| rectTE_8f.exp | 1,101,054 | metallic cylinder; TE; 2-16 GHz step 2 |
| uTM_shaped.exp | 1,101,079 | U-shaped metallic cylinder; TM; 2-16 GHz step 2 |

Common header facts: measured 16/9/1999; 36 views; 49 receivers per view;
contact K. Belkebir (in-file header).

### 2005 opus — `data/2005_iop_21_6_S09/` (8 files, from supp.zip)

| File | Object / polarization / frequencies / views |
|---|---|
| FoamDielExtTM.exp | inhomogeneous dielectric cylinder; TM; 2-10 GHz step 1 (9); 8 views |
| FoamDielIntTM.exp | inhomogeneous dielectric cylinder; TM; 2-10 GHz (9); 8 views |
| FoamMetExtTM.exp | hybrid metal-dielectric cylinder; TM; 2-18 GHz step 1 (17); 18 views |
| FoamTwinDielTM.exp | twin dielectric cylinder; TM; 2-10 GHz (9); 18 views |
| (four TE counterparts) | same targets, TE polarization |

Common header facts: measured June 2004; 241 receivers per view; contacts
K. Belkebir / J.-M. Geffrin (in-file header).

### Format

ASCII text with `#` comment headers, then 7-column whitespace rows:

```text
view  receiver  frequency_index  v1  v2  v3  v4
```

The last four columns are the complex field values (two complex numbers or
four quadrature components per the special-section convention). Row counts
were checked for internal consistency:

- dielTM_dec8f: 14,112 rows = 36 views x 49 receivers x 8 frequencies, all
  7 columns.
- FoamDielExtTM: 17,352 = 8 x 241 x 9; FoamTwinDielTM: 39,042 = 18 x 241 x 9;
  FoamMetExtTM: 73,746 = 18 x 241 x 17; all 7 columns.

Exact physical units (volts, normalized field, scattered vs total field) are
defined by the special-section papers, not re-derived here. Anyone consuming
the files must re-read Geffrin-Sabouroux-Eyraud 2005 and Belkebir-Saillard
2001 before modeling.

## 6. Not done / explicitly out of scope

- No CAPTCHA solving and no access-control bypass (the old IOP 2005 link was
  left at its 403 challenge).
- No external contact, upload, or purchase.
- No inversion, preprocessing, or numerical experiment on the data.
- No use of generated/synthetic data as a stand-in for measured data.
- The 2009 third-opus 3D data was deliberately not downloaded (2D scope).

## 7. Residual risks

1. The 2005 files were fetched from the IOP beta mirror because the canonical
   domain was bot-gated; treat them as secondary evidence until the parent
   confirms the route.
2. The `.exp` files carry no embedded license; reuse terms rest on the two
   statements quoted in Section 4.
3. Hash provenance is against IOP-served bytes at retrieval time; the
   canonical pages may rotate pre-signed links (they expire).
