# Existing48 primary-view gallery observations

This authorized offline visual review covers all 48 cases on all 12 saved diagnostic-row gallery pages, including visibly distorted and similar-looking cases. It compares A+ with DR and records CIE, Dual-CIE and G19 context. It does not compute new scores, thresholds, success bins or success rates, and does not explain the algorithmic cause of any image difference. These are observations of the **saved primary estimates**, not raw native updates or a fresh reconstruction experiment.

The existing comparison authorization and three original/extension seals were validated before the allowed generator metadata was read. No GT NPZ or candidate NPZ was opened for this review. Family/variant/acquisition/grid/noise labels come from the authorized original `sealed_truth/fresh48/GENERATOR_MANIFEST.json`; geometry descriptions come from the actual galleries. The already-running physical-prefix analysis was not repeated.

## Reading the evidence

Each linked page contains four rows, top to bottom. Columns are **GT, CIE, Dual-CIE, G19, DR, A+**. Every row shares its GT-derived physical-χ color range, so within-row comparisons are meaningful; colors across different rows do not share an absolute range. The original displayed above-row-max counts indicate clipping, not a newly chosen threshold. “Bright,” “dim,” “hole,” “bridge” and “retains” below are qualitative display observations, not a segmentation rule or a new success classification. Weak halos can make connectivity ambiguous, and saturated yellow conceals peak amplitudes.

All 12 PNGs were directly inspected at the tool's standard resized display (approximately 1886–1888 by 1343–1344 pixels from approximately 2890–2895 by 2059 source pixels). No replacement figures or edited images were made. Very small features require the stated uncertainty; the original linked PNG remains the evidence. An earlier two-page full-resolution display attempt was truncated; the subsequent one-page-at-a-time review covered all pages. Its cost and human/model review CPU/wall were not measured and remain unknown. No new forward, adjoint, LU, RHS, scoring or inversion work was done.

## Descriptive family notes

Family labels are metadata, not inferred from the reconstructions. Each row below includes every one of that family's six cases; it is not a family success ranking.

| Original family | All case suffixes | Description and exceptions |
| --- | --- | --- |
| `mixed_scale_disconnected` | 0000, 0004, 0011, 0019, 0036, 0043 | A+ and DR show very similar four-island layouts in 0000, 0004, 0011, 0019 and 0043, with diffuse weak right-side components. Case 0036 differs from those contrast-separated targets: both methods make the small upper-left/lower-right bright GT islands dim and diffuse. Recognizable island positions therefore do not imply recovered amplitudes. |
| `multi_island_ellipses` | 0001, 0013, 0025, 0027, 0034, 0047 | A+ and DR have closely similar broad lobe locations in 0001, 0013, 0025, 0027, 0034 and 0047. Small/weak lobes and narrow gaps blur, often leaving diffuse bands between brighter regions. The images do not firmly establish component separation where GT lobes already approach or overlap. |
| `open_notched_porous_body` | 0002, 0006, 0007, 0014, 0022, 0040 | The family contains different visible outcomes. In 0002 DR retains the open arms and left hole while A+ forms displaced loops; 0007 is distorted in both. In 0006, 0022 and 0040 both retain the main notch/arms with smoothing. In 0014 both retain the notch but fill the small upper hole into a depression. A single family-wide topology statement is unsupported. |
| `continuous_random_field` | 0003, 0010, 0012, 0023, 0031, 0032 | In 0003 A+ retains a smoother broad field than the pocketed DR image; 0012 and 0023 are broadly similar across methods. Case 0010 is heavily clipped on its low row scale. In 0031 and 0032 A+ retains the full continuous body while DR becomes sparse points or a false frame. These images concern graded fields, not binary pore topology. |
| `perforated_slab_closed_holes` | 0005, 0008, 0016, 0024, 0028, 0033 | A+ retains the slab/hole organization in 0005 and 0024 while DR forms displaced fragments. In 0016 neither retains the slab. Both broadly retain the other slabs (0008, 0028, 0033), but small holes become shallow/enlarged depressions, and 0033's upper hole nearly opens toward the blurred top boundary. Exact hole counts after thresholding were not computed. |
| `eccentric_multilayers` | 0009, 0015, 0018, 0030, 0035, 0042 | In 0009, 0015, 0018, 0030 and 0035 A+ and DR broadly retain outer layers while tiny central openings/depressions fill or broaden. DR has more visible mottling in 0030. In 0042 DR retains a layered oval whereas A+, CIE, Dual-CIE and G19 show shifted arcs and false blocks. Outer-ring recognition does not establish recovery of the inner small feature. |
| `random_level_patchwork` | 0017, 0029, 0037, 0039, 0041, 0044 | A+ and DR retain broad sector patterns in 0017 and 0039, with A+ smoother and DR more mottled; 0041 is strongly saturated. Cases 0029, 0037 and 0044 have false rings/fragments or near-empty interiors in both A+ and DR. In 0029 A+ fills more area than sparse DR but still places bright material and holes incorrectly; greater occupancy is not treated as a recovery claim. |
| `thin_bridge_network` | 0020, 0021, 0026, 0038, 0045, 0046 | A+ and DR broadly retain nodes and the stronger right bridge in 0020, 0021, 0026, 0045 and 0046, but the weak left bridge remains diffuse or apparently interrupted. In 0038 DR/G19 retain the main network while A+ has displaced blocks. Several CIE/Dual images are fragmented even where A+/DR/G19 retain main geometry. No thresholded bridge-continuity score was introduced. |

## Per-case observations

### Page 01

[Original diagnostic gallery](../comparison48/FIGURES/diagnostic_row_scale_page01.png) — rows 0000–0003, top to bottom.

**r13_fresh_0000 — page 01, row 1.** `mixed_scale_disconnected`, variant 4; acquisition `half12`; relative noise RMS 0.01; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page01.png)

- **Displayed GT:** Four separated islands: two bright left islands, larger dim upper-right island and small dim lower-right island.
- **A+:** Retains all four islands; left bright components and right dim components remain separated, with softened edges.
- **DR:** Very similar four-island layout; small lower-right component is also diffuse.
- **Amplitude:** Both retain the strong left/right contrast difference; both have clipped bright pixels on the row scale.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 give nearly the same broad island layout.
- **Uncertainty:** Fine boundary and peak-amplitude differences are small at this display scale.

**r13_fresh_0001 — page 01, row 2.** `multi_island_ellipses`, variant 1; acquisition `full12`; relative noise RMS 0.02; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page01.png)

- **Displayed GT:** Three unequal islands: bright vertical left ellipse, bright lower-right island and faint upper-right island.
- **A+:** Recovers the three locations, including the faint upper-right island; diffuse low-level structure appears between the bright islands.
- **DR:** Same three-part arrangement, with a slightly broader bright-left footprint and diffuse central/background structure.
- **Amplitude:** Strong and faint island hierarchy retained by both; saturated peaks prevent precise visual amplitude comparison.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 are similarly close in overall geometry.
- **Uncertainty:** Faint bridges may be low-amplitude blur rather than actual recovered connectivity.

**r13_fresh_0002 — page 01, row 3.** `open_notched_porous_body`, variant 4; acquisition `half12`; relative noise RMS 0.01; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page01.png)

- **Displayed GT:** Bent open body with two long arms, a left-side hole and a broad lower dim extension.
- **A+:** Replaces the body with displaced high-valued loops and scattered blocks; the original open-arm geometry is poorly retained.
- **DR:** Retains the open-arm outline and left hole, with rounded edges and a broad dim lower extension.
- **Amplitude:** A+ contains extensive above-row-max saturation; DR also clips some bright body pixels.
- **CIE / Dual-CIE / G19 context:** CIE and Dual-CIE also form spurious loops; G19 is sparse and fragmented.
- **Uncertainty:** Saturation conceals the magnitude of A+ overshoots, but the geometry difference is clear.

**r13_fresh_0003 — page 01, row 4.** `continuous_random_field`, variant 3; acquisition `full6`; relative noise RMS 0.02; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page01.png)

- **Displayed GT:** Rounded continuous body with a smooth dim upper region and brighter lower/right interior.
- **A+:** Retains the rounded support and broad upper-to-lower intensity trend, with some shallow interior mottling.
- **DR:** Retains the support but introduces stronger interior dark pockets and a bright lower rim.
- **Amplitude:** A+ has a smoother amplitude pattern closer to the displayed broad trend; DR has a few clipped pixels.
- **CIE / Dual-CIE / G19 context:** CIE and Dual-CIE are smooth with mild interior pockets; G19 is much more mottled.
- **Uncertainty:** No binary topology is asserted for a continuous field; fine texture is display-dependent.

### Page 02

[Original diagnostic gallery](../comparison48/FIGURES/diagnostic_row_scale_page02.png) — rows 0004–0007, top to bottom.

**r13_fresh_0004 — page 02, row 1.** `mixed_scale_disconnected`, variant 5; acquisition `half12`; relative noise RMS 0.02; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page02.png)

- **Displayed GT:** Four separated islands with two bright left components and two faint right components.
- **A+:** Preserves the four locations; lower-left island has a soft halo, and the faint right islands are broadened.
- **DR:** Nearly the same four-island geometry; faint upper-right island is slightly more irregular.
- **Amplitude:** Large left/right amplitude difference survives in both; bright peaks clip.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 produce similar footprints.
- **Uncertainty:** Differences are subtle relative to pixelation and low-amplitude background.

**r13_fresh_0005 — page 02, row 2.** `perforated_slab_closed_holes`, variant 5; acquisition `half12`; relative noise RMS 0.02; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page02.png)

- **Displayed GT:** Nearly square slab with three closed holes: one upper-right and two lower holes.
- **A+:** Retains slab and three holes, with rounded/distorted boundaries and an extra shallow upper-left depression.
- **DR:** Produces displaced saturated bars and branching fragments; the square slab and three-hole arrangement are lost.
- **Amplitude:** A+ retains a coherent bright slab with some clipped pixels; DR has extensive clipped blocks.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 retain the three-hole slab much like A+.
- **Uncertainty:** The small upper-left depression is not a GT hole; row-scale clipping prevents exact peak comparison.

**r13_fresh_0006 — page 02, row 3.** `open_notched_porous_body`, variant 3; acquisition `full6`; relative noise RMS 0.02; forward/inverse grid 36/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page02.png)

- **Displayed GT:** Rotated open U-like body with a small closed hole near its lower vertex and a dim band beside the upper arm.
- **A+:** Retains the open notch, lower hole and two main arms; boundaries are softened.
- **DR:** Retains the same topology with somewhat more mottled bright arms and thicker lower-hole rim.
- **Amplitude:** Both preserve bright arms against the dim side band, with saturated bright pixels.
- **CIE / Dual-CIE / G19 context:** CIE and Dual-CIE are close to A+; G19 shows somewhat stronger interior mottling.
- **Uncertainty:** Hole size and arm thickness cannot be judged precisely at the coarse pixel scale.

**r13_fresh_0007 — page 02, row 4.** `open_notched_porous_body`, variant 1; acquisition `full12`; relative noise RMS 0.02; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page02.png)

- **Displayed GT:** Open rectangular C-like body with a small upper-left hole and a dim lower band.
- **A+:** Shows a partial left vertical/rim structure but displaced top/right blocks and broken lower arm; original C geometry is not recovered.
- **DR:** Mostly sparse bright fragments, with little continuous left or lower arm.
- **Amplitude:** A+ and DR both have clipped fragments; A+ has more extended high-valued structure.
- **CIE / Dual-CIE / G19 context:** CIE and Dual-CIE form large spurious loops; G19 is sparse like DR.
- **Uncertainty:** A+ retains more left-side outline than DR, but neither provides reliable hole/notch geometry.

### Page 03

[Original diagnostic gallery](../comparison48/FIGURES/diagnostic_row_scale_page03.png) — rows 0008–0011, top to bottom.

**r13_fresh_0008 — page 03, row 1.** `perforated_slab_closed_holes`, variant 1; acquisition `full12`; relative noise RMS 0.02; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page03.png)

- **Displayed GT:** Diamond-shaped slab with three closed holes of unequal size.
- **A+:** Retains the diamond and all three hole locations, with enlarged/rounded holes and softened perimeter.
- **DR:** Retains the three holes and slab, but the upper/right material is more mottled and the left hole less sharply dark.
- **Amplitude:** Both have clipped bright slab pixels; A+ looks more spatially uniform within the material.
- **CIE / Dual-CIE / G19 context:** CIE and Dual-CIE are close to A+; G19 has a somewhat irregular bright lower region.
- **Uncertainty:** Small differences in hole area are limited by coarse pixels and blur.

**r13_fresh_0009 — page 03, row 2.** `eccentric_multilayers`, variant 3; acquisition `full6`; relative noise RMS 0.02; forward/inverse grid 36/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page03.png)

- **Displayed GT:** Eccentric layered circular body with a bright outer ring, dim inner region and a small central dark opening.
- **A+:** Retains the outer annulus and dim interior, but fills the small central opening with a brighter spot.
- **DR:** Very similar layered outline and filled central opening.
- **Amplitude:** Both retain a bright lower/left rim and weaker upper-right rim; peaks clip.
- **CIE / Dual-CIE / G19 context:** CIE and Dual-CIE are similar; G19 adds stronger rim/interior fluctuations.
- **Uncertainty:** The displayed center is unresolved as a hole by either A+ or DR.

**r13_fresh_0010 — page 03, row 3.** `continuous_random_field`, variant 0; acquisition `full12`; relative noise RMS 0.01; forward/inverse grid 36/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page03.png)

- **Displayed GT:** Rounded filled body with very weak visible internal variation on this row scale.
- **A+:** Retains the circular footprint, with shallow interior depressions and a diffuse halo.
- **DR:** Very similar footprint and depressions.
- **Amplitude:** Many pixels exceed the low GT-derived maximum in all methods, masking detailed amplitude differences.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 look very similar in support and broad amplitude.
- **Uncertainty:** Saturation makes this page unsuitable for judging subtle internal amplitude recovery.

**r13_fresh_0011 — page 03, row 4.** `mixed_scale_disconnected`, variant 3; acquisition `full6`; relative noise RMS 0.02; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page03.png)

- **Displayed GT:** Four separated islands: two strong left components, larger dim upper-right and small dim lower-right.
- **A+:** Retains all four islands, with mild smearing around the lower-left and upper-right components.
- **DR:** Same four components, with a slightly longer rightward tail on the lower-left island.
- **Amplitude:** Both preserve the bright/dim hierarchy; a few bright pixels clip.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 provide nearly the same geometry.
- **Uncertainty:** Small shape differences are within a few coarse pixels; no clear topological difference is visible.

### Page 04

[Original diagnostic gallery](../comparison48/FIGURES/diagnostic_row_scale_page04.png) — rows 0012–0015, top to bottom.

**r13_fresh_0012 — page 04, row 1.** `continuous_random_field`, variant 1; acquisition `full12`; relative noise RMS 0.02; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page04.png)

- **Displayed GT:** Rounded continuous body with a bright upper-right region, a smaller bright left region and a dim central vertical swath.
- **A+:** Retains both bright regions and dim central swath, with smooth support and moderate boundary blur.
- **DR:** Nearly the same broad field, with slightly more visible small interior fluctuations.
- **Amplitude:** Both preserve the nonuniform amplitude pattern without row-scale clipping.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 are close in broad structure.
- **Uncertainty:** Small numerical differences are not reliably separable from this display; no binary topology claim is appropriate.

**r13_fresh_0013 — page 04, row 2.** `multi_island_ellipses`, variant 5; acquisition `half12`; relative noise RMS 0.02; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page04.png)

- **Displayed GT:** Isolated bright lower-left island and a clustered upper/right group containing bright and dim overlapping lobes.
- **A+:** Preserves the isolated island and the clustered lobe locations; weak lobes blur together.
- **DR:** Very similar arrangement and blurred upper cluster.
- **Amplitude:** Both retain bright small cores and dim broad lobes, with clipped core pixels.
- **CIE / Dual-CIE / G19 context:** All three auxiliary comparators show nearly the same clustered structure.
- **Uncertainty:** Connectivity of the dim lobes is visually ambiguous because their GT footprints already approach or overlap.

**r13_fresh_0014 — page 04, row 3.** `open_notched_porous_body`, variant 0; acquisition `full12`; relative noise RMS 0.01; forward/inverse grid 36/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page04.png)

- **Displayed GT:** Inverted U-like body with a narrow central open notch and a small closed upper-right hole.
- **A+:** Retains the two legs and open notch but broadens the top bridge; the small closed hole becomes a shallow depression.
- **DR:** Same broad inverted-U shape, with the small hole also filled and the notch blurred.
- **Amplitude:** Both have extensive clipping relative to this low-amplitude GT row, obscuring top-bridge variation.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 exhibit the same missing small-hole feature.
- **Uncertainty:** The notch is visible, but its width and the degree of hole filling cannot be quantified from the gallery.

**r13_fresh_0015 — page 04, row 4.** `eccentric_multilayers`, variant 4; acquisition `half12`; relative noise RMS 0.01; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page04.png)

- **Displayed GT:** Vertically elongated eccentric layers with bright outer ring, dim inner layer and a tiny central opening.
- **A+:** Retains outer ring and dim interior but fills the tiny central opening with a compact brighter region.
- **DR:** Similar layers; the left/lower rim is somewhat more irregular.
- **Amplitude:** Both preserve the bright outer/lower ring and weaker upper-left ring, with clipped peaks.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 have almost identical broad layered geometry.
- **Uncertainty:** The central opening is not resolved by any method at this scale.

### Page 05

[Original diagnostic gallery](../comparison48/FIGURES/diagnostic_row_scale_page05.png) — rows 0016–0019, top to bottom.

**r13_fresh_0016 — page 05, row 1.** `perforated_slab_closed_holes`, variant 4; acquisition `half12`; relative noise RMS 0.01; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page05.png)

- **Displayed GT:** Rectangular slab with three closed holes: two upper holes and a larger lower-left hole.
- **A+:** Shows sparse arcs and bars rather than a filled slab; the three-hole arrangement is absent.
- **DR:** Also loses the slab, producing a different set of bright diagonal and vertical fragments.
- **Amplitude:** Both contain clipped bright fragments and large near-empty regions where GT material is present.
- **CIE / Dual-CIE / G19 context:** CIE and Dual-CIE are almost identical to A+; G19 is also sparse and incomplete.
- **Uncertainty:** Fragment-level differences do not establish reliable recovered topology for any method.

**r13_fresh_0017 — page 05, row 2.** `random_level_patchwork`, variant 3; acquisition `full6`; relative noise RMS 0.02; forward/inverse grid 36/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page05.png)

- **Displayed GT:** Rounded patchwork body with distinct bright upper-left and lower-right sectors, a dim upper-right sector and an intermediate lower-left sector.
- **A+:** Retains the four broad sector arrangement and outer footprint, with smoothed boundaries and mild interior variation.
- **DR:** Retains the sectors but introduces stronger mottling and a bright lower rim.
- **Amplitude:** A+ is smoother within the sectors; both have some bright pixels above the row maximum.
- **CIE / Dual-CIE / G19 context:** CIE and Dual-CIE contain spurious rings and fragmented high values; G19 retains only disconnected arcs.
- **Uncertainty:** Internal sector boundaries are blurred, so exact regional amplitudes are not inferred.

**r13_fresh_0018 — page 05, row 3.** `eccentric_multilayers`, variant 1; acquisition `full12`; relative noise RMS 0.02; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page05.png)

- **Displayed GT:** Tilted eccentric layered body with bright outer ring and a tiny dark center inside an inner layer.
- **A+:** Retains the outer ring and broad inner layers; the tiny central opening is filled, with a separate dim patch to its right.
- **DR:** Very similar ring and filled center, with slightly more interior variation.
- **Amplitude:** Both reproduce the bright lower rim and dim upper rim with clipped peaks.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 are close in support and layered pattern.
- **Uncertainty:** The tiny central hole is below the visible recovery resolution here.

**r13_fresh_0019 — page 05, row 4.** `mixed_scale_disconnected`, variant 1; acquisition `full12`; relative noise RMS 0.02; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page05.png)

- **Displayed GT:** Four separated islands, with small bright upper-left and elongated bright lower-left islands plus two faint right islands.
- **A+:** Retains both bright islands and the upper-right weak component; the lower-right component is very faint and diffuse.
- **DR:** Similar bright components; the two weak right islands are also low contrast and diffuse.
- **Amplitude:** Both preserve large amplitude separation; bright peaks clip while the weakest component is barely visible.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 show nearly the same pattern.
- **Uncertainty:** Presence and shape of the weakest lower-right component cannot be judged firmly from this scale.

### Page 06

[Original diagnostic gallery](../comparison48/FIGURES/diagnostic_row_scale_page06.png) — rows 0020–0023, top to bottom.

**r13_fresh_0020 — page 06, row 1.** `thin_bridge_network`, variant 0; acquisition `full12`; relative noise RMS 0.01; forward/inverse grid 36/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page06.png)

- **Displayed GT:** Three bright node regions joined by two thin sloping bridges in an open triangular arrangement.
- **A+:** Retains the node arrangement and right sloping bridge, but broadens/weakens the left bridge with an apparent gap near its upper end.
- **DR:** Very similar nodes and bridges, also with a weak/interrupted upper-left connection.
- **Amplitude:** Bright nodes clip; thin GT bridges become broader lower-amplitude bands in both.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 show the same broadening and weak left connection.
- **Uncertainty:** Visual continuity depends on low-intensity pixels; no thresholded connectivity claim is made.

**r13_fresh_0021 — page 06, row 2.** `thin_bridge_network`, variant 2; acquisition `full6`; relative noise RMS 0.01; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page06.png)

- **Displayed GT:** Two bright nodes, a dim lower-right node, bright upper-to-right thin bridge and a much weaker left bridge.
- **A+:** Retains both bright nodes and the right sloping bridge; the weak left bridge and dim lower-right node are diffuse.
- **DR:** Almost identical layout and weak-bridge visibility.
- **Amplitude:** Both retain bright-node versus dim-node contrast, with some clipping at the bright nodes.
- **CIE / Dual-CIE / G19 context:** Dual-CIE and G19 are similar to A+/DR; CIE produces displaced saturated fragments.
- **Uncertainty:** The low-contrast left bridge is not reliably resolved as continuous by either A+ or DR.

**r13_fresh_0022 — page 06, row 3.** `open_notched_porous_body`, variant 2; acquisition `full6`; relative noise RMS 0.01; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page06.png)

- **Displayed GT:** Slanted open C-like body with a small upper-left closed hole and a dim lower-side band.
- **A+:** Retains the open notch, closed hole and sloping arms, with rounded edges.
- **DR:** Retains the same geometry but shows more mottled arm intensities and a broader lower-side halo.
- **Amplitude:** Both have clipped bright body pixels; A+ is smoother across the arms.
- **CIE / Dual-CIE / G19 context:** CIE and Dual-CIE are close to A+; G19 has much stronger interior fluctuations.
- **Uncertainty:** Small changes in hole size and arm thickness remain display-limited.

**r13_fresh_0023 — page 06, row 4.** `continuous_random_field`, variant 2; acquisition `full6`; relative noise RMS 0.01; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page06.png)

- **Displayed GT:** Rounded continuous field with bright side lobes, a dim central vertical region and weaker upper/lower regions.
- **A+:** Retains the broad side lobes and dim center with a smooth footprint.
- **DR:** Retains the same field but with slightly stronger local fluctuations and a few clipped side pixels.
- **Amplitude:** A+ preserves the broad amplitude ordering without clipping; DR clips a few high pixels.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 are also close in broad field structure.
- **Uncertainty:** This continuous field has no discrete pore/bridge topology to infer from the image.

### Page 07

[Original diagnostic gallery](../comparison48/FIGURES/diagnostic_row_scale_page07.png) — rows 0024–0027, top to bottom.

**r13_fresh_0024 — page 07, row 1.** `perforated_slab_closed_holes`, variant 3; acquisition `full6`; relative noise RMS 0.02; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page07.png)

- **Displayed GT:** Slanted rectangular slab with three closed holes: a large upper-right hole and two lower holes.
- **A+:** Retains the slab and three depressions/holes, with a diffuse left boundary and distorted hole shapes.
- **DR:** Produces disconnected arcs and high-valued blocks; the filled slab and three-hole layout are absent.
- **Amplitude:** A+ remains spatially filled but nonuniform; DR has large clipped fragments and empty interior regions.
- **CIE / Dual-CIE / G19 context:** CIE and Dual-CIE retain the slab with more pronounced internal/boundary distortion; G19 is sparse.
- **Uncertainty:** The two smaller holes are blurred; their exact sizes cannot be inferred visually.

**r13_fresh_0025 — page 07, row 2.** `multi_island_ellipses`, variant 2; acquisition `full6`; relative noise RMS 0.01; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page07.png)

- **Displayed GT:** Compact cluster containing two bright lobes along a rising diagonal with dim material between/above them.
- **A+:** Retains the two bright lobes and connecting dim region, with softened edges.
- **DR:** Nearly identical cluster and lobe spacing.
- **Amplitude:** Both preserve relative brightness and have clipped core pixels.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 look similar in cluster geometry.
- **Uncertainty:** Overlap in the GT cluster makes component-count or connectivity claims inappropriate.

**r13_fresh_0026 — page 07, row 3.** `thin_bridge_network`, variant 5; acquisition `half12`; relative noise RMS 0.02; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page07.png)

- **Displayed GT:** Open triangular network with two bright nodes, a dim lower-right node, strong right bridge and weaker left bridge.
- **A+:** Retains node positions and strong right bridge; left bridge is weak/diffuse and lower-right node broad.
- **DR:** Very similar structure and weak left-bridge visibility.
- **Amplitude:** Both preserve bright/dim node contrast; a few bright pixels clip.
- **CIE / Dual-CIE / G19 context:** G19 is close to A+/DR; CIE and Dual-CIE show displaced saturated blocks and broken bridges.
- **Uncertainty:** Continuity of the weak left bridge remains uncertain at displayed resolution.

**r13_fresh_0027 — page 07, row 4.** `multi_island_ellipses`, variant 4; acquisition `half12`; relative noise RMS 0.01; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page07.png)

- **Displayed GT:** Cluster with a dim upper-left lobe, bright lower-left lobe and a larger irregular bright right group.
- **A+:** Retains the principal lobes but introduces a diffuse band between lower-left and right regions.
- **DR:** Similar principal lobes and diffuse connecting band.
- **Amplitude:** Both retain the dim-versus-bright hierarchy and clipped bright regions.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 yield almost the same smoothed cluster.
- **Uncertainty:** Apparent joining across small gaps is blur-dependent; individual overlapping right lobes are not separately resolved.

### Page 08

[Original diagnostic gallery](../comparison48/FIGURES/diagnostic_row_scale_page08.png) — rows 0028–0031, top to bottom.

**r13_fresh_0028 — page 08, row 1.** `perforated_slab_closed_holes`, variant 0; acquisition `full12`; relative noise RMS 0.01; forward/inverse grid 36/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page08.png)

- **Displayed GT:** Rotated slab with three closed holes, including a small left hole and larger upper/lower holes.
- **A+:** Retains slab and larger hole locations; the small left hole becomes a broad shallow depression.
- **DR:** Almost identical blurred three-hole pattern.
- **Amplitude:** Both clip many bright slab pixels relative to the low row maximum.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 have similar geometry and small-hole smoothing.
- **Uncertainty:** The small GT hole is not sharply resolved; saturation limits interior-amplitude interpretation.

**r13_fresh_0029 — page 08, row 2.** `random_level_patchwork`, variant 2; acquisition `full6`; relative noise RMS 0.01; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page08.png)

- **Displayed GT:** Rounded patchwork body with dim left sector, intermediate upper-central sector and bright right/lower-right sectors.
- **A+:** Forms a roughly body-sized but misplaced bright left/lower shell with spurious internal holes and bright grid-like fragments.
- **DR:** Leaves only a few isolated points/arcs; most of the body is absent.
- **Amplitude:** A+ has extensive clipped bright structure with incorrect amplitude ordering; DR is mostly near background.
- **CIE / Dual-CIE / G19 context:** CIE and Dual-CIE contain displaced bars and loops; G19 is nearly as sparse as DR.
- **Uncertainty:** A+ occupies more of the image but does not recover the displayed sector organization; occupancy alone is not treated as recovery.

**r13_fresh_0030 — page 08, row 3.** `eccentric_multilayers`, variant 2; acquisition `full6`; relative noise RMS 0.01; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page08.png)

- **Displayed GT:** Eccentric layered circular body with bright upper/right/lower rim, weak left rim and a small dark inner opening.
- **A+:** Retains broad layers and asymmetric rim, but blurs the weak left side and fills the tiny inner opening.
- **DR:** Retains outer footprint but has stronger dark pockets and interior mottling.
- **Amplitude:** A+ shows a smoother bright rim; both clip some peaks.
- **CIE / Dual-CIE / G19 context:** CIE and Dual-CIE look similar to A+; G19 has strong internal bright/dark fluctuations.
- **Uncertainty:** The weak left rim and central hole are not sharply resolved; apparent gaps should not be read as exact topology.

**r13_fresh_0031 — page 08, row 4.** `continuous_random_field`, variant 4; acquisition `half12`; relative noise RMS 0.01; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page08.png)

- **Displayed GT:** Rounded continuous field with bright upper/side regions and a long dim central/lower region.
- **A+:** Retains full rounded support and broad amplitude pattern.
- **DR:** Collapses to a few isolated bright points; broad support and internal pattern are absent.
- **Amplitude:** A+ has no row-scale clipping and retains the smooth intensity range; DR has a few saturated points.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 retain the broad field much like A+.
- **Uncertainty:** The visible contrast is clear, but the image alone cannot explain the differing trajectories.

### Page 09

[Original diagnostic gallery](../comparison48/FIGURES/diagnostic_row_scale_page09.png) — rows 0032–0035, top to bottom.

**r13_fresh_0032 — page 09, row 1.** `continuous_random_field`, variant 5; acquisition `half12`; relative noise RMS 0.02; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page09.png)

- **Displayed GT:** Rounded continuous body with bright upper/left regions and a dim central/lower region.
- **A+:** Retains the broad support and amplitude pattern, with mild boundary blur.
- **DR:** Forms a saturated polygonal frame and large lower block rather than the smooth filled field.
- **Amplitude:** A+ has only a few clipped pixels; DR has extensive saturation and empty regions.
- **CIE / Dual-CIE / G19 context:** CIE and Dual-CIE preserve the smooth field; G19 retains support but is more mottled.
- **Uncertainty:** The gallery shows a clear final-state difference, without identifying its cause.

**r13_fresh_0033 — page 09, row 2.** `perforated_slab_closed_holes`, variant 2; acquisition `full6`; relative noise RMS 0.01; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page09.png)

- **Displayed GT:** Upright rectangular slab with three holes: small upper-left, large right-middle and small lower-left.
- **A+:** Retains the slab and right hole; upper hole enlarges toward the top edge, and lower-left hole becomes a shallow depression.
- **DR:** Same broad geometry with more interior mottling and an indistinct lower-left hole.
- **Amplitude:** Both have extensive clipped bright material; A+ is smoother in the slab interior.
- **CIE / Dual-CIE / G19 context:** CIE and Dual-CIE look similar to A+; G19 is more mottled.
- **Uncertainty:** Whether the blurred upper hole remains closed cannot be established without choosing a new threshold.

**r13_fresh_0034 — page 09, row 3.** `multi_island_ellipses`, variant 3; acquisition `full6`; relative noise RMS 0.02; forward/inverse grid 36/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page09.png)

- **Displayed GT:** Two bright separated upper/right islands with a dim lower-left island near the upper-left bright one.
- **A+:** Retains the three principal locations but makes the dim island more compact and adds a diffuse horizontal band between bright islands.
- **DR:** Very similar three-region pattern and diffuse band.
- **Amplitude:** Both preserve bright cores but shift the weak-island pattern; bright peaks clip.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 have almost the same geometry.
- **Uncertainty:** The diffuse band does not establish a true recovered connection between GT islands.

**r13_fresh_0035 — page 09, row 4.** `eccentric_multilayers`, variant 0; acquisition `full12`; relative noise RMS 0.01; forward/inverse grid 36/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page09.png)

- **Displayed GT:** Elongated layered oval with a bright broad outer body and a small dark eccentric inner depression/opening.
- **A+:** Retains the oval but broadens and fills the inner dark feature into a larger green region.
- **DR:** Almost identical oval and softened interior depression.
- **Amplitude:** Both clip many outer-body pixels on the low-amplitude row scale and make the center visibly brighter than GT.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 show the same general interior smoothing.
- **Uncertainty:** Exact layer boundaries and peak amplitudes are obscured by blur and clipping.

### Page 10

[Original diagnostic gallery](../comparison48/FIGURES/diagnostic_row_scale_page10.png) — rows 0036–0039, top to bottom.

**r13_fresh_0036 — page 10, row 1.** `mixed_scale_disconnected`, variant 0; acquisition `full12`; relative noise RMS 0.01; forward/inverse grid 36/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page10.png)

- **Displayed GT:** Four separated bright islands of unequal size, including small upper-left and lower-right components.
- **A+:** Retains four broad locations, but suppresses the small upper-left and lower-right islands into diffuse dim patches; larger islands are reshaped.
- **DR:** Nearly the same pattern of two bright large islands and two weaker diffuse small islands.
- **Amplitude:** Both misrepresent the similar bright GT amplitudes by creating a strong large/small amplitude split; large-island peaks clip.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 show the same amplitude/size distortion.
- **Uncertainty:** Positions remain recognizable, but this should not be described as preserving the GT amplitude hierarchy.

**r13_fresh_0037 — page 10, row 2.** `random_level_patchwork`, variant 4; acquisition `half12`; relative noise RMS 0.01; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page10.png)

- **Displayed GT:** Rounded asymmetric patchwork body, mostly bright, with dim right and lower sectors.
- **A+:** Produces displaced saturated fragments, loops and a lower arc; the filled body and sector arrangement are lost.
- **DR:** Produces an oversized lower/right arc and isolated upper-left blocks, also losing the body.
- **Amplitude:** Both contain extensive clipping and large false empty regions.
- **CIE / Dual-CIE / G19 context:** CIE and Dual-CIE form large spurious loops; G19 is sparse and fragmentary.
- **Uncertainty:** No reliable regional topology or amplitude recovery can be inferred for any method from these final images.

**r13_fresh_0038 — page 10, row 3.** `thin_bridge_network`, variant 1; acquisition `full12`; relative noise RMS 0.02; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page10.png)

- **Displayed GT:** Open triangular node-and-bridge network with bright top/lower-left nodes, strong right bridge and weak left bridge.
- **A+:** Produces displaced bright loops/blocks and loses the node-and-bridge layout.
- **DR:** Retains bright node positions and the strong right bridge; weak left bridge remains diffuse.
- **Amplitude:** A+ has extensive saturated false structure; DR retains the bright-node versus dim-background pattern with few clipped pixels.
- **CIE / Dual-CIE / G19 context:** G19 is similar to DR; CIE and Dual-CIE form spurious large loops/blocks.
- **Uncertainty:** DR's weak left connection is not sharply resolved, despite recognizable main geometry.

**r13_fresh_0039 — page 10, row 4.** `random_level_patchwork`, variant 1; acquisition `full12`; relative noise RMS 0.02; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page10.png)

- **Displayed GT:** Rounded patchwork body with bright central/left material, dim right sector and an irregular dim inner/lower-left region.
- **A+:** Retains outer body and regional amplitude ordering, with smoothed internal boundaries.
- **DR:** Retains broad regional structure but adds stronger small-scale mottling and an irregular bright lower rim.
- **Amplitude:** Both clip bright pixels; A+ is smoother within broad regions.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 look close to A+ in broad geometry.
- **Uncertainty:** Internal patch borders are diffuse; no new segmentation is applied.

### Page 11

[Original diagnostic gallery](../comparison48/FIGURES/diagnostic_row_scale_page11.png) — rows 0040–0043, top to bottom.

**r13_fresh_0040 — page 11, row 1.** `open_notched_porous_body`, variant 5; acquisition `half12`; relative noise RMS 0.02; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page11.png)

- **Displayed GT:** Slanted open body with two long arms, a broad central notch and a small upper-right closed hole.
- **A+:** Retains notch, arms and upper-right hole, with softened edges and a shallow lower-arm depression.
- **DR:** Retains the same geometry but shows more mottled arms and a more irregular notch boundary.
- **Amplitude:** Both clip bright material; A+ is somewhat smoother inside the arms.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 preserve the same main features.
- **Uncertainty:** Small hole area and notch width cannot be assessed precisely from the gallery.

**r13_fresh_0041 — page 11, row 2.** `random_level_patchwork`, variant 0; acquisition `full12`; relative noise RMS 0.01; forward/inverse grid 36/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page11.png)

- **Displayed GT:** Rounded tilted patchwork body that appears nearly uniform and bright at the displayed row scale.
- **A+:** Retains the footprint but introduces shallow upper and lower interior depressions and a blurred halo.
- **DR:** Nearly the same footprint and depressions.
- **Amplitude:** All methods have many above-row-max pixels, obscuring subtle GT/internal amplitude differences.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 are visually close to A+/DR.
- **Uncertainty:** The row scale is strongly saturated; apparent uniformity is a display observation, not a claim about the underlying array.

**r13_fresh_0042 — page 11, row 3.** `eccentric_multilayers`, variant 5; acquisition `half12`; relative noise RMS 0.02; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page11.png)

- **Displayed GT:** Tilted layered oval with bright asymmetric outer rim, dim interior and a small central opening.
- **A+:** Produces a shifted, partly broken outer arc, a misplaced bright central bar and spurious lower blocks; inner layers are lost.
- **DR:** Retains the oval outer ring and a dim inner ring/opening, with softened and somewhat uneven rim.
- **Amplitude:** A+ has many clipped false fragments; DR has fewer clipped pixels and retains graded layers.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 show similar distorted arcs/blocks to A+.
- **Uncertainty:** DR's small central structure is blurred, but the broad layered geometry difference is visible.

**r13_fresh_0043 — page 11, row 4.** `mixed_scale_disconnected`, variant 2; acquisition `full6`; relative noise RMS 0.01; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page11.png)

- **Displayed GT:** Four separated islands with two bright left components and two faint right components.
- **A+:** Retains left bright islands and a diffuse upper-right patch; lower-right patch is weak and spread.
- **DR:** Very similar component layout and weak right-side visibility.
- **Amplitude:** Both retain a strong bright-left/dim-right difference, with clipped left peaks.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 show nearly the same island pattern.
- **Uncertainty:** The faint lower-right component is difficult to distinguish from low-level background at this scale.

### Page 12

[Original diagnostic gallery](../comparison48/FIGURES/diagnostic_row_scale_page12.png) — rows 0044–0047, top to bottom.

**r13_fresh_0044 — page 12, row 1.** `random_level_patchwork`, variant 5; acquisition `half12`; relative noise RMS 0.02; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page12.png)

- **Displayed GT:** Rounded patchwork body with a broad bright right region, dim left region and intermediate upper region.
- **A+:** Forms spurious rings/holes and misplaced bright lower/left blocks; the filled sector pattern is not retained.
- **DR:** Produces a large broken outer arc with isolated bars and little filled interior.
- **Amplitude:** A+ and DR both contain clipped false structures and large near-empty areas.
- **CIE / Dual-CIE / G19 context:** CIE and Dual-CIE are very similar to A+; G19 is sparse with an incomplete ring.
- **Uncertainty:** The different false patterns do not support a reliable recovered regional arrangement.

**r13_fresh_0045 — page 12, row 2.** `thin_bridge_network`, variant 3; acquisition `full6`; relative noise RMS 0.02; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page12.png)

- **Displayed GT:** Open triangular network with two bright nodes, dim lower-right node, strong right bridge and weaker left bridge.
- **A+:** Retains node locations and right bridge; weak left bridge is blurred and dim lower-right node spreads downward.
- **DR:** Similar network, with somewhat broader halos and a more diffuse dim lower-right node.
- **Amplitude:** Both preserve bright-node contrast and contain a few clipped pixels.
- **CIE / Dual-CIE / G19 context:** G19 is similar to A+/DR; CIE and Dual-CIE produce displaced saturated fragments.
- **Uncertainty:** Left-bridge continuity remains uncertain without imposing a threshold.

**r13_fresh_0046 — page 12, row 3.** `thin_bridge_network`, variant 4; acquisition `half12`; relative noise RMS 0.01; forward/inverse grid 24/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page12.png)

- **Displayed GT:** Open triangular network with elongated bright top node, bright lower-left node, dim lower-right node and unequal-strength bridges.
- **A+:** Retains main node locations and strong right bridge; weak left bridge is broad/dim and lower-right node is diffuse.
- **DR:** Nearly identical broad network and weak-bridge visibility.
- **Amplitude:** Both retain the bright/dim node ordering, with a few clipped bright pixels.
- **CIE / Dual-CIE / G19 context:** Dual-CIE and G19 are close to A+/DR; CIE has broken displaced fragments.
- **Uncertainty:** The weak bridge and dim node are not sharply resolved as geometrically precise structures.

**r13_fresh_0047 — page 12, row 4.** `multi_island_ellipses`, variant 0; acquisition `full12`; relative noise RMS 0.01; forward/inverse grid 36/24. [Figure](../comparison48/FIGURES/diagnostic_row_scale_page12.png)

- **Displayed GT:** Irregular left/lower cluster with a narrow downward-right extension and a separate rounded upper-right island.
- **A+:** Retains the large cluster, narrow extension and upper-right island, with rounded corners, broadened extension and diffuse material in the intervening gap.
- **DR:** Very similar arrangement, with somewhat brighter/broader extension and background halo.
- **Amplitude:** Both clip many bright pixels relative to the low row maximum.
- **CIE / Dual-CIE / G19 context:** CIE, Dual-CIE and G19 give a similar smoothed arrangement.
- **Uncertainty:** Connectivity across the upper gap is ambiguous because of blur; no thresholded component count is assigned.

## Provenance and limits

The companion CSV has 48 unique rows and exact page/row links; no aggregate winning count is derived from the descriptions. Rows were reviewed in page/case order without dropping any case. This review leaves primary choices, paid prefixes, comparisons, numerical outputs, plots and seals unchanged. It supplies descriptive evidence for the parent's scientific and visual review; it does not establish a causal mechanism, independent physical validation, continuum accuracy, or equal-cost superiority. Native-state and physical-gradient evidence are separate analyses.

Measured numerical-work counts for this review are zero (no numerical physics was called). Review/inspection CPU and wall time are unknown because no complete timer was running; file-writing or tool latency is not substituted for review cost. The earlier truncated display attempt is retained as a known unmeasured attempt. Existing plot-generation/scoring costs remain in their original receipts.
