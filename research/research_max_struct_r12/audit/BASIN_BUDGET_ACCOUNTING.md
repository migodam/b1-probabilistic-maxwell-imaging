# Basin v1 timing correction, before material scoring

The first basin runner gives a nominal144 CPU-second proposal window and96 CPU-second common corrector window. It writes the complete projection archive between those phases. That serialization is charged in actual CPU, but the v1 corrector deadline was expressed as proposal CPU after serialization +96. Consequently some final endpoints exceed240 total CPU (first remote DR ring263.1s). They MUST NOT be presented as240-CPU endpoints.

Retain complete v1 runs and their actual costs. For a common-budget comparison, independently select the lowest original Phi among history candidates whose recorded whole-body CPU<=30,60,120,240. These are actual paid candidates, never interpolated materials. Candidate arrays are indexed exactly by history candidate IDs. This public objective rule is registered before any basin material scoring. Startup imports are outside these particular internal clocks; external launcher timing remains separate, and the comparison is labeled body-CPU prefix, not whole-process CPU equivalence.

The independent-confirmation wrapper uses inherited WholeChildBudget to avoid this accounting gap. No historical v1 result is removed or overwritten. All later archive cost still counts in total execution cost even when it is not used as an algorithmic improvement claim.
