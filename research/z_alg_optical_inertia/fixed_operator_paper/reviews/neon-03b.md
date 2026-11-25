**Final Verdict: Accept.**

The 23-page revision ([paper-text.txt]) and its interval receipts address all five Major Revision points with no remaining blocking issues:

- **(a) Occupied 16/32 wide-class result (Verified):** [occupied_wide_control.json] uses the interior full-wave witness `c_solid` (relative residual $\le 4.18\times 10^{-33}$); both $s=16$ and $s=32$ saturate ($\nu=k_o=13$) while $R_y$ has inertia $(25,7)$ ($h_{\text{src}}=7$, six positive dictionary modes vs $q=13$).
- **(b) Coherent vs ordinary directions at equal fluence (Verified):** [coherent_control.json] and [energy.json] compare 11 coherent eigenbeams against 11 angular plane waves at identical fluence $1.25E_{\inf}\in[0.7385349228,0.7385349229]$, witness, and noise. Both saturate ($\nu=k_o=0$), yet coherent beams yield $R_y\succ 0$ ($\lambda_{\min}>0.00266$) whereas plane waves leave $R_y$ at $(27,5)$.
- **(c) PSD and uniform-acquisition claims (Verified):** No advantage is claimed over unrestricted PSD ([positive_tail_upper_audit.json] is scoped to the 72-weight bank; [summary.json] matches 14/32 exclusions and compares runtime only) or over uniform acquisition in realized separation (the text explicitly notes uniform completion sufficed on the first beam).
- **Independent check:** From `r=0.038890131339625975` and the 11 verified roots, recomputing $r^2\sum\beta_i^{-1}$ yields $[0.5908279382413251,0.5908279382414432]$ (matching `disk_e_inf_interval`), $1.25E_{\inf}$ matches `common_target_fluence_interval`, and $r^2/4=0.0003781105789033$ matches the coherent `source_test` spectrum $[0.0003781105789025,0.0003781105789046]$.

**Remaining material issues:** None blocking. Non-blocking scope limits—single fixed scalar weights and a fluence-independent error majorant $r^2 I_s$ (conservative at $s=16$, $\|E_{16}\|_F\le 0.02475<r$)—are already disclosed.
