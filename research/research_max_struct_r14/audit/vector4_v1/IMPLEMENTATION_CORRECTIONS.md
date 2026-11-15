# Implementation corrections and audit boundaries

The preserved attempt01 source snapshot is the exact code used for all24 original trajectories. The first runner's raw-view assignment used the best proposal-phase candidate rather than min(initial,proposal). Primary/secondary/corrector starting points already used the correct global paid incumbent, so optimization and primary selection were unaffected. Before any GT access, all24 histories/raw metadata were audited: zero saved raw predictions differ from the intended paid initial+proposal minimum. No prediction arrays, raw files or original seal were rewritten; no numerical reruns occurred. `RAW_SELECTION_AUDIT.json` records every candidate/Phi identity.

After the24 children terminated, the live runner received the explicit one-line future raw-view fix (`raw=self.best`). Original attempt01 source and outputs remain preserved. The post-seal scorer independently resolves raw from paid history before calling inherited truth loaders.

Parent explicitly approved post-seal offline scoring in `PARENT_SCORING_AUTHORIZATION.json`. Scoring validated every sealed artifact's recorded bytes/mtime before reading four old truth files. The original seal remains unchanged. Separate scoring and full-state inspection receipts record their own costs and zero forward LU/RHS.
