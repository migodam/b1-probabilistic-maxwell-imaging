# Accidental private-text search incident

Recorded UTC: 2026-11-07T04:42:00.176256+00:00. Exact command UTC was not captured in the tool output; it occurred after the root scaled-map correction request and before the worker's initial disclosure in the preserved conversation. Do not substitute this record time as the exact incident time.

The intended task was a public source/environment lookup. The actual command was:

```sh
rg -n 'projgr|2834|2841|lbfgsb' research/research_max_struct_r11 -g '*.md' -g '*.json' -g '*projection*' -g '*stop*' | head -n 40
```

Working directory: `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r11_worktree`. The search was mistakenly scoped to the entire R11 tree, including private text metadata. The retained tool output confirms one printed private line: `sealed_truth/crossgrid8_paired/generator_metadata.json:54`, field `inner`, containing one floating-point geometry parameter. This audit deliberately does not repeat the number. No case ID, seed, material array or truth image was printed in that matched line. No GT NPZ array was loaded. **The internal rg scan may have read additional private text within its broad scope; its exact read extent is unknown, so this is not a claim that only one private file was accessed.** No private file has been reopened or re-searched to prepare this audit.

The inverse worker, online inputs, algorithms and selection rules did not use this metadata. The paired8 contract was already frozen and its inversions were running; no observation, geometry, seed, method, budget or rule was changed, and no reroll/redraw occurred. The scalar was not used for any calculation or choice. Root was notified promptly and directed disclosure of this process deviation at paired8 delivery. It is incorrect to describe the whole research process as having zero private-metadata contact.

Subsequent work uses only explicit public paths: existing main48 scored CSV, corrected endpoint bookkeeping JSON, and those JSON-listed public cap/saved-gradient files. Named environment receipts may be read; no full-tree search or runtime probe will be used. This audit adds no physical computation, scoring, or science change.
