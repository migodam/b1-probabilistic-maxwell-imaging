# BGGA publication audit

Status: **READ_ONLY_AUDIT_COMPLETE**. Scope: local publication scripts, ignore rules, saved metadata, and the existing public W-SLAM packaging convention. No script, experiment, physics, publication, network readback, Git mutation, or gate edit was performed. Paths below are relative to `research/research_w_bgga/` unless prefixed `../` or `../../`. The main thread is concurrently updating sources; findings describe the inspected snapshot.

## Stale preflight claims and overwrite points

| Source | Exact stale text or behavior | Required disposition before a new publication |
|---|---|---|
| `tools/build_site.py:40` | Navigation label `G0待授权方案` | Label the resource document as the historical request, or link the current scope/accounting record. Do not present an old request as the current authorization state. |
| `tools/build_site.py:41` | Title `Research W-BGGA｜实现与资格前检查`; banner `准备工作已执行；新增全波0。G0资源授权待确认，未宣布物理/概率/速度通过。` | Rebuilding the HTML currently reinstates pre-G0 status. Use current reviewed status with separate G0 physical work, zero-new-PDE G1 replay, and unestablished speed/continuum claims. |
| `tools/prepare_publication.py:27` | Research-index card `固定公共网格与区间概率代码、32项测试完成；新增全波0，G0额度待授权。` | This is hardcoded, so rerunning the helper overwrites any manually corrected research card. Replace its stage summary after main review. |
| `tools/prepare_publication.py:32` | Repository README entry `静态Green吸收实现与资格前检查` and `32项测试通过；新物理待独立资源授权` | The 32-test count belongs to the saved preflight receipt (`CODE_CHECK.json:20–21`), not a fresh G0/G1 test attestation. Keep its historical scope or cite the actual current receipt. |
| `tools/prepare_publication.py:17–24` | Refreshes three CPU/wall/RSS fields but leaves old `COST_SUMMARY.json` status and FE/RHS fields untouched; writes `accounting_gap=... no new physical work` and a permanently pending current-report-job note | `COST_SUMMARY.json:2–4` still says `PREFLIGHT_ONLY`, `new_full_fe:0`, `new_complex_rhs:0`; line 18 still labels preflight/readback jobs. A rerun alone does not reconcile physical accounting. Close out metadata against completed wrapper receipts by phase, retaining failed jobs and timing gaps. Describe zero physics as this metadata action's scope, not the cumulative campaign. |
| `START_HERE.md:1–3,7,13,15` | `零PDE preflight`, `新全波G0等待独立资源授权`, `新FE=0、新物理RHS=0`, and pending G0 request framing | These are historical preflight assertions. They conflict with the saved G0 decision and completed G1 artifacts. |
| `REPRODUCE.md:23,25` | `当前这条物理命令会因缺少BGGA独立额度而拒绝` and `后两阶段方法/速度实验未运行` | Explain the actual closed/current authorization scopes and the completed exposed G1 replay. Keep G2/new physical work separate. Do not instruct readers to reauthorize or rerun G0 merely to read saved results. |
| `REPRODUCE.md:8–10` | Default sequence includes `tools/report_preflight.py` | That tool now explicitly refuses when physical reservations exist (`tools/report_preflight.py:12–14`). The documented default is no longer a valid closeout workflow and must not be used to regenerate current gates. |
| `GATE_DECISION.json` | `overall=G0_NUMERICAL_QUALIFIED_G1_LOGIC_NEXT`; `FINITE_PRIOR_POSTERIOR` and `JOINT_DEPENDENCE` still contain `NOT_RUN` | Main-thread scientific review must reconcile these with the saved G1 outputs. This audit deliberately does not edit or replace these decisions. |

Concrete contradiction sources, without assigning new scientific PASS: `results/G0_DECISION.json` records `QUALIFIED_DISCRETE_FULLWAVE_FOR_TESTED_CONTRACT` and physical usage 482 FE / 3934 RHS / 350.36966 CPU seconds. `results/g1_replay_001/RUN_RECEIPT.json` records 168 case-method rows with actual new FE/RHS 0; its 495600 online and 70800 reference queries are cached logical accesses. Neither a publication helper nor a zero-PDE replay erases the paid G0/historical setup.

## Git versus Release

The existing W-SLAM public convention is explicit in `../research_w_slam/REPRODUCE.md:18–24,92–113`: Git contains source, documents, figures, compact JSON/CSV, original ledgers and receipts/logs; NPZ/NPY and raw evidence are separate Release assets. `git ls-files` found zero tracked W-SLAM NPZ files while its job receipts/logs are tracked. Its `RELEASE_LAYOUT.json` lists 3847 saved files, and `RELEASE_ASSETS.json` records six stage-relative ZIPs, sizes, missing/layout changes, failed-array retention, and `scientific_completeness_claim=false`. Packaging completion is separate from scientific PASS and network upload.

BGGA's `.gitignore:1–5` currently excludes only Python bytecode, `results/.mpl/`, and locks. The enclosing repository `.gitignore` has no BGGA/WSLAM or `release_assets` rule. Do not assume the exclusion text recorded by an old packager is enforced in this checkout.

Recommended stage-local patterns, **not applied by this audit**:

```gitignore
# Raw arrays remain in the BGGA evidence Release.
*.npz
*.npy
/release_assets/
/research/delegated/g1_figures/.matplotlib/

# Dense G1 transcripts go in the same evidence Release, not the compact Git view.
/results/g1_replay_001/parent_*/M[0-6].json
/results/g1_replay_001/parent_*/OFFLINE_SCORES.json
```

Use an explicit publication inclusion list; ignore rules alone are not a packaging manifest. Do not blanket-ignore `results/`, `research/delegated/`, `*.json`, or `*.log`: that would hide compact gate evidence, source audits, failures and accounting. Suggested split:

| Git (compact, directly readable) | Release (actual saved evidence, preserved in place) |
|---|---|
| Reports, `START_HERE`, reproduction instructions, requirements, contracts and reviewed gate decisions | All actual G0/G1/normalization NPZ/NPY, including failed/late evidence and correction arrays |
| `src/`, `tools/`, substantive tests, fixed vendor source and provenance | Full per-case `M0`–`M6` JSON and `OFFLINE_SCORES.json` if excluded above: they contain dense matrices/traces, not merely scalar summaries |
| Root G1 method/parent tables and factual summary; detailed scalar CSVs; per-phase decisions, plans, seals and run/input receipts | Full original checkpoint/query transcripts and original raw signed comparisons, plus normalization repair as an additional distinct group |
| G0 compact results, source/cache map, source-access attestations and ownership/authorization records | Any later offline-extra raw arrays only if actually produced and explicitly included; no fabricated members for NOT_RUN stages |
| Actual PNG/PDF figures and their compact data/manifest; ledger/failure records and reasonably sized wrapper receipt/log evidence | Large raw job streams if necessary, with their locations and recorded scope visible in Git |

The dense-data split matters: `results/g1_replay_001` occupies about 126 MiB in the inspected checkout. In one case, `M4.json` is 909 KiB, `M5.json` 925 KiB, and `OFFLINE_SCORES.json` 2.3 MiB. Excluding only `.npz` leaves substantial array-like JSON in Git. However, `tools/g1_figures.py` needs these JSON transcripts, so the Release must include them and the public instructions must restore their exact paths before table/figure regeneration.

Retain the old convention of stage-relative members, source/archive byte counts and explicit missing entries; no added checksum step. A proposed BGGA archive should unzip directly into `research/research_w_bgga/`, without another enclosing directory. Keep the existing B1 repository/site and preserve the historical W-SLAM Release. Reuse its online bank as an explicit dependency instead of silently regenerating or duplicating paid physics.

Do **not** directly run W-SLAM's packager against BGGA without adaptation: `../research_w_slam/tools/package_release.py:15` selects only NPZ/NPY/PNG, so it omits BGGA's dense JSON transcripts; its final prerequisites require W-SLAM bank/inference paths, group names are W-SLAM-specific, and its archive names/unpack root are hardcoded (`:135–147`). Reuse the convention, not its unchanged executable.

## Public reproduction and readback gaps

1. **Old bank location.** `CACHE_MAP.json:4–7,13,...` preserves absolute `/Volumes/...` provenance. Current `tools/g1_replay.py:81–100,370–384` now supports explicit `--old-bank-root`; its default still uses the original machine path. Document the flag rather than rewriting historical provenance or using automatic fallback discovery. The recorded `research-w-slam-controlled_E2E3.zip` contains `results/e2e3/finite_bank/online/{FINITE_GRID,FINITE_BANK,OBSERVATIONS}.npz` and sibling `offline/FINITE_PRIOR_TRUTHS_AND_NOISE.npz`. Extract it into `research/research_w_slam/` before supplying that public online root. Replay uses the three online inputs; old truth remains confined to the separately sealed offline audit.
2. **Normalization relocation remains incomplete.** `tools/g1_normalization_repair.py:16` directly uses `ROOT.parents[2]/research-w-slam/.../FINITE_GRID.npz`. In a normal B1 checkout this points outside `research/research_w_slam/`. It lacks the replay/offline-audit relocation argument. A documented, explicit old-grid path is needed before claiming this command portable; no source-search fallback or new physical generation is needed.
3. **Two distinct reproduction paths are missing from the current guide.** Saved-result inspection/table-figure regeneration needs restored BGGA JSON and no old GT or physics. Fresh cached replay needs the old online bank, a new output directory, its original scope/accounting conditions, and the explicit relocation flag. Directed normalization uses sealed M4 NPZ plus old prior grid. Do not present these as a single `all` command or as a cold PDE benchmark. Many outputs use exclusive creation; use a fresh reproduction destination and retain the original ledger/receipts.
4. **Import/host requirements.** `prepare_publication.py:4` imports `bgga` without adding `src`; the documented wrapper supplies `src/vendor` (`tools/run_job.py:17–21`), while standalone `python tools/prepare_publication.py` lacks an installed package entry point. There is no BGGA `pyproject.toml`/setup file. Keep invocation explicit or make the helper self-contained. The wrapper also imports Unix `fcntl/resource` and uses `preexec_fn`; native Windows execution has not been established. Its optional local sibling `.deps` preference should be removed only in a documented public copy, as the W-SLAM `PUBLICATION_FILES.json.public_copy_adaptations` convention records, or explained as absent in a clean checkout.
5. **Publication checks validate a local tree, not the public package.** `prepare_publication.py:11–15` checks only root Markdown links with local existence. A locally present but ignored Release-only file passes. It does not check nested reports, generated HTML URLs, archive restoration, or anonymous availability. `build_site.py:9` maps only selected suffixes to GitHub; fragment/query suffixes and other relative targets follow a different path. Images remain relative and therefore must be committed or explicitly hosted. Its heading renderer adds no IDs, so heading-fragment links are not supported by the generated HTML.
6. **Current readback helper still proves the old preflight, not a new package.** `tools/publication_readback.py:22` uses the old `266,326` marker; it checks report, resource request, repository page and Pages HTTP status, with byte equality only for raw report. It has no current G1/figure/CSV/Release membership checks. Saved `PUBLICATION_READBACK.json` is a historical PASS for commit `6fd504785f4e4be74d919f2363e0b2d1ef111534`; this audit performed no new HTTP verification. Use current reviewed content/commit and actual figure/data paths for a later readback, while keeping publication success separate from scientific approval.

Handoff: fix the hardcoded preflight generators and cumulative cost summary first; publish compact Git files plus a manifest of the actual raw Release package; document exact old-bank restoration and normalization relocation; then do the already-authorized publication/readback workflow under the main thread. This audit neither authorizes additional experiments nor changes any scientific gate.
