# Research publication conventions

Updated by the user on 2026-11-18. The standing preference is to keep research in the existing project repository and report site; explicit task restrictions override it. A direct “publish” authorizes the normal commit, push, Release upload, and existing-site deployment workflow without repeated confirmation.

## Default workflow

1. Publish the completed report, relevant code/evidence, and requested package to the existing project. Preserve scientific limits and negative results. A new phase does not need a new repository.
2. Use the push/upload result and, when the task includes a website deployment, the provider's completion status. Do not call a queued or failed deployment completed. Wait for the necessary status with reasonable intervals, not frequent unchanged polling.
3. Return the public links and any actual failure. No obligatory publication receipt, second metadata commit, or second package is needed merely to mark completion.

## No automatic post-publication test suite

- Do not download the released ZIP again, compare bytes or hashes, run SHA-256/CRC, perform anonymous download tests, crawl links, repeatedly probe Pages/CDNs, take browser screenshots, or rerun scientific/smoke tests just because publication occurred.
- Do such work only when the user explicitly requests it or there is a concrete upload/deployment/correctness problem that the check will resolve. Keep the diagnosis targeted; stop when the issue is resolved.
- Existing archive manifests may be retained as provenance. Do not generate new manifests or hashes by default. Ordinary documentation-only edits need no research experiment or deployment acceptance exercise.
- If publication status is successful, finish the task; do not invent an additional “acceptance gate.” If the provider reports failure, fix that specific failure or report the real blocker.

## Repository and evidence boundaries

- B1 stages stay in research/ with a concise current entry in START_HERE and the existing site. C1 continues in its own existing repository; A17 remains its separate theory/publication project.
- Preserve historical URLs, releases, frozen data, failures, costs, and license status. Keep navigation useful rather than adding a large new process for each round.
- Do not rewrite frozen scientific results to improve a release. Clearly distinguish post-result development from blind evidence, and publication from scientific validation.
- Publish source and evidence within the requested scope. Exclude credentials and irrelevant sensitive runtime dumps if encountered. This does not require rescanning the entire project for every upload.
- Publication does not authorize new physics runs, locked-data access, training, hardware changes, or autonomous research frameworks.
