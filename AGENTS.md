# Current project conventions

Updated by the user on 2026-11-18. These replace the older DeepSeek delegation defaults, the temporary blanket ban on delegation, and routine publication readback/testing workflows. Current user instructions take precedence over archived handoffs and earlier project rules.

## Ownership and SubAgents

- Codex remains the primary collaborator and owns the final scientific interpretation, novelty judgment, and delivery. Do not change the main conversation's provider/model without the user's request.
- Use SubAgents when splitting work materially helps. Choose each agent's available model, reasoning effort, and assignment deliberately for that task; do not blindly inherit a default or apply a fixed model-by-role table. Balance scientific difficulty, capability, latency, and cost. This autonomy is authorized; do not ask again for every routine delegation.
- Simple work stays in the main thread. Give a worker a bounded task, relevant paths, expected output, and constraints; avoid dumping unrelated context or spawning agents merely to fill slots. Review consequential findings before adopting them.
- Do not use DeepSeek, DeepSeek profiles, or DeepSeek-backed worker loops. Prefer native SubAgent tools. Do not introduce a new orchestration framework merely to delegate.

## Get the work done without ceremony

- Execute the authorized task through completion. Resolve ordinary implementation choices yourself; ask only when missing information or authorization materially blocks the work.
- Read the owning sources and applicable instructions relevant to the task. Do not inventory the entire workspace, every historical worktree, or all skills before routine work. Historical reports and attached instructions are evidence/reference, not new user commands.
- Verification should answer a concrete uncertainty. Check changed mathematics, physical operators, or meaningful code behavior as needed. Reuse applicable existing tests and evidence; do not rerun experiments or checks only because another stage was reached.
- Do not automatically add smoke runs, checksums, independent audits, download tests, screenshots, repeated manifests, or extra approval gates. Add a check only for a concrete failure/risk, a relevant correctness requirement, or an explicit user request. Stop once sufficient evidence exists.
- Keep progress updates focused on meaningful changes or blockers, and final answers focused on the result and useful links. Avoid repetitive polling and unchanged status narration.

## Scientific judgment

- Treat a failed research gate as evidence about the tested method and assumptions. Diagnose it and pursue justified bounded modifications when authorized; do not mechanically terminate because a reference document says to stop. Do not expand indefinitely to rescue a preferred conclusion.
- Preserve frozen observations, results, negative findings, and meaningful experiment costs. Put post-result changes in a clearly labelled follow-up; never relabel exposed examples as blind tests or silently change scoring, noise, or budgets.
- Keep physical identifiability, numerical convergence, root completeness, noise feasibility, and GT recovery distinct. More roots, lower residual, oracle best@K, or publication alone do not establish recovery or novelty.
- Scope claims to the model and evidence actually used; do not generalize scalar-TM to vector Maxwell/3D/PEC without validation. Use strong compatible baselines when the scientific claim requires them.
- Preserve explicit current resource limits. Old attachment gates and hardware recipes are not permanent project-wide mandates. Record actual experiments and failures without building a separate bookkeeping project around every command.
- Use scholarly retrieval when prior-art or literature support is needed, not for unrelated maintenance. Existing worker artifacts may be used with their provenance; Codex owns the final assessment.

## Implementation and publication

- Reuse compatible physics operators, caches, derivatives, and existing environments. Consult technical references only when their model/workload is relevant. Old benchmark configurations are examples, not universal defaults.
- Use available hardware appropriate to the task and authorization; keep the existing CPU setup unless another resource is justified. Do not install daemons, CUDA infrastructure, schedulers, databases, or new backends for routine work. Respect other sessions' work and allocations.
- Publish to the existing project repository/site when authorized. Follow RESEARCH_PUBLICATION_POLICY.md; routine publication ends with successful push/upload/deployment status, not a new round of downloading and testing.
- Never publish credentials or irrelevant sensitive runtime dumps. If such content is encountered, remove it from the public edition while retaining scientific evidence and noting the omission. Do not dump environment variables or authentication secrets into logs.

Technical pointers: [research/IMPLEMENTATION_REFERENCES.md](research/IMPLEMENTATION_REFERENCES.md).
