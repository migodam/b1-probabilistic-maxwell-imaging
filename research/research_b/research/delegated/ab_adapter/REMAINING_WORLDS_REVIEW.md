# Remaining worlds 导出说明

代码 `research/research_b/code/export_remaining_worlds.py` 已运行，输出 `results/REMAINING_WORLDS.npz` 和 `.json`。不调用 solver/evaluate、不新抽噪声、不修改 E1/E2 原结果。

完整保留 625 材料 catalogue，另存固定3.2、robust、material-risk 三策略 ×32 parent×2 noise 的兼容 mask、全部逐world profile χ²和合法gain。原始材料不因mask为false而从交付删除。

兼容规则仅为逐world bounded shared-gain joint profile χ²≤best+9。这是 heuristic data-compatibility set；不是posterior、置信区域、覆盖保证或完整支持保证。不同world可以有不同gain，但每个world跨整个自身history必须同一个gain。所有概率none、outside mass null，family/event数不能当概率。

Simulator只读取已保存outcomes/calibration；原始old1.6只出现一次，新频率/重复测量按原step index取得同一个potential outcome。筛选函数只接收全625预测、观测和sigma，不接收真实材料、真实gain或offline MAE。三策略使用不同action history，不能把策略比较当新独立parent。

768个已保存trajectory step逐一重放，最小profile χ²与原记录最大绝对差为0。最终192个episode-policy中160个仍保留多个posthoc结构family；这只是有限catalogue内的已发现歧义，不排除漏掉其他世界。每条噪声重复仍属于同一parent。

读取示例：

```python
with np.load("REMAINING_WORLDS.npz", allow_pickle=False) as q:
    worlds = q["catalogue_latent"]       # (625,4), full original catalogue
    mask = q["compatible_mask"][0,0,0]  # policy,parent-index,noise-index
    gains = q["per_world_shared_gain"][0,0,0]
    count = int(mask.sum())
```

精确材料可用已交付 `src/b1_joint/material.py:decode` 从latent重建，支持R=1同心径向材料。原full response bank `E1_RESPONSE_BANK.npz` 是同Release的另一个证据资产。JSON保存policy轴顺序、posthoc event/family定义、原观测ID和block长度；history数组需要按valid_real_channels切片，padding不是观测。

计费在 `Cost.job('e2-export-remaining-worlds-v1','E2')` 中完成，具体CPU/wall见同目录 `REMAINING_WORLDS_REPLAY_AUDIT.json`。
