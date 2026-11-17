# R5 measurement 复核摘要

对象：`research/delegated/z_alg_r2_publication/repository/research/research_z_alg_r5_measurement`
（只读复核；推理已冻结，允许读 GT 仅用于核对；未运行/未改动任何实验或历史文件。）

## 结论

- 证据链闭合：公式/实现自查、盲测、预算、证书与评分全部通过；无伪造、无结果级错误、无科学过度声明。
- 区间证书的序列化盒处理**正确**：qualification 用“向内一 ULP 安全框 + 向外 K + 严格包含”重建证明链，
  16/16 认证根复算通过，26/26 菜单响应包围成立。正式域是 `certified_safe_binary64_box`，
  原始 Q.json 的 X 端点只是最近舍入、不可单独当证书域（LOW，纯序列化建议：向外加宽落盘）。
- 局部 vs 全局覆盖区分正确：全局/噪声可行集覆盖 12/12 全为 UNKNOWN，无越界声明；
  `material_boxes_disjoint` 是 safe 盒级不交，足以证明两根不同（CERTIFIED_AMBIGUOUS 成立），
  但不等于完整认证盒不相交（LOW 措辞）。
- 被拒匹配方向已显式审计：feed>4 拒绝并落盘；T0 的拒绝偏边界（6.4 vs 4），强超指向只在 T3/T4（109–721）。
- 模型失配真实超过 allowance 2e-5 于 4/6 目标（如实报为 old_model_bound_ok=False）；
  噪声臂设计优势基本消失（maxmin≈fixed），direct completion（3 测量）系统压过单测量设计。

## Gate2 判定

**不支持。** strict coverage 不成立（全部 UNKNOWN），且 noise-margin 测量优势未在 direct completion
与 multistart 归属下存活。建议 NOT_RUN、NO-GO for expansion，与本包 THEORY_DERIVATIONS.md §7/§8 一致。

## 严重度

- 无 HIGH。
- MEDIUM（科学限制，非缺陷）：覆盖未证；allowance 被真实失配超过；噪声臂无测量优势。
- LOW：score_v2 冻结后一行修复（已记录）；X 序列化非向外；safe 盒级不交措辞；噪声相位哈希；
  T0 匹配方向拒绝边界；Q 无显式 30s 墙钟。
- 仅记账/序列化类可改项：X 端点向外加宽落盘；盒级不交改用向外盒测试；score_v2 sha 并入 freeze。

详细证据与行号见 `REVIEW.md`。
