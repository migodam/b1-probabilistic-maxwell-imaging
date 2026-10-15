# Research A — 交付与执行入口

**Maxwell Compensation-Driven Multi-Solution Imaging**  
研究日期：2026-10-14  
当前裁决：**GO for bounded falsification；尚无新的成像优势或强 TAP 成果主张。**

## 先读什么

| 文件 | 用途 |
|---|---|
| `RESEARCH_REPORT_ZH.md` | 完整中文研究报告：证据审计、连续/离散 Maxwell、条件命题、八架构、主算法、替代路线、MAE/选择下界、文献矩阵、阶段预算、13问裁决 |
| `CODEX_EXECUTION_PROMPT_ZH.md` | 直接交给 Codex 的下一轮实现与执行任务。默认仅 E0/E1；二维/三维/训练均有独立禁用开关 |
| `MaterialHypothesisPackage.schema.json` | Route A 到 Route B 的机器可读候选交接格式 |
| `validate_package.py` | 检查 JSON Schema、跨字段一致性、候选权重语义和可选资产 SHA-256 |
| `validation/verify_algebra.py` | 可独立运行的八矢量偶极代数检查 |
| `validation/results.json` | 上述检查本轮实际结果 |
| `validation/test_contract.py` | 纯合成的接口正确/错误案例检查 |
| `validation/contract_test_results.json` | 接口测试的实际结果 |
| `validation/environment.json` | 本轮运行环境 |
| `MANIFEST.sha256` | 交付文件完整性清单，不含清单自身 |

## 本次实际做了什么

已核对报告中列出的项目文件与相关文献；对无法取得全文的文献明确记录访问层级。推导和比较了材料可实现性、局部流形、有限相对反馈、局部编辑 Schur 核、受限补偿夹逼、选择风险下界以及条件集合风险。

实际执行的数值工作仅为八矢量偶极 CPU 小矩阵检查和合成接口测试。没有新 Mie 成像试验、二维/三维重建 campaign、GPU 工作、网络训练或新 response bank；没有改写、推送任何上游仓库。

**当前包不包含一套新产生的真实 `MaterialHypothesisPackage` 候选数据。** Schema 与 validator 是下一轮的实现契约，不是已经生成并验证了新材料世界的证据。`validation/test_contract.py` 中的合成 fixture 不能向 Route B 当作真实观测交接。

## 复现小型检查

Python 3.10+；NumPy、SciPy、jsonschema。依赖见 `requirements.txt`。在本目录执行：

```bash
python validation/verify_algebra.py
python validation/test_contract.py
```

这些命令覆盖对应测试结果 JSON，不接触外部数据。保留交付结果时，应先复制整个目录到新的工作副本。

验证将来生成的候选包：

```bash
python validate_package.py /path/to/manifest.json
python validate_package.py /path/to/manifest.json --check-assets
```

`--check-assets` 以 manifest 所在目录为资产根目录，逐文件核对 SHA-256，拒绝逃逸路径。本版按**文件**校验资产；Zarr 目录应先打包成文件，或将其逐文件清单实现为独立扩展，不能把目录 hash 假装成普通文件 hash。

Validator 验证的是格式与已声明字段的一致性：它不会重新求解 Maxwell，不证明数值误差界、网格收敛、被动性检查本身的正确性或全部材料覆盖。正式科学验收仍须调用冻结的物理 backend。

## 关键阅读提醒

最推荐的算法是 material-first finite edits + nonlinear compensation，而不是 raw TSOM current-nullspace sampling。TSOM、ROM、MDN 均可被删除而不破坏主研究问题。

原始 2009/2011/2015 TSOM 全文访问有缺口；部分论文只取得原始摘要或书目元数据；PDF 截图失败处未用于视觉事实判断。报告不宣称穷尽文献或已经证明“首次”。2026 年 PNGF 的 full-wave 局部低秩更新是必须正面比较的计算先例；MCL、null-space exploration、nonlinear data-consistent learning 与集合反演也不是本任务的发明。

所有“拟议 GO”“预算”“合格研究 MAE”均是下一轮规范，不是本轮已达到的结果。工程合格 MAE 仍需应用容差；不可将研究性相对改善当作工程认证。

## Codex 使用方式

将整个目录放在能访问既有 B1/C1/A17–A30 资产的工作区，把 `CODEX_EXECUTION_PROMPT_ZH.md` 作为任务输入，并让其先读研究报告。建议先用 **Plan mode** 核对仓库路径、资产可用性、budget wrapper 和 truth isolation，再只执行默认 E0/E1。既有 541 次剩余额度不会被继承，本包也未授权自动进入大型三维确认试验。
