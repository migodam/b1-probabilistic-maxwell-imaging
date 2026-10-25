# 来源与依赖许可边界

`src/wslam/fem_backend.py`来自B1研究仓库提交`483d260c6b5b111db9fb814b096685fccb362448`中`research/research_w_continuation/src/wq/fem.py`的薄适配，改动及范围见`research/delegated/assets/ADAPTER.md`。历史联合见证来源提交为`958a478739d9553269316ded8125fda9df225c0c`。`legacy_inputs/p0a/analytic_witness.py`保留旧Research W一维Maxwell控制来源，不是新增二维/三维理论或当前宏网格先验。

本目录没有断言一份新项目许可。用户原始研究包、历史源摘录及各第三方材料保留原许可状态；公开可读本身不授予新的复用许可。参阅现有B1仓库根目录`LICENSE_STATUS.md`及相关上游来源。

NumPy、SciPy、Matplotlib、scikit-fem与Triangle作为外部安装依赖，不随此目录复制其安装分发。各库保持原有许可，尤其不得把Triangle自行重标为MIT。可选研究包PDF构建需要外部Pandoc、XeLaTeX和CJK字体，也不在此声明中被重新授权。
