# 来源与依赖

本目录的新增研究实现使用Python、NumPy、SciPy、Matplotlib、scikit-fem和Triangle。依赖由包管理器安装，不在仓库或数组Release中复制其安装目录。各库保持原有许可证；尤其不把Triangle重新声明为MIT许可。复现或商业使用时请检查其上游许可。

`vendor/p0/`是本仓库旧Research W P0的源文件快照，来源commit `db7146cae2d74ba231577d42f225dbfc9a389268`，保留原研究来源与本仓库许可状态。没有导入未经审核的A17/A26接口，也没有将旧四块恒等式说成新增理论。

原Pro研究包、52页报告及一维脚本不重复打包，见[旧P0来源与审计](https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022)。本轮的输入文件列表在`docs/INPUTS.json`。Git版本号用于标识来源；未执行额外SHA256校验。

本目录不是对所有第三方材料的重新授权。参见仓库根目录`LICENSE_STATUS.md`。
