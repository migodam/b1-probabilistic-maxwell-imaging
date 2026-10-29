# 公开交付与版本边界

本轮复用现有B1，不创建新仓库。

- [可读中文报告](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/p3l_r11_stable_structure_digits/)
- [匿名原始Markdown](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/p3l_r11_stable_structure_digits/RESEARCH_REPORT_ZH.md)
- [原始数组、代码、完整物理账本与阅读包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p3lr11-stable-digits-20261030)
- [158条公式](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/maxwell_formula_atlas/)

执行数据固定在7d4d1fa2930afacdd8f683c3025e5538b4f3be8a。初版完整tar.zst与reading kit未覆盖，两个payload均经匿名全字节流式BLAKE2b回读匹配，见PUBLIC_READBACK.json；581个科学文件的本地无损归档验证见audit/PACKAGE_ROUNDTRIP.json。

视觉资格补充单独发布R11_report_qualification_update.zip与R11_REPORT_UPDATE_MANIFEST.json，明确P3本轮decoder构造中横缺失。该修订不改任何GT、y、预测、物理solver、原数值fit门、主MSE或费用。169份数组/solver/配置身份保持不变，见audit/REPORT_UPDATE_DATA_IMMUTABILITY.json。

Source_manifest的files列表指初版固定执行身份，报告修订另由manifest固定。不要将旧归档数值PASS误当作九个P对象全部视觉通过。整轮资源协议仍为FAIL_CAP_EXCEEDED。

历史R10路径核的原输入由[原R10 Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p3lr10-bilateral-20261030)提供；R11完整包包含本轮全部输出，以及可独立分析数字数据的固定算子与原prior输入，不复制整个历史TRAIN。
