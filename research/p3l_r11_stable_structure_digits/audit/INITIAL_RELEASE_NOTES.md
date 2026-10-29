R11-S3：三个原生 Sol/ultra 并行，完成N/P两套各九个数字、四个后端共72次真实非线性Maxwell反演。含原始truth、测量、预测、全部费用、失败和盲视觉记录。

- A精确活动支撑：稀疏数字完整Jacobian静态核观察节时63.54%，整成像优势尚未普遍建立。
- B固定总场：原192节点anchor128仅14通过；数字trial全部精确回退。
- C统一耗散界：b*=0.009210948658>0，但未带来加速。
- 数字图像并未全部恢复：8双孔失败；Full66 posterior与材料多解没有认证。
- **资源协议FAIL_CAP_EXCEEDED：28,727全阶LU超过自设12,000；363,240保守RHS低于400,000。发现后停止新增物理，不重定义预算。**

完整中文报告：https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/p3l_r11_stable_structure_digits/

完整tar.zst为无损原始证据；reading kit只作轻量阅读，不能冒充全部数组。manifest提供BLAKE2b-256和字节数。完整包按仓库相对路径解包到隔离目录；含可移植的零新增Maxwell分析入口。历史R10路径原输入见原p3lr10-bilateral-20261030 Release，本轮未覆盖历史报告。
