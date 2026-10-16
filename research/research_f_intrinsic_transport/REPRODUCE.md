# 复现与数据

源码冻结commit：938be85331aa36fa70954337b7e307f2e3cdad59（实际F1执行）；F0审计后增加了缺失说明和已存数组分析，不重算物理。所有路径相对B1仓库。

原后端复用相邻A30仓库717d3d26a290248bd2b9aeb55f3d03694b407b84及A26生成器2603f3bd16c48e48587fedc4363975779def17f5。原physics、constitutive、noise与source/receiver保持不变。资料包内文件实际名称为pro_source/report/RESEARCH_REPORT_ZH.md与CODEX_NEXT_PROMPT.md，并无单独名为Research_F_Codex_Prompt.md的文件；按内容对应读取。

F0无需任何物理求解：既有环境Gaussian/.venv_nn/bin/python运行 tools/research_f_geometry.py；需要上轮Release中的CHART/JACOBIAN/inputs，提供--cache-root，或本地相邻b1-finite-feedback。输出覆盖时应指向新分析目录以保留本轮证据。

F1执行命令曾为：
```
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 ../../.venv_nn/bin/python -m research_f.anchor
```
本轮RUN_CLOSED与stage-start ledger会拒绝无授权重复执行。不要删除账本绕过。要复现付费反演，需另行授权的新复现目录/账本，不能叫零物理重放。

证据Release包含全部24状态的chi/current/field/prediction/LU/pivots、24完整Jacobian、在线观测、离线GT、原已支付背景数组与相关代码。只读统计/绘图不使用LU或调用solver。每endpoint先fsync落盘，后续才计算Jacobian/下一proposal。

原ChatGPT对话“非线性材料解集探索”已取得一个含用户研究任务和Pro总结的完整返回页；hasMore=false。这是服务返回范围，不声称取得其他会话全史。未公开原私人对话，只有所附用户研究包及本轮来源说明。

既有代数脚本在独立副本运行；不是新增2D/3D物理实验。不做SHA256复核，不安装新运行时，不启动GPU训练。当前本地环境SciPy/NumPy已存在。
