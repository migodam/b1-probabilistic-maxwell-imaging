# 原2019 CIE公式转录纠正（新3D运行前）

主线程2026-11-06直接核对[开放原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8320915/)的数学公式，发现本轮初版经典案例审计漏写两个平方根：

- Eq.(16)：`delta = sqrt(mean_source(||gamma_new-gamma_old||² / ||gamma_old||²))`，阈值为1e−3。
- Eq.(17)：`Err = sqrt(mean_voxel(|epsilon_est-epsilon_true|² / |epsilon_true|²))`。

最初审计文字将它们写成平方量，主线程第一次分派原尺寸3D时也曾继承该错误。已在该扩展的任何反演执行之前发送纠正：原尺寸实验须使用带根号的停止量；小尺寸固定240循环实验只记录此量，不用它提前停止。漏根号若不纠正，会把停止要求放宽，不能称公平复现。

此事件属于文献公式转录／指令错误，未发生使用错误停止量的R9反演，物理求解费用为0。R8冻结证据不改写，R9相关审计直接纠正并保留本记录。R9的contrast幅度NRMSE、relative-contrast L2和原论文Err分别命名，不混分母。

原文Eq.(15)的source-normalized数据项和本构项相对权重为1；原尺寸论文条件适配采用rho=1，无额外材料prior。小尺寸R8继承合同rho=10另列，不称作者原参数。
