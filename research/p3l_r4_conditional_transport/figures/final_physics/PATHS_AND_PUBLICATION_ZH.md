# 发布目录与原运行路径说明

发布应包含本目录：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r4-conditional-transport/research/p3l_r4_conditional_transport/figures/final_physics`。

原生成目录确实位于仓库根的 `figures/final_physics`，绝对路径为 `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r4-conditional-transport/figures/final_physics`；此前与阶段目录 OUT 的区别未说清楚。该原目录及所有已封账输出完整保留。本次使用 `copytree(copy_function=copy2)` 复制全部文件，并逐文件直接比较字节、文件大小和纳秒修改时间；没有重新生成图、读取新数据、重新拟合或调用物理求解器。副本核验明细见 COPY_PROVENANCE.json。

原生成代码在仓库根 `tools/plot_r4_final_physics.py`，排版修正代码为 `tools/refine_r4_final_physics_layout.py`。原代码的 DEST 指向仓库根图目录，这是实际运行路径，因此没有改写。运行时所用的两份代码副本也随本目录的 `source/` 完整复制，可用于发布源码与图证据对应。

本目录内 FACTS.json、VISUAL_QA.json、COST_PROFILE.json 和 source 代码中的绝对路径继续保留原运行、原来源或原封账收据路径；这些字段是出处，不是对本次副本存放位置的声明。不能用搜索替换改写已封账证据。发布图像请使用本目录四组 PNG/PDF；图说明用 DRAFT_ZH.md；紧凑数据用同目录 CSV。原 Stage IV 绘图账为 `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r4-conditional-transport/research/p3l_r4_conditional_transport/IV/jobs/final_physics_figures_v1/CLOSED.json`，本次复制费用独立记录，不重计原物理成本。

发布工具若只收集阶段输出，应以 `research/p3l_r4_conditional_transport/figures/final_physics` 为入口，无需重跑生成器。
