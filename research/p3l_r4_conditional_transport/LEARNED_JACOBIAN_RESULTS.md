# Learned-Jacobian与物理动作：真实训练后的裁决

共享几何teacher取得459,712RHS，含10,063个独特配对；公共physics anchors另计。不是每seed重新付同一teacher，复用不等于冷启动免费。

六次局部前向/JVP surrogate拟合（value/jvp×3seeds）都运行1,000updates。原最大白化L2≤.05与max|ΔPhi|≤.05门全部失败；实际最大白化误差约328–348，Phi误差约9.6e4–1.1e5。快速网络推理不能抵消必要Full回退。此处仅是小型局部DINO-like控制，不是完整LazyDINO文献复现。

三个动作误差排序器各6,000updates，误差约.00845–.00895，弱于普通total-field更新+exact decoder约.00724；frozen-J约.07474、Broyden约.02286。排序器还需要候选Full场特征8RHS，不是免费预测。

普通物理复用在材料action fidelity上有价值，但未证明提高最终风险、混合或可靠likelihood决定。停止扩大完整J回归，保留exact decoder、可核验残差与决定驱动Full refresh作为更简单基线。

完整算法、每seed选择update、训练与算子测试门和费用来源见AUDIT/PHYSICS_MECHANISM_FACTS.json与physics_mechanism_tables。原预发布summary可能写“等待locked”；最终已执行结果以本报告与results/IV/interval_chain_locked_v1为准，不改写中期封存文本。
