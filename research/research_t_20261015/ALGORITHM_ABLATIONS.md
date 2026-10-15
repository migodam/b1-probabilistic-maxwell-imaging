# 实际消融

1. 同N3、同物理的full与Haar64：数据阈值两臂通过，但full MAE更差；这是材料表示/正则自由度的重要负结果。
2. 40个完全相同编辑的无补偿版本均失败；普通nonlinear补偿有32条通过，不能归因于只做一个免费几何编辑。
3. v1与v2相同raw和三步预算，只增加非core正材料相保护。额外新增near-air均值29.44→11.06，而通过15→10；非core自由度参与补偿。
4. v3更换minimum-change base，其余corrector不变，固定8个提案7个通过。它保留现象但不是独立盲测，也没有提高材料MAE。
5. constrained screen只是付费轨迹顺序重放；A17独立增量、G的adaptive/QD、真正3D均NOT_RUN，不能从合成代数通过推断这些增量。
