# proposal → correction → full validation → archive

## 已执行与原E0的区别

原二维保护格缺陷0；31/32有限编辑未生出声明尺度数字洞，唯一新洞保留却不相容。原有限编辑invalid→valid=0。不能用“填洞”解释全部失败。

G重新生成72个真实field，raw full-valid0，endpoint full-valid18，不相容54；70个raw结构编辑、14个终点保留（含invalid），14个合法终点折返prefix结构邻域。不要把保留而invalid的结构算成功，也不要把材料contrast reward算结构。

144QP：120收敛、24未收敛flag。每个raw/两步/独立端点都保存χ、current、field与32×8复预测。QP未收敛不抹掉；端点合法性独立判断，未建立全局优化保证。

新结构world只有p250073 bandit一次：raw89.9456→25.7212→22.3180<23.78545，独立重算预测差0，solver和c=χE缺陷~1.4e-15；retention=.138013。support42→52，孔数0→0，是边界/支撑变化，不是新孔；dS=.081767，.16敏感性下不计数。两轮同world不算复现。Best-of-K仍由base取得。

新发现不能全归生成器，也不能说corrector毫无作用：这条提案是实际invalid→valid；其余大量提案由有限corrector折返或仍失败。区分物理不可行与32实chart/两步优化不足需要T的共同机制试验，本轮没有实现第二套corrector。

E0质量替换重放有369次替换，但G当前新轨迹中替换未触发；两轮结果相同。修复一个archive缺陷不自动带来新结构或MAE改善。
