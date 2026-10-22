# 四方法见证搜索

逆算法读取测量、已知背景、噪声、完整48×32材料像素和固定数学任务，不读取wall_values、真实mask或孔洞数。初值为在线background Born/adjoint解及公开prior；所有源/频率共享同一contrast与conductivity。投影统一满足box、passivity和最多3面积，精确空气允许组件生灭。像素/平滑先验仍限制可表达的最小尺度，未建立任意拓扑搜索覆盖保证。

|参考seed / 带宽|Full-J|随机匹配|双初始化FWI|Task-aware|
|---|---:|---:|---:|---:|
|2026102201 / narrow|29.942024|19.912259|32.065442|26.713549|
|2026102201 / wide|45.354386|45.354386|45.354386|45.354386|
|2026102202 / wide|49.892531|48.699200|49.892531|47.129469|


共同正确器用完整离散Full-J、同一ridge/trust幅度和最多4个标准半步trial；接受只依据cal全波失配和同prior。N3、truth最优步、未知墙体Green均没有作为输入。multistart在本pilot仅为共同BP与一个投影到box/面积的随机扰动初始化，不是假称充分多启动。

最小特征缺口：原0.1限制未被该表示/投影强制。本文的“合法”仅在明确的box、损耗、面积和共享频率律层面核验，不能读成全部冻结几何先验获认证。新一轮必须先实现该条件；本轮不事后改变材料类，也不重新搜索。

任务方向来自归一化材料坐标下K与J的正则广义比，只是local proposal。它没有稳定超过简单控制，也没有达到非线性数据相容。此轮不能判断局部危险方向能否沿真正的近等价纤维有限延拓。

原始[修正版轨迹](results/p0b/search-corrected/RESULTS.json)、[排除的错误版本](results/p0b/search/RESULTS.json)、[离线truth评分](results/p0b/SCORED_RESULTS.json)。后者在在线输出保存后运行，不参与选择。表中不同带宽属于同两parent，不是三独立对象；没有伪造bootstrap显著性或posterior概率。

费用表的method wall只覆盖增量proposal/corrector；共享setup、BP、F/J、任务方向与数学求解另包含在整个作业wall/CPU。RHS standalone归因为普通控制+80共享cal RHS、task-aware+160共享cal/task RHS。实际共享费用只在总账计一次；不能拿增量wall当端到端加速。
