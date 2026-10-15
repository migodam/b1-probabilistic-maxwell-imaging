# QD、bandit与强基础搜索

| 方法（首轮） | 新states | 复RHS | 合法端点 | 材料新颖events | 结构新颖events | Best-of-4复MAE | 旧定义Coverage@4 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| multistart | 36 | 864 | 2/9 | 0 | 0 | 0.04224802 | 33.3% |
| static_scale | 36 | 864 | 2/9 | 1 | 0 | 0.04224802 | 33.3% |
| QD | 36 | 864 | 2/9 | 2 | 0 | 0.04224802 | 33.3% |
| bandit | 36 | 864 | 3/9 | 1 | 1 | 0.04224802 | 33.3% |

表为首轮，第二轮相同。每方法9proposal/36state/864复RHS，公共3base/24RHS全局记一次。所有方法同合法支持/动作库/原32实Haar corrector/full validator；static分连续尺度与幅度，不分孔数。multistart为base-centered相关全域扰动，不冒称全材料域prior。

bandit描述性命中一个结构，QD没有；不是稳定超越static/multistart的证据。仅3个已曝光合格parent、每arm每parent3步，T-guided未运行，五臂不完整。三步beta更新不足以判断学习机制。

Best-of-2/4/8各方法全相同，Top1都是base；原Coverage@K=1/3，几何coverage=1/3。没有用G新归一化coverage=3/3替换旧目标。短K不填充。完整逐parent、prefix、实际RHS/wall见E1_DISCOVERY_SUMMARY和E1_OFFLINE CSV。

当前不能称QD/CEM/RL比强baseline更好，亦不能由小pilot证明它们永远无效。优先接入真正共同T handoff并扩大有效轨迹证据，未获新授权不启动下一campaign。
