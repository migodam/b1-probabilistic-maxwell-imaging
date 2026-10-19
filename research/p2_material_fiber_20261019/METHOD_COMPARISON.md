# 同物理目标、完整费用比较

| 方法 | 合格held-out数 | 不同替代/对象 | Best-of-4 MAE | 该panel基础MAE |
| --- | --- | --- | --- | --- |
| B:full_J | 4 | 0.00 | 0.048404 | 0.048404 |
| B:hard_prior | 4 | 1.75 | 0.046260 | 0.048404 |
| B:native | 4 | 0.00 | 0.048404 | 0.048404 |
| B:random_matched | 4 | 0.00 | 0.048404 | 0.048404 |
| B:soft_Maxwell | 4 | 0.00 | 0.048404 | 0.048404 |
| B:soft_adaptive | 4 | 1.50 | 0.045905 | 0.048404 |
| B:soft_dynamic | 4 | 1.50 | 0.045928 | 0.048404 |
| B:soft_frozen | 4 | 1.25 | 0.047154 | 0.048404 |
| D:flow | 4 | 0.50 | 0.048404 | 0.048404 |
| D:gaussian | 4 | 0.50 | 0.048404 | 0.048404 |
| D:nonlearned | 4 | 5.50 | 0.046163 | 0.048404 |

主局部表为4个合格held-out，不能隐去另外4个基础恢复失败。soft_dynamic相对hard的Best-of-4 MAE差−.000333，95%配对区间约[−.002487,.002063]，未确立优势。hard 136 RHS，dynamic504 RHS，adaptive440 RHS（相同已执行held-out blocks）。

## 补充强基线

Full-J sphere-only发现held-out10个替代，hard-prior sphere-only9个。该审查发生在主评分后，明确为exposed；原角向progress约束可能限制Full-J，不能据主对照零成功证明其能力不存在。

## A17与标准评价器

| 评价器 | query+setup RHS | fallback /48 | cold中位ms | warm中位ms | setup摊销ms |
| --- | --- | --- | --- | --- | --- |
| A17_exchange | 296 | 17 | 38.660 | 47.572 | 59.490 |
| Krylov | 88 | 7 | 12.307 | 24.522 | 23.563 |
| exact_T | 1056 | 0 | 7.639 | 7.700 | 9.487 |
| full | 384 | 0 | 12.441 | 1.622 | 5.346 |
| goal_ROM | 192 | 4 | 13.635 | 28.911 | 38.249 |

所有setup纳入，full允许自有factor重用，普通ROM允许增量basis，固定材料不被corrector改变。

E fallback采用已付费offline reference误差审查，非可部署在线误差界；耗时表标为oracle-assisted benchmark，不作为在线速度保证。
