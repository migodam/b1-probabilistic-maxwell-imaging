# R2 资源审计

审计状态：**PASS_INTERNAL_ACCOUNTING**；发现 0 个账本/回执不一致。

快照范围：sequence 866738，最后账本时间戳 1792693656.32776；读取时活动预留 `{'final_resource_audit': 1800.0}`。
本次自身分析作业的结束时间及 JOB_END 尚未进入此快照；本报告不是最终封账。

| stage | reserved RHS | state attempts | reserved LU | controlled wall s | physics wall s | GPU s |
|---|---:|---:|---:|---:|---:|---:|
| A | 2109800 | 86961 | 86961 | 6277.660459583392 | 941.7836741455831 | 0.0 |
| B | 1334472 | 150169 | 150169 | 1769.7762113341596 | 1259.243868453661 | 0.0 |
| C | 311800 | 7795 | 7795 | 1981.660588498693 | 88.35886803478934 | 0.0 |
| D | 147456 | 18432 | 18432 | 191.92124658473767 | 120.29703528922983 | 0.0 |
| L | 2052800 | 59360 | 59360 | 1390.600012539886 | 585.3946993213613 | 0.0 |

RHS 角色、后端实际 forward/LU/receiver/gradient 计数及其覆盖作业数见 RESULTS.json。预算预留和后端完成操作分别列出，不把两者当成同一个计数；JOB_END 与 CLOSED 重复的 counters 只累计一次。

失败作业数：2；未消耗 forward credit 回执逐作业保留。

## 每父样本上下文实测成本

### first256

| component | RHS total | RHS median | RHS P95 | wall total s | wall median s | wall P95 s |
|---|---:|---:|---:|---:|---:|---:|
| first256/native | 565840.0 | 1976.0 | 3720.0 | 331.1115254559554 | 1.2275867915013805 | 2.178482770686969 |
| first256/probit | 454368.0 | 1844.0 | 2368.0 | 279.9021674140822 | 1.1256213129963726 | 1.4682468539103866 |
| first256/final_context_J | 0.0 | 0.0 | 0.0 | 0.01984279160387814 | 7.091660518199205e-05 | 0.00012125051580369473 |
| first256/matched_extra | 55792.0 | 128.0 | 536.0 | 75.98713745758869 | 0.24998195841908455 | 0.48665096901822835 |

first256 最终 probit 上下文的描述性 noise-ball：256 / 256 个有 Phi 回执的父样本；缺失 0。阈值 282.8737921622988。

### second256_confirmation

| component | RHS total | RHS median | RHS P95 | wall total s | wall median s | wall P95 s |
|---|---:|---:|---:|---:|---:|---:|
| second256_confirmation/native | 529864.0 | 1668.0 | 3712.0 | 298.71108254091814 | 0.9916648539947346 | 1.983830270881299 |
| second256_confirmation/probit | 432600.0 | 1704.0 | 2368.0 | 256.9790967425797 | 1.02893385454081 | 1.3655763230635785 |
| second256_confirmation/final_context_J | 0.0 | 0.0 | 0.0 | 0.02085983077995479 | 6.866652984172106e-05 | 0.00011713505955412984 |
| second256_confirmation/matched_extra | NOT_RUN | NOT_RUN | NOT_RUN | NOT_RUN | NOT_RUN | NOT_RUN |

second256_confirmation 最终 probit 上下文的描述性 noise-ball：252 / 256 个有 Phi 回执的父样本；缺失 0。阈值 282.8737921622988。

第二批为 native96 + probit64；matched16 为 NOT_RUN（设计未执行），不计为缺失失败。
两批成本与 stationarity 分开报告；native KKT 与 probit 梯度分别统计。优化器未成功、预算到限或缺失回执均保留，不能宣称所有点已达到 MAP 收敛。

## 测量边界

已付背景数组导入、缓存读取、setup 与 CPU 训练/分析的受控墙钟由所属账本作业包含。已有付费缓存不是从零免费建立；每轮训练的 CPU/wall 只作分解，不重复加到总账。
GPU 累计：0.0 秒；这是快照值，不是未结束作业的最终资源证明。GPU 本行增量缺字段的范围见 JSON。
单物理 worker 与 receiver transfer cache≤4 是源码策略；本工具只审计活动预留，不能据此证明操作系统层面没有并发。
成本是原执行中保存的增量测量，不包含新冷启动/热启动重复实验；首父样本可能承担背景 setup，后续复用缓存。未记录的 cold/warm 分拆标为不可用，不能由总量反推。
只读 R2 ledger/CLOSED、上下文 costs/stationarity JSON 与训练 RESULTS；未打开 GT、观测数组、LU 或模型，未重复实验。

缺字段/缺文件类别：0；详见 RESULTS.json。

## 交付范围

上表是分析审计时快照。最终追加式账本和打包作业结束后的封账记录作为Release独立资产交付；最终数值以 `CAMPAIGN_CLOSED.json` 为准。原始流水不删除或重置，打包/上传不增加任何Maxwell RHS。历史P3-L的8,374,472 RHS保持封账，不属于本轮。
