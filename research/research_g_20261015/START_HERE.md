# Research G：Generative Maxwell Solution Search

[中文网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_g_20261015/) · [完整报告](RESEARCH_G_REPORT_ZH.md) · [原始Markdown](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_g_20261015/RESEARCH_G_REPORT_ZH.md) · [全部证据](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-g-20261015) · [实际图库](MULTISOLUTION_GALLERY/index.html)

已执行零物理E0与历史同一A30后端的小型E1，291/384states、6936/12288复RHS、GPU0。三个观测的material代表3/3/1、structure代表1/2/1，新增结构只有一个近门槛边界world；两轮重复不加算。全部方法Best-of-4=.04224802，未改善base；p250077真洞仍未覆盖。

判断为NO_GO_SEARCH：不扩大当前算法优越性主张，不训练PPO/GRPO，不进入新3D。保留合法witness、全部失败及未来T接口入口。T-guided尚未在G运行；本轮不是完整五臂T/G比较。

先读RESEARCH_G_REPORT_ZH.md首屏八问，再看E0_REPORT_ZH.md、SEARCH_HISTORY_SIGNAL.md、QD_VS_BASELINES.md和CORRECTOR_ATTRIBUTION.md。GATE_DECISION仅科学判断；NEXT_STAGE_PLAN中更大factorial/RL/3D均NOT_RUN。

原始新候选引用在CANDIDATE_MATERIALS/index.json，72份侧车在GSearchEpisode/。21个valid endpoint ID包含跨轮重复，不能写成21个不同world。全部full预测和内部fields/current为真实已执行数组，candidate-only、概率质量null。
