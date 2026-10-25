# 训练是否真的推进、Spline 与 CFM 怎样比较、费用由谁承担

本页从已封存训练/评价表解释结果，没有新模型、推断或 Maxwell 查询。原始来源为 `AUDIT/independent_probability_review/CLOSED_FITS.csv`、`CLOSED_FAMILY_SUMMARY.csv` 和 `AUDIT/final_probability/FRAME_AND_DATA_INFO.csv`。

## 不是遇到 update0 就关闭神经路线

|实际训练条件|验证指标|选定 update|解释|
|---|---:|---:|---|
|Gaussian，4096 TRAIN，lr2e−4|NLL93.9817|0|首轮确实无泛化增量，保留失败|
|同模型与数据，lr8e−4|NLL93.9817|0|提高学习率没有救回|
|同模型与数据，lr5e−5|NLL92.0594|300|较低学习率出现有限训练收益|
|Gaussian，32768 TRAIN、三seed|NLL83.0943–83.2619|1200|数据与优化预算增加后能学习，不是全家族update0|
|Spline，4096 TRAIN|NLL87.9311|300|有限数据仍有增量，但明显较弱|
|Spline，16384 TRAIN|NLL81.5435|4500|同一family的真实规模实验|
|Spline，32768 TRAIN、三seed|NLL77.6174–77.9154|5700–6000|更充分的联合监督泛化|

这里的数据/优化条件并非完整全因子设计。Gaussian 的4k与32k还改变了更新预算，因此不能把全部差异唯一归于数据规模。Spline三档是相同family/LR与相同6000更新上限的实际比较；训练样本重复次数不同，最佳checkpoint仍由独立validation选择。

先前R2“15次都选update0”没有被解释为神经概率方法不可能；R4实测反驳了这种过度概括。另一方面，validation NLL变好也不代表原噪声尺度的完整posterior已正确。

## 相同原始测量输入下，Spline 与 CFM 的可比部分

下表均为32,768 TRAIN、512维原始测量输入、约0.683M参数、三个seed、固定K128/NFE32模型动作。不是更昂贵anchor-conditioned CFM与raw-y Spline混比。

|两批新parent的描述均值|Spline|raw-y CFM|
|---|---:|---:|
|512：material energy score|0.103851|0.111351|
|256：material energy score|0.111692|0.120202|
|512：模型样本均值材料MSE|0.022780|0.026389|
|256：模型样本均值材料MSE|0.026549|0.030960|

这些可比描述指标偏向Spline，但这不是新预注册的显著性检验。CFM的数值NLL未合格，不能虚构一个统一NLL冠军。另有anchor-conditioned有界物理CFM，对同容量随机控制的energy-score改善约21%；它使用更少的1024付费physics contexts，必须与上述raw32k线路分开解释。

## 有限样本均值与模型偏差

原P1材料有界，使直接K128模型样本均值的额外总体MC风险项至多0.002840625。观测到的样本均值与好的CNN风险差约0.016–0.019，不能主要靠增加K解释或修复。严格条件和非逐对象含义见[附加数学诊断](THEORY_CALIBRATION_ERROR_BOUND.md)。没有在测试后换K重新报告主要成绩。

## 费用与可比较范围

- 全部风险动作共享同一份实际付费的native96+probit64基础反演与所需receiver-J缓存；provider/CNN融合没有额外Maxwell查询，不能因此省略共同基础解费用。
- strong CNN使用1024/256 context，历史R2只有512/128；不同训练数据成本单独记，不把同函数类别重新训练的收益叫新架构。
- fusion系统除6418参数head外，还需约0.68M provider、K32生成/decoder与特征写盘。主要测试固定canonical batch64；三个seed是loss平均，不是免费集成。
- 旧32768/2048联合模拟与R2模型是真实历史复用；它们的冷生成不被重计为R4新物理，也不被声称为从零部署免费。新增context及所有新对象physics计入R4独立账本。
- 风险45方法整面板、provider特征、density26面板的实际作业profile保存。没有单模型隔离时延时，不把面板总时间简单除45冒充每模型时延；也没有证明端到端风险系统比单独R2 CNN更快。
- 正面的计算节省只主张固定核的安静Full/interval/Krylov重放，包含setup与fallback。它不覆盖网络训练、概率混合或更大二维系统。

完整拥有成本以追加式科学ledger和每作业CLOSED为准，打包/上传另记delivery成本，不参与算法性能排名。未执行的冷启动部署或分布外扩展仍为NOT_RUN。
