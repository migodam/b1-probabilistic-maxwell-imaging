# 复现与实际执行边界

从现有B1仓库根运行。需要已有NumPy、SciPy、PyTorch、Matplotlib、Pillow和threadpoolctl；本轮使用原Python3.13环境、M4/MPS网络与单线程CPU complex128完整物理。没有安装CUDA或新调度系统，没有占用远端4060。权重和数组恢复见[ARTIFACTS.md](ARTIFACTS.md)。

所有下列命令是可复算入口；已发布结果由原回执产生，不表示打包时再次重复运行。不要在冻结目录覆盖结果，复算请在另一个克隆中执行。

## 依赖与数据

保留仓库research目录层级。R5复用R4的`data/operators.npz`、`data/manifests/evaluation.json`、`audit/OPTIONAL_FT_CONTRACT.json`、TRAIN/VAL材料与观测、原严格`anchors/`。固定谱来自`research/research_gn_green_error_20261106/results/fixed_spectra.npz`。R5归档含新`inputs/noiseaware_anchors/`以及实际11套权重／全部预测与候选；不重复打包整个旧研究。

R4公开数据、离线评分材料、严格GN与固定谱分别在其已发布Release的data-public、data-evaluator、anchors、fixed-projectors附件；原冻结Python依赖附件包含所需相对导入和包初始化器。精确下载命令见归档页。离线GT只用于TRAIN监督及封存后的开发评分，公共预测不读取它。

## 从头训练：实际12轮上限

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python research/research_gn_nn_feedback_wavelet_r5/experiments/train_development.py --anchor-mode strict --epochs 12 --patience 6 --batch 8 --device mps --seed 7
python research/research_gn_nn_feedback_wavelet_r5/experiments/prepare_noiseaware.py
python research/research_gn_nn_feedback_wavelet_r5/experiments/train_development.py --anchor-mode noiseaware --epochs 12 --patience 6 --batch 8 --device mps --seed 7
python research/research_gn_nn_feedback_wavelet_r5/experiments/train_development.py --anchor-mode noiseaware --arms haar_closed_coupled --noise-ball --epochs 12 --patience 6 --batch 8 --device mps --seed 7
```

最后一条在独立`haar_closed_noiseball_coupled`目录输出，不覆盖前五主臂。命令默认20轮／patience8是一般入口；**本轮实际均显式传12／6**，由各模型config/summary记录。CPU可用`--device cpu`；没有对其他平台承诺本机计时。

## 公共GN输入预测

```bash
python research/research_gn_nn_feedback_wavelet_r5/experiments/predict_public.py --checkpoint research/research_gn_nn_feedback_wavelet_r5/models_noiseaware/haar_closed_coupled/best.pt --operators research/research_gn_nn_preservation_r4/data/operators.npz --anchors PUBLIC_ANCHOR_1.npz PUBLIC_ANCHOR_2.npz --output new-run/public_predictions.npz --steps 4 --device mps
```

输入只含公开m0/j0/Q0/Y、算子与噪声，不接收类别、幅度、字体或GT。学习网络内部迭代无完整LU；输出材料仍用原F验收。观测→GN→CNN→物理校正的统一包装入口：

```bash
python research/research_gn_nn_feedback_wavelet_r5/experiments/infer_observation.py --observation PUBLIC_OBSERVATION.npz --operators research/research_gn_nn_preservation_r4/data/operators.npz --checkpoint research/research_gn_nn_feedback_wavelet_r5/models_noiseaware/haar_closed_coupled/best.pt --output new-run/inversion --device mps
```

这份统一包装代码已实现；本轮实际运行了各组成环节和批次入口，**没有额外重新运行一个未知观测的单体包装测试**。观测文件须有Y和公开noise。它会依据checkpoint的anchor_mode选择GN停止规则，准备公共参考并记录冷启动费用。

## 校正、评分、绘图

训练入口只训练／选择checkpoint，不会自动写公共预测。恢复归档已有完整预测；若从头重训，先按以下相同VAL公共锚点执行批量预测，再评分。这里仅读清单的split／obs_id，不读GT数组：

```python
from pathlib import Path
import json
from research.research_gn_nn_feedback_wavelet_r5.experiments.predict_public import infer
root = Path('research/research_gn_nn_feedback_wavelet_r5')
r4 = root.parent/'research_gn_nn_preservation_r4'
records = json.loads((r4/'data/manifests/evaluation.json').read_text())['records']
ids = [r['obs_id'] for r in records if r['split'] == 'val']
arms = ['pixel_closed_coupled', 'haar_closed_coupled', 'db2_closed_coupled',
        'haar_closed_no_fullfield', 'haar_free_coupled']
for mode in ['strict', 'noiseaware']:
    models = root/('models' if mode == 'strict' else 'models_noiseaware')
    anchors = r4/'anchors' if mode == 'strict' else root/'inputs/noiseaware_anchors'
    results = root/'results' if mode == 'strict' else root/'results/noiseaware'
    chosen = arms if mode == 'strict' else arms+['haar_closed_noiseball_coupled']
    for arm in chosen:
        infer(models/arm/'best.pt', [anchors/(i+'.npz') for i in ids],
              r4/'data/operators.npz', results/'VAL_PREDICTIONS'/(arm+'.npz'),
              steps=4, batch_size=8, device='mps')
```

```bash
python research/research_gn_nn_feedback_wavelet_r5/experiments/score_development.py --anchor-mode strict
python research/research_gn_nn_feedback_wavelet_r5/experiments/score_development.py --anchor-mode noiseaware
python research/research_gn_nn_feedback_wavelet_r5/experiments/evaluate_corrector.py --phase strict --workers 2 --execute
python research/research_gn_nn_feedback_wavelet_r5/experiments/evaluate_corrector.py --phase noiseaware --anchor-root research/research_gn_nn_feedback_wavelet_r5/inputs/noiseaware_anchors --workers 2 --execute
python research/research_gn_nn_feedback_wavelet_r5/experiments/evaluate_corrector.py --phase noiseaware --anchor-root research/research_gn_nn_feedback_wavelet_r5/inputs/noiseaware_anchors --arms haar_closed_noiseball_coupled --tag noiseball --workers 2 --execute
python research/research_gn_nn_feedback_wavelet_r5/experiments/score_noiseball.py
python research/research_gn_nn_feedback_wavelet_r5/experiments/collect_costs.py
python research/research_gn_nn_feedback_wavelet_r5/experiments/summarize_development.py
python research/research_gn_nn_feedback_wavelet_r5/experiments/plot_development.py --phase strict --require-results
python research/research_gn_nn_feedback_wavelet_r5/experiments/plot_development.py --phase noiseaware --require-results
```

校正入口默认只列计划；`--execute`才运行。先保存公共候选／排名，之后离线读VAL GT。原封存回执包含原机绝对路径，作为历史证据保留；移机复算在新目录生成新回执，不能编辑旧记录冒充原运行。

直接恢复归档即可重画图，不需要为重画重复GN或训练。额外noiseball模型结果在单独文件中，五主臂表和图不偷偷加入这个后续对照。`results/COST_LEDGER.csv`将实际新执行和独立部署参考分开；相同上限不是相同实际LU/RHS。

新鲜TEST、多种子和高对比度为NOT_RUN。可运行的代码不能替代未运行的科学证据。
