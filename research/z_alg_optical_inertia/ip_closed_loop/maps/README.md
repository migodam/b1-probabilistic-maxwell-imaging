# 固定 7 发 9 收的全视场材料定位

本轮从单频 400 MHz、边长 2 m 的二维离散 TM 观测出发，在每个整数格中心检查一个物理材料类。32² 和 96² 分别覆盖 1,024 和 9,216 个中心，包含四边和四角。候选以格中心平移，超出成像域的图案部分截去。数据在决策时只给原始复数端口、噪声球半径、已校准的离散算子和材料类。空间图不依赖重建图或真值掩码。一个物理见证只标记通过独立全波验算的具体材料及其中心。

[读者图及结果](index.html) · [逐对象 CSV](results/centres_summary.csv) · [四类逐对象 CSV](results/library_per_object.csv) · [宽环逐对象 CSV](results/broad_centres.csv) · [机器结果](results/summary.json) · [先验和生成样本冻结](freeze.json) · [公开观测索引](public/acquisition.json)

## 观测和材料类

- 标量 TM，时间约定 $e^{-i\omega t}$，等面积盘脉冲单元。两网格均为 $2.667\lambda$ 的正方形。32² 采用与存档 96² 相同的发收角和发射归一，接收端口为未加权物理坐标，重新构造相应单元积分。96² 原样复用固定算子论文实验 B 的 Fourier 符号、七束平面波和九行接收算子。发射跨 135°，接收跨 146.25°。两网格用保存的二进制系数定义精确的离散合同。跨网格比较先对齐原生接收归一，随后单独审查 [真实材料对的网格转移](model_transfer/README.md)。
- 新生成 16 件独立 32² 对象，实盘、环、双夹杂、Austria 型各四件，另有两件原生 96²（环和右边界裁剪双夹杂）。每件有一对物理背景：外部真空与实际逐格扰动背景。后者在每一外部格独立生成于 $[0,1/32]$，推断允许区间内的所有逐格实值。每次观测用相同相对范数 3% 的复噪声规范，保护半径比投入噪声大 $2\cdot10^{-6}\|Y_{\rm clean}\|_F$。同一对象的置信域、噪声及见证回放是配对条件。
- 主要形态合同中，局部轮廓 $p_i$ 在核心取 1，边缘取 0.5。每一活动格独立满足 $3.75p_i\le \chi_i\le 4.25p_i$。所有外围和孔内均逐格允许 $[0,b]$。四形态库条件把四种允许类的并集赋给每个中心，四类库同时使用没有形态与中心标签的代号观测，四批证书保存后再解封对照表。宽环合同在核心使用 $[3,5]$，可选壳层使用 $[0,5]$，其他格仍独立允许 $[0,b]$。
- 真值及生成电流在 [offline](offline/) 保存，在线脚本加载不含生成标签的 [材料和观测合同](contract.json) 及 [public](public/) 所列数组。评分阶段读取两边。形态方案、独立中心、散射幅度、随机种子、窗口规则先写入冻结文件。物理输出可按 [复现命令](#复现命令) 重生。

## 证书机制

在周期填充网格上，给原域之外的传感算子列补零。任一原问题的物理电流以零延拓。外围扇区 $|j|^2-b\Re(\bar j e)\le0$ 允许该零延拓，因此在更大的空间排除也是对原域的保证。统一外围权重取 1，固定状态块及其逆为

$$
K_0=I-b\operatorname{Re}_{H}G,\qquad
K_0^{-1}=\mathscr F^{-1}\operatorname{diag}\big[(1-b\Re\lambda(G))^{-1}\big]\mathscr F .
$$

平移只改变活动格的本构行 $T=[P_A;P_A G]$。以 $D$ 表示局部扇区相对外围扇区的差，$\Sigma=TK_0^{-1}T^*$，$W=D^{-1}+\Sigma$。每个不同的剪裁图案复用自己的 $W$ 及经验证的惯性。相同图案的不同平移共享它们。在线以公共的 Fourier 协方差提取各位置对应的源、接收耦合。

$$
J_b=K_0^{-1}{b B\over2},\quad V_0=TJ_b+[0;B_A],\quad C=SK_0^{-1}T^*.
$$

用 $W^{-1}$ 消去局部更新后，接收块为 $P=SK_0^{-1}S^*-CW^{-1}C^*$，中心为 $Y_c=SJ_b-CW^{-1}V_0$，源块半径记 $\Theta=\tfrac b2 B^*J_b-V_0^*W^{-1}V_0$。对任一正权 $\tau$，带数据的状态负指数为 $\nu=n_+(W)-n_+(D)-n_-(I+P/\tau)$。当这些逆或符号不能验证时，当前权重不作决定。其余计算在 9 端口、7 源的小矩阵完成。

$$
H=(Y-Y_c)^*(I+P/\tau)^{-1}(Y-Y_c)-\tau\Theta .
$$

取 $H$ 的最大 $\nu$ 个特征值之后的正特征值之和作为平方残差下界。严格超过 $r^2$ 时排除。该关系保留了非零负指数。常值正权扫过的局部更新可用这里的 $D^{-1}$ 公式，遇到未验证的奇异局部更新保留未决。同轮完整的非奇异和奇异固定迹算法在 [reduction](../reduction/README.md) 单独给出。

数值验算使用 128 位 Arb 逆 Fourier 缓存、Arb 检验复指数根的直接 DFT 和有保守累加界的线性代数。令 $\widehat X$ 为近似解 $W X=V$。用
$$
V^*W^{-1}V=
V^*\widehat X+\widehat X^*V-\widehat X^*W\widehat X+
(V-W\widehat X)^*W^{-1}(V-W\widehat X)
$$
控制源、接收和交叉协方差的误差。原先用整个逆的范数误差放大所有协方差，32² 实盘弱背景每件仅排除 15–32 个中心。保留原公式并改用上述残差二次项后，每件为 1017–1021。独立的 [Arb 回放](replay.py) 从原始符号重建协方差，并直接用区间线性求解和惯性检验核查抽取的结论。

## 对照、见证和适用范围

同合同方向性 Neumann 使用局部中值的全波解、实方向一阶梯度及其支持函数，独立传播活动区和全域背景的非线性余项。另记录旋转阻抗正性得到的不同逆范数界。OSM 由同一复观测和归一后的采样原子直接相关，只报告采样峰及真中心位次。宽环的 PSD 对照保留相同待查中心、相同五个相对权重与十七个比例，要求带数据的状态块正定。更宽连续权重搜索与独立任意精度复查发现至少一个该有限权重库漏掉的 PSD 排除，详见 [PSD 复查](../reduction/psd_maps/README.md)。比较结论限于实际运行的权重和预算。

局部见证从端口匹配的中值全波模板提出，允许随后用局部最小二乘调整，必要时允许在合法弱背景区间内全域调整。每个已验收材料用已知弱参考解析地约化全波方程，通过参考收缩、局部 Schur 逆缺陷和全域 VIE 残差核验 $\|SJ-Y\|_F<r$。弱背景见证可显式是非均匀的。材料区间和剪裁检查与物理端口复算均在结果中保留。

在一个经定位的中心，幅度分档检验 $p_i[a,b]$ 的公共档位并集。一般逐格变化的单点峰值是不同的问题：对整盒加全局上界给出峰值下界，对“任一格达到阈值”逐格查并，需要全部子盒排除才得到上界。条件幅度分档缩窄不会替代这一峰值检验。依据单个高值物理见证可以标记相容，其他未有证据的单点保持未决。

合成数据用于检验带噪离散模型内的集合判定。完整视场、四形态库与未知弱背景这些结论适用于列明的材料先验和仪器模型，独立几何库外形态、跨网格误差、损耗、三维和校准偏差需要另设合同。

## 复现命令

从项目根执行，使用带 [requirements.txt](requirements.txt) 的 Python。在已有公开数据上可从算子界直接重放。

```bash
M=research/z_alg_optical_inertia/ip_closed_loop/maps
python -B "$M/prepare.py"
python -B "$M/operator_bounds.py"
python -B "$M/gallery.py" prepare

for shape in solid ring double austria; do
  python -B "$M/certified_map.py" 32 "$shape" primary 3.75 4.25
  python -B "$M/sector_map.py" 32 "$shape" primary 3.75 4.25
  python -B "$M/sector_map.py" 32 "$shape" library_all 3.75 4.25 --blind --fast
  python -B "$M/witness.py" 32 "$shape"
  python -B "$M/contrast.py" 32 "$shape"
done

python -B "$M/gallery.py" reveal
python -B "$M/sector_map.py" 32 ring broad 3 5 --broad
python -B "$M/sector_map.py" 32 ring broad_psd 3 5 --broad --psd
python -B "$M/replay.py" 32 ring M32-R01 broad --samples 3

for shape in ring double; do
  python -B "$M/certified_map.py" 96 "$shape" primary 3.75 4.25
  python -B "$M/sector_map.py" 96 "$shape" primary 3.75 4.25 --fast
  python -B "$M/witness.py" 96 "$shape"
  python -B "$M/contrast.py" 96 "$shape"
done
python -B "$M/sector_map.py" 96 ring broad 3 5 --broad --fast
python -B "$M/replay.py" 96 ring M96-R01 primary
python -B "$M/replay.py" 96 double M96-P01 primary
python -B "$M/replay.py" 96 ring M96-R01 broad
python -B "$M/replay.py" 96 ring M96-R01 broad --interior
python -B "$M/sector_map.py" 96 ring broad_psd_interior 3 5 --broad --psd --fast --interior
python -B "$M/model_transfer/run.py"

python -B "$M/score.py"
python -B "$M/stress.py" --prepare
for shape in solid ring double austria; do python -B "$M/stress.py" 32 "$shape"; done
for shape in ring double; do python -B "$M/stress.py" 96 "$shape" --fresh-native-witness; done
python -B "$M/plot.py"
python -B "$M/build_report.py"
```
