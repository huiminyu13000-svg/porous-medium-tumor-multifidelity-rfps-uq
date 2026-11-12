# Test1a：BF、RFPS 与 TRF 完整实验

本包保留原 `Test1_a_lambda_bi_vivo.m` 对应的营养不确定性实验。
m=8 为默认，m=64 是同一实验的原参数对照；不包含 Test1b、Test2 或其他实验。

## 运行

在 MATLAB 中将当前目录设为本文件所在目录。LS 求解需要 Image Processing Toolbox。

```matlab
run_experiment                    % Test1a，m=8；重新生成 LF/粗网格和所选 HF 基底
run_experiment('fresh',64)         % 同一 Test1a 的 m=64 对照
run_experiment('fresh_reference',8) % 连 100 个 HF 参考样本也重新计算
run_experiment('cached',8)         % 复用快照，重新执行选点、补点、投影和误差计算
run_experiment('pde_check',8)       % 单个原样本的完整时间 PDE 计算
plot_saved                        % 直接画原 BF 补点误差及 TRF 均值、标准差
```

`fresh` 保留原 100 个昂贵 HF 样本作为比较参考和补点节点值。
`fresh_reference` 是包含参考样本的完整重算。运行输出写到 `results/`，不覆盖输入数据。
直接画图仅保留 Test1a 必要的数值缓存，不包含原大型 workspace 中的图形对象。

## 论文实验：Test 2(a) / Test 2(b)

H. Yu, L. Liu, Y. Feng, T. Qi, *On multi-fidelity methods for a tumor growth model with uncertainties*，§5.2，图 3–9；方法见 Algorithm 3.1–3.2 和 Appendix A。
原代码文件编号 `Test1a` 对应这一个营养不确定性实验，包含 BF 和 TRF 两个部分。

- 随机维数 5，各分量独立服从 U[-1,1]。
- λ=50(1+Σ z_k/(2k))，cB=20(1+Σ z_k/(2k))，G0=0.5(1+z1)，R0=0.45(1+z1)。
- 初值为半径 R0 内密度 0.95 的圆盘；区域 [-2.5,2.5]^2，T=1。
- 细网格 dx=0.05，粗网格 dx=0.1；m=8 时 HF/粗 PME 步长分别为 0.0015/0.006，m=64 时为 0.0005/0.002；LS 步长为 0.03。
- 训练集 N=500；独立统计集 Ns=1200。Stage I checkpoint K0=15，K1=0，Stage II 最多 K2=5 个补点。
- RFPS 权重 ω=0.4；εtol=τtol=1e-3、χtol=2、κtol=3，稳定化参数 1e-12。

Test 2(a) 比较 LS-BF 的 Cholesky 基线与 RFPS 补点；Test 2(b) 使用
CMPME–LS–PME 的分层结构。论文在相同 HF 预算下比较重构精度，
并在独立 MC 样本上展示均值与标准差。误差定义见稿件式 (5.1)。

## 本包运行

本包保留原提供程序及其保存数据的计算入口；输入参数保存在
`data/bf_compute_m08.mat` 和 `data/bf_compute_m64.mat`。
`plot_saved` 直接绘制随原程序提供的保存场和误差；`run_experiment` 执行计算。
`run_smoke` 对两个 m 值在 T=0.05 调用所有 PDE solver，再执行缓存快照上的
BF/RFPS/TRF 重构和直接画图。短时测试不用于生成论文的完整统计结果。

## 文件与依赖

| 文件 | 功能 | 对应原程序 |
|---|---|---|
| `run_experiment.m` | 完整 Test1a 流程 | Test1_a_lambda_bi_vivo / larger_m 与 tri step3 |
| `solver/solve_pme.m` | HF PME | solver_PME_vivo_modified |
| `solver/solve_heleshaw.m` | LF level-set | solver_HeleShaw_vivo_modified |
| `solver/solve_cmpme.m` | 粗网格 PME | solver_PME_vivo_modified_medium_tri |
| `solver/select_basis.m` | 初始基底选点 | cholesky_node_selection_ip_data |
| `solver/enrich_basis.m` | 原 RFPS 补点 | greedy_LFres_radius |
| `plot_saved.m` | 原数据画图总入口 | 调用 plot_error、plot_fields |
| `plot_error.m` / `plot_fields.m` | 误差 / 均值标准差图 | 原画图逻辑 |
| `data/bf_compute_m08.mat` / `bf_compute_m64.mat` | 原参数、候选、LF/HF/粗快照、基底 | 原 Test1a 两个 workspace 与 TRF workspace |
| `data/bf_enrichment_paper.mat` / `trifidelity_paper_fields.mat` | 精选画图数值 | 原保存误差与场 |

营养场、压力场、通量、距离函数等依赖保留在对应 solver 内。
整数时间参数先转为浮点再计算步数，保证原 T=1 对应 500 步 PME；
粗网格插值和细网格算子分别保留。基础函数只简化文件名，不删数学步骤。

## 论文结果与运行检查

论文报告 RFPS 在所测工况下改善同等 HF 预算的 BF 重构；TRF 利用中间层信息进一步改善重构。
本次检查范围与运行结果见 `TESTS.md`，不将短时检查视为完整论文数值实验。

GitHub 可直接保留 `solver/`、`data/` 子目录。将本目录内容作为仓库根目录，
保留相对目录关系即可运行；`results/` 自动生成并由 `.gitignore` 排除。
