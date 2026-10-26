# GPR各类别：分别估计外生路径与自身峰值 IRF

2026-10-26，从 [Caldara–Iacoviello官网](https://www.matteoiacoviello.com/gpr.htm) 下载原始月度工作簿，官网标注月度版本更新于2026-10-01，数据到2026-09。原文件、SHA256、原始字段说明、季度化输入、完整峰值清单均冻结在 `data_20261026/`。

## 运行范围：14个独立模型

| 字段 | 含义 | 原始数据范围 |
|---|---|---|
| GPR / GPRT / GPRA | 总指数 / 威胁 / 行动 | Recent，10家报纸 |
| GPRH / GPRHT / GPRHA | 总指数 / 威胁 / 行动 | Historical，3家报纸 |
| SHAREH_CAT_1 | War Threats，战争威胁 | Historical文章占比 |
| SHAREH_CAT_2 | Peace Threats，和平威胁 | Historical文章占比 |
| SHAREH_CAT_3 | Military Buildups，军事扩张 | Historical文章占比 |
| SHAREH_CAT_4 | Nuclear Threats，核威胁 | Historical文章占比 |
| SHAREH_CAT_5 | Terror Threats，恐怖威胁 | Historical文章占比 |
| SHAREH_CAT_6 | Beginning of War，战争爆发 | Historical文章占比 |
| SHAREH_CAT_7 | Escalation of War，战争升级 | Historical文章占比 |
| SHAREH_CAT_8 | Terror Acts，恐怖行动 | Historical文章占比 |

这8个细分列是3家报纸历史口径下的**文章占比**，不能误称为10家报纸的近期指数，也不能当作互不重叠、相加恰等于总指数的分解。Recent与Historical的基准期和报纸覆盖不同，分别估计，不拼接。

完整工作簿还有44国各两种国家指标、6个检索定义变体、文章总数和重复量纲的总文章占比。它们在 `official_series_inventory.csv` 全部列明，本轮未把这些当成额外的事件类别。没有估计文章数N10/N3H，没有凭空造出不存在的新版细分类。

## 外生冲击的具体含义

每次仅选一个类别 `z_t = ln(GPR类别的完整季度均值)`。国家模块为：

`x_i,t = c_i,t + A_i,t x_i,t-1 + B_i,t x*_i,t + D_i,t x*_i,t-1 + beta_i,t z_t + gamma_i,t z_t-1 + u_i,t`

- 14国的内生变量始终仅为五类宏观变量，共70个。GPR和其滞后是每个国家VARX的外部控制。
- 为计算未来IRF路径，另估计一个不含任何国家宏观变量、外国变量和贸易权重的自主标量TVP-AR(1)：`z_t = a_t + phi_t z_t-1 + eta_t`。
- 程序用独立GL路径行进行矩阵堆叠。它不是新增的国家内生GPR方程。GL行对任何宏观变量的当期和滞后系数都强制为0，并在每个IRF抽样中复核。
- 注入外部路径的创新，国家创新设为0，归一化到原指标当期上升10%，即对数变量冲击 `ln(1.10)`。
- 本轮没有同时加入其他GPR类别、EPU、WTI或VIX。不同类别的结果不解释为控制其他类别后的净效应。

**外生路径是施加的模型限制，不等于已通过外部工具变量或准实验识别出因果效应。** 条件OLS残差检查也不是外生性检验。

## 样本与处理

保持上一轮宏观候选快照，以免同时改变宏观数据和GPR：

- 14国：AU BR CA CH CN EA UK JP KR NO SG TR US ZA。
- 每国：y=100Δln(GDP)、dp=100Δln(CPI)、r=Δ利率、de=100Δln(REER)、deq=100Δln(股价)。
- p=1，2000Q1仅用作初始滞后，实际估计为2000Q2–2025Q3，共102季度。
- 新GPR历史保留到2026Q3，但不能把宏观共同样本之后的峰值当作已估计状态下的IRF时点。
- 月度0值是真实观测，保留后再计算季度均值。必须三个月齐全，不能用两个月平均。实际使用的全部季度均值均大于0，因此取自然对数，不加epsilon、不把月度0改成缺失。
- 四条利率代理定义差异（CN/US/SG/TR）及原贸易权重来源问题依然存在，随每个模型的结果输出，不假装已统一或解决。

## 自身峰值

所有类别使用同一预先固定规则，但**对各自季度原始水平单独筛选**：局部高点、达到本组75%分位、四季度邻域突出度至少本组0.5标准差、峰间至少四季度。平台取首季度，较高峰优先，等高取较早。端点不强行视为双侧峰。

只在实际估计窗口筛选。共123个类别—峰值组合，各组数量见 `run_manifest.csv`，日期与数值见 `peak_dates.csv`。不固定每5年选一个，不按IRF显著性反选日期，也不把所有类别统一强制放在2003/2008/2020/2022。

图题标出类别及自身峰值季度。对应历史事件的因果解释尚未核验，文件明确标记，不自动把同季度事件称为峰值原因。

## 抽样、诊断与输出

首轮默认每组 burn=1000、保存阶段=1000、thin=0.5，预计实际保留500个后验抽样。相较上一轮250抽样有所加长，但仍是**每组一条独立链**，不宣称已通过多链收敛。

1. 先验证冻结数据、0值处理、每组自身峰值、15个模块的实际设计满秩。
2. 验证解析解下的外部AR传播与禁止宏观反馈，实际标量采样器做短抽样功能测试。
3. 分别估计各组，保存原始后验后才计算IRF。
4. 输出选定参数的单链ESS、条件OLS残差BG筛查、每个峰值的谱半径、稳定比例、G条件数及失败抽样日志。
5. 同时保留 `valid_all` 和 `stable_only`。后者是稳定条件后验，不能代替前者或称为MCMC修复。
6. 每图固定一个类别、一个自身峰值、一个响应变量，14国分面。深色68%、浅色90%，不叠加不同冲击。
7. 利率变动与利率水平响应分别输出。水平响应先在**每次后验抽样内**累计，再求可信区间，不累加区间边界。

每个模型的附件：

| 文件 | 内容 |
|---|---|
| posterior.rda | 原始模型后验，供重新作图/诊断 |
| irf_draws.rda | 全部选定日期的逐抽样IRF和稳定掩码 |
| irf_summary.csv | 中位数、68%及90%上下限 |
| IRF_68_90.pdf、figures/ | 全部自身峰值、五变量及利率水平图 |
| significance_counts.csv | 各时点、变量、区间、后验筛选的正/负/含零数量与国家数 |
| country_horizon_significance.csv | 各国显著和不显著的具体预测期 |
| stability.csv、posterior_ESS_screen.csv | 数值稳定性与抽样精度筛查 |
| settings.R、LIMITATIONS.txt | 设定、解释与尚未解决的问题 |

国×预测期单元数每组182（14×13），国家数每组14，分开统计。可信区间为逐点贝叶斯区间，未进行多重比较修正。

工作流一次最多并行4组，其余排队。正式估计由提交信息 `[run-gpr-categories]` 或手动运行明确触发，普通代码检查不会重跑14组。每组结果单独上传，另有轻量FIGURES附件，避免只看图也必须下载大后验。

## 来源与复现

Caldara, Dario and Matteo Iacoviello (2022), *Measuring Geopolitical Risk*, American Economic Review, 112(4), 1194–1225。数据下载于2026-10-26，原网页标注CC BY许可。官网与数据URL见 `source_manifest.json`。

```bash
pip install pandas numpy xlrd openpyxl scipy statsmodels matplotlib
python3 8.12/gpr_categories_exogenous/build_data.py --verify
python3 8.12/gpr_categories_exogenous/test_data.py
python3 8.12/gpr_categories_exogenous/diagnostics.py
Rscript 8.12/gpr_categories_exogenous/test_architecture.R
# 安装原采样器依赖后：
GPR_SERIES=SHAREH_CAT_1 Rscript 8.12/gpr_categories_exogenous/run.R
python3 8.12/gpr_categories_exogenous/plot.py --series SHAREH_CAT_1
```
