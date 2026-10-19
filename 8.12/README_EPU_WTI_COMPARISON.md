# 美国/中国 EPU + 官方 WTI：模块位置比较

## 四组主规格
|模型|EPU经济含义|方程归属|风险变量递归排序|
|---|---|---|---|
|US_GLOBAL|美国新闻型EPU|独立GL模块|EPU → WTI|
|US_US|美国新闻型EPU|美国模块|EPU → WTI → y → dp → Δr → de → deq|
|CN_GLOBAL|中国SCMP EPU（Factiva）|独立GL模块|EPU → WTI|
|CN_US|中国SCMP EPU（Factiva）|美国模块|EPU → WTI → y → dp → Δr → de → deq|

“中国EPU放进美国模块”仅改变方程归属，不把指标变成美国EPU，也不等于中国模块。
全球规格GL方程没有国家宏观解释变量；各国接收5个外国宏观变量及2个共同风险变量。
美国规格美国7个内生变量接收5个外国宏观变量；其他国家接收5个外国宏观变量及2个共同风险变量。
因此这不是重命名实验，而是宏观反馈和弱外生性限制的敏感性检验。

## 数据（冻结下载日期 2026-10-19 UTC）
- 美国：https://www.policyuncertainty.com/us_monthly.html
  文件 US_Policy_Uncertainty_Data.xlsx；Main News Index；News_Based_Policy_Uncert_Index。
  不使用旧的Three Component Index，也不使用全球GEPU。
- 中国：https://www.policyuncertainty.com/scmp_monthly.html
  SCMP_China_Policy_Uncertainty_Data.xlsx；China EPU Index_Factiva；China News-Based EPU。
  不与ProQuest序列混接，也不混用大陆报纸版。该文件Factiva最后一月为2026-01；完整季度截至2025Q4。
- WTI：美国能源信息署 EIA RWTC，Cushing OK WTI Spot Price FOB，美元/桶。
  https://www.eia.gov/dnav/pet/hist/LeafHandler.ashx?n=PET&s=RWTC&f=M
  https://www.eia.gov/dnav/pet/hist_xls/RWTCm.xls
  使用名义现货价格，不是Brent，不是WTI期货。

原始下载文件保存在 epu_wti_sources/；manifest记录原始URL、工作表、SHA256、覆盖期。
epu_wti_monthly_audit.csv保留月度值；epu_wti_quarterly.csv保留季度均值、月份数和自然对数。
每季度必须完整3个月；先算月值等权算术平均再取log。
EIA月均值的等权季度平均与按每日观测数量加权季度平均并不相同。
中国EPU原始月值中的真实0保留，季度均值为正才取log；无填0、插值、外推、跨版本拼接。

EPU引用：Baker, Bloom and Davis (2016), Measuring Economic Policy Uncertainty, QJE；
中国SCMP另引用 Baker, Bloom, Davis and Wang (2013)，及 PolicyUncertainty.com。
原始文件保留网站署名和使用说明。WTI来源署名 U.S. Energy Information Administration。

## 宏观与估计口径
承接当前 compare-independent-shocks 的GDP对数水平口径，并使用新加的14国利率一阶差分。
其他宏观变量及贸易权重沿用现有流程；历史权重镜像处理仍沿用，不宣称已撤销该处理。
共同样本由同一份宏观CSV决定，四规格均严格覆盖，任何缺失即中止，不独立缩短。
当前预期样本2000Q3–2025Q3，101期；以运行 aligned_input.csv 和 settings.R 为准。
默认 p=1、外国滞后q=1、save=500、burn=500、thin=0.5；固定先验和种子。
初步500/500运行不构成正式收敛性验证。
同一个联合模型分别识别EPU和WTI两种冲击，统一归一化为自身当期+10%。
WTI创新混合供求与信息因素，不能直接称为“石油供给冲击”。
排序固定EPU在WTI前，这是递归识别假设，非数据证明的因果先后。
GDP为log水平，IRF不累计；r为差分，另输出逐抽样累计的利率水平响应。

## 运行
Actions → EPU WTI - Global versus US Block。
首次提交自动运行：数据/结构检验 → 四规格MCMC → 合并比较。
每个模型后验、抽样IRF与失败记录独立保存；重新作图不必重跑MCMC。
preflight验证真实样本的设计矩阵，并用已知AR(1)系统验证两种位置、p=1/2、两类冲击及归一化。

## 输出与解释
- aligned_input.csv / design_audit.csv / settings.R / provenance.txt
- posterior.rda / irf_draws.rda
- irf_summary.csv：valid_all及stable_only，68%与90%区间；
  response_kind=raw为原变量响应，rate_level_cumulative为逐抽样累计Δr。
- stability.csv：总抽样数、有效数、失败数、稳定数、稳定比例、G条件数。
- failed_draws.csv（如有）：记录数值失败，不把失败抽样当稳定；IRF阶段不强行加ridge。
- placement_comparison.csv：同EPU、同日期、同冲击、同变量的中位数差、符号及显著性变化。
- US/CN_EPU_WTI_placement_comparison.pdf：两种模块位置叠加图与90%区间。
- all_irf_summaries.csv / all_stability.csv / model_completion.csv。

不同规格分别估计，不能把同编号抽样当成配对抽样；中位数之差不是差异的可信区间。
稳定抽样筛选改变条件分布；必须结合stable_share和valid_all结果解释。
若某时点没有稳定抽样，不伪造该时点稳定IRF；看稳定性报告。
若任一模型失败，比较任务保存部分结果并显式失败，不宣称四模型齐全。
