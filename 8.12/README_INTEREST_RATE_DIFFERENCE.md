# 利率一阶差分版

本版本以现有 compare-independent-shocks.yml 的数据口径为基础：GDP 为重建的对数水平，CPI、REER、股价保持现有处理，仅改变 14 个经济体的利率。原始 Excel 和现有流程保留。

## 处理规则
- 对新构建的 RATE_LEVEL 映射列 *_r，计算 r[t] - r[t-1]；不取对数、不计算增长率、不乘除 100。
- 单位沿用源利率。若源数据以百分数记录，差分单位是百分点。
- 必须季度连续、升序、无重复，数值完整；异常即终止。
- 删除第一期；保留期间的所有非利率变量逐字段不变。
- 外国利率仍由固定贸易权重加权；现在加权的是利率差分。
- 只向 data/rate_diff 写入新文件；不可把差分输出再作为输入。始终从现有构建器重新生成水平输入。
- 这不会把原始 Excel 的 RATE_LEVEL 列改名或覆盖。

## 数据下载
Actions → Prepare First-Differenced Interest Rates → 最新成功运行 → 下载 interest-rate-first-difference-*。
压缩包包含：
1. model_input_rate_diff.csv：可直接用于模型的利率差分数据；
2. model_input_rate_level_aligned.csv：同一样本期的利率水平对照；
3. rate_difference_audit.csv：各国每期前值、现值、差分；
4. RATE_TRANSFORM.json：实际起止期、删除期、样本数及校验结果。

## 重新估计
Actions → Compare Independent Shocks - First-Differenced Rates → Run workflow。
此流程保留原四组独立 GPR/VIX/EPU/Oil 模型、滞后阶数选择、先验及识别方式，读取差分数据，输出到 results/independent_shocks_rate_diff。
输出中内部变量名仍为 r，但经济含义是 Δr；现有通用图程序若显示“利率水平”，应改成“利率变动”。
利率差分 IRF 不能直接解释为利率水平 IRF。需逐个后验抽样从 h=0 累计后再求可信区间，不能简单累计区间上下界。
所有模型必须重新估计；已有后验及旧 IRF 不适用于新数据。
比较差分和水平规格时，用同样本对照文件重估水平版；不要把样本缩短的影响混入规格比较。

## 本地运行
从仓库根目录执行：
```bash
Rscript --vanilla 8.12/02b_build_14e_input_gdp_loglevel.R
python3 8.12/difference_interest_rates.py
```
GDP 构建器需要沿用现有流程的 TVPGVAR_GDP_DLOG_SCALE=100 等环境配置。新流程已显式设置。
