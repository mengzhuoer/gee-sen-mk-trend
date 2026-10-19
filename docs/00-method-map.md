# 00 · 方法地图：原文 → 本仓库

本仓库不是对原文的转载，而是把原文涉及的方法逐项拆开、重新实现并验证。
下表用于核对“原文讲了什么、对应本仓库哪里”，也标出了我们做的扩展与修正。

## 1. 原文方法逐项对照

| 原文部分 | 方法点 | 本仓库对应 |
|---|---|---|
| 先说结论 | Sen 给幅度、MK 给显著性，二者互补 | [README](../README.md) 第 1 节 |
| 01 Sen's slope | β = 全部年份对斜率的中位数；β 的符号与量级 | [01-sen-slope.md](01-sen-slope.md) |
| 01 读法 | 单位 = 变量单位/年；NDVI·a⁻¹ | [01-sen-slope.md](01-sen-slope.md) 第 4 节 |
| 02 示例数据 | MOD13A1，500 m，16 天，比例因子 0.0001 | [04-annual-metric.md](04-annual-metric.md) 第 2 节 |
| 02 年度合成 | SummaryQA ≤ 1、NDVI ∈ [-2000,10000]、4—10 月最大值 | [04-annual-metric.md](04-annual-metric.md) 第 2 节 |
| 03 Step 1 | 年度 NDVI 集合；`toInt()` 保留 ties | [05-gee-implementation.md](05-gee-implementation.md) 第 2—3 节 |
| 03 Step 2 | `ee.Reducer.sensSlope()`；年份波段从值波段派生 | [01-sen-slope.md](01-sen-slope.md) 第 5 节 |
| 03 Step 3 | MK 的 S、ties 修正方差、Z | [02-mann-kendall.md](02-mann-kendall.md) 第 2—6 节 |
| 03 Step 4 | 五级分类 | [03-trend-classification.md](03-trend-classification.md) |
| 04 分级表 | 1—5 的判定条件与配色 | [03-trend-classification.md](03-trend-classification.md) 第 1 节 |
| 完整代码 | 地图、图例、年际曲线、面积统计、导出 | [code/gee/sen_mk_ndvi.js](../code/gee/sen_mk_ndvi.js) |
| 完整代码 | `setDefaultProjection`、`no unmask(0)`、每张影像 clip | [05-gee-implementation.md](05-gee-implementation.md) 第 2—3 节 |

## 2. 本仓库新增的内容

| 扩展 | 说明 | 位置 |
|---|---|---|
| GPP 口径改造 | 状态量取最大 → 通量取年累加；46 时段完整才累加 | [04-annual-metric.md](04-annual-metric.md) 第 3 节 |
| 山西案例 | 2001—2022，GAUL 边界，2023 版省界叠加 | [08-case-shanxi-gpp.md](08-case-shanxi-gpp.md) |
| 分阶段 GEE | 检查 / 导出 / 下载 / 指纹 / 补缺，规避容量限制 | [code/gee/gpp_shanxi/](../code/gee/gpp_shanxi/) |
| 本地 Python | 纯 numpy 向量化，7—10 秒完成 22 年趋势 | [06-local-python.md](06-local-python.md) |
| 精确 ties | 本地实现统计包括 0 在内的所有重复值 | [02-mann-kendall.md](02-mann-kendall.md) 第 6 节 |
| 数值自测 | 6 项测试，逐项对照教科书/暴力法 | [tests/test_algorithms.py](../tests/test_algorithms.py) |
| 文件指纹 | 哈希文件名认领年份，独立仓库 | [yaogan-zhiwen](https://github.com/mengzhuoer/yaogan-zhiwen) |
| 踩坑复盘 | 现象 → 根因 → 解决 → 教训 | [07-pitfalls.md](07-pitfalls.md) |
| FAQ | 15 个常见问题 | [09-faq.md](09-faq.md) |

## 3. 与原实现的有意差异

1. **`ALPHA` 参数化**：原文硬编码 1.96 / 0.05；本仓库用 `ALPHA`，便于敏感性分析。
2. **精确 ties**：原文 GEE 数组技巧会把真实 0 值并入哨兵组；
   本地实现精确统计，含 0 数据更准确（见 02 第 6 节）。
3. **完整序列掩膜**：原文已有 `count == N_YEARS`；本仓库把它作为硬性前置条件，
   并在本地脚本入口校验网格。
4. **年度指标分离**：原文在 GEE 内完成年度合成；本仓库把“年度口径”与
   “趋势算法”分开，本地脚本只接收同网格年度栅格。
5. **容量策略**：原文单体脚本；本仓库提供分阶段 GEE + 本地计算两条路。
6. **面积口径**：本地用球面像元面积按行计算；GEE 用 `ee.Image.pixelArea()`。

## 4. 原文未覆盖、但复刻必须知道的点

这些不在教程正文里，却是实际落地时最容易卡住的：

- GAUL 省名是 `Shanxi Sheng`，不是 `Shanxi`；
- GEE Code Editor 的“认证失败”可能只是工作区没绑 Cloud Project；
- `getDownloadURL` 的 `filename` 不生效，下载文件是哈希名；
- 同一批年度文件可能混入不同投影；
- 缺测补 0 会制造虚假趋势；
- 自相关会让 MK 显著性偏高；
- GPP 上升不等于碳汇上升。

逐条见 [07-pitfalls.md](07-pitfalls.md) 与 [08-case-shanxi-gpp.md](08-case-shanxi-gpp.md) 第 6 节。
