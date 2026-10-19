# 06 · 本地 Python 实现：容量不够时的可靠退路

## 1. 为什么把趋势计算搬到本地

山西案例在 GEE 上一次请求里同时做 22 年累加、231 个年份对、全省统计和出图，
触发 `Earth Engine capacity exceeded`。把“取数”和“算趋势”拆开后：

- GEE 只负责数据筛选、年度合成、导出（可分年、可重试）；
- 本地负责 Sen/MK/分类/统计/出图；
- 22 年 × 1,326,180 像元，全套 7—10 秒；
- 结果完全可离线复现，不依赖 GEE 当时的负载。

## 2. 输入要求

- 每年一个栅格文件，文件名含年份，例如 `GPP_raw_annual_2001.tif`；
- 所有年份**同一 CRS、同一 transform、同一尺寸、同一 nodata**；
- 每年一个波段；
- 缺测用 nodata 表示，不能补 0；
- 年份连续，按升序传入。

脚本会逐文件校验网格，任一不一致立即退出，不会“带病叠栈”。

## 3. 运行

```bash
pip install numpy scipy rasterio matplotlib
# 可选：叠加行政边界
pip install geopandas

python code/python/sen_mk_local.py \
  --data-dir "D:/Zcode/GEE/导出" \
  --pattern "Shanxi_GAUL_GPP_raw_annual_{year}.tif" \
  --start-year 2001 --end-year 2022 \
  --scale 0.1 \
  --boundary "data/shanxi.geojson" \
  --out results
```

参数说明：

| 参数 | 含义 |
|---|---|
| `--pattern` | 文件名模板，必须含 `{year}` |
| `--scale` | 原始值 → 物理值比例（GPP 用 0.1，已是物理值用 1.0） |
| `--alpha` | 显著性水平，默认 0.05 |
| `--chunk` | 每块像元数，默认 60,000 |
| `--prefix` | 输出前缀，多变量时避免覆盖 |
| `--no-plots` | 只算不画图 |

## 4. 计算流程

```text
读取 22 年栅格
  → 校验网格
  → nodata → NaN
  → 完整序列掩膜（22 年全有效）
  → 分块：
       Sen 斜率（231 对的中位数）
       MK S
       ties 精确修正
       Z（连续性修正）
       p（正态近似）
       五级分类
  → 写出 4 个栅格
  → 球面像元面积统计
  → 年际均值 CSV
  → 三张图
```

## 5. 关键实现

### 5.1 完整序列掩膜

```python
flat = stack.reshape(n_years, -1)
complete = ~np.isnan(flat).any(axis=0)
cols = np.flatnonzero(complete)
```

只对 `cols` 计算，省内存也避免对无效像元做无意义的中位数。

### 5.2 Sen 斜率

```python
ii, jj = np.triu_indices(n_years, 1)
dt = (jj - ii).astype(float)
slopes = (x[jj] - x[ii]) / dt[:, None]
slope = np.median(slopes, axis=0)
```

`x` 形状 `(n_years, k)`，`slopes` 形状 `(231, k)`。n=22 时 231 对。

### 5.3 MK 与 ties

```python
S = np.sign(x[jj] - x[ii]).sum(axis=0)

eq = (x[:, None, :] == x[None, :, :])   # (n, n, k)
m = eq.sum(axis=1)
tie_sum = ((m - 1) * (2 * m + 5)).sum(axis=0)

base = n * (n - 1) * (2 * n + 5) / 18.0
var_s = np.maximum(base - tie_sum / 18.0, 0.0)
```

这个写法**把真实 0 值也当作普通重复值统计**，比 GEE 数组技巧更精确。
`tests/test_algorithms.py` 在含 0 的小整数序列上与教科书公式逐一比对。

### 5.4 Z 与 p

```python
safe = np.where(var_s > 0, np.sqrt(var_s), 1.0)
z = np.where(var_s > 0,
             np.where(S > 0, S - 1,
                      np.where(S < 0, S + 1, 0.0)) / safe,
             0.0)
p = np.where(var_s > 0, 2.0 * norm.sf(np.abs(z)), 1.0)
```

`VarS=0` 的常量序列：$Z=0$、$p=1$、分类 3。

### 5.5 面积统计

经纬度网格按行计算球面像元面积：

```python
RE = 6371007.181
dlon = np.deg2rad(abs(transform.a))
lat_edges = transform.f + np.arange(height + 1) * transform.e
area_row = RE**2 * dlon * (sin(lat_hi) - sin(lat_lo))
```

投影网格退化为仿射行列式 `abs(a*e - b*d)`，脚本会打印提示。
不要用“像元数 × 常数”，否则高纬度会有百分比级误差。

## 6. 输出

| 文件 | 内容 |
|---|---|
| `SenSlope.tif` | Sen 斜率（物理单位/年），float32，NaN=nodata |
| `MK_Z.tif` | MK 标准化统计量 Z |
| `MK_P.tif` | 双侧 p 值 |
| `TrendClass.tif` | 五级分类 1—5，uint8，0=无数据 |
| `trend_area_stats.csv` | 五类面积与占比 |
| `annual_mean.csv` | 完整序列像元上的逐年均值 |
| `fig_annual_curve.png` | 年际曲线 + 线性参考 |
| `fig_trend_map.png` | 五级趋势图（可叠边界） |
| `fig_sen_slope.png` | Sen 斜率连续色带图 |

## 7. 与 GEE 结果对照

山西案例用本脚本复算，得到：

- 完整序列像元 767,787 / 1,326,180（57.9%）；
- 有效面积 151,728.7 km²（全省 96.8%）；
- Sen 斜率中位数 134.6 原始单位/年 = **13.46 g C/m²/yr²**；
- p2 = 2.95，p98 = 22.91 g C/m²/yr²；
- 显著像元占比 96.9%；
- 五类面积与 `results/trend_area_stats.csv` 完全一致。

这些数字与项目早期（未重构前）的运行结果逐位一致，说明本地算法的重构没有改变口径。

## 8. 性能与内存

- 每块 `(22, 60000)` 的 float64 约 10.6 MB；
- ties 的 `eq` 是 `(22, 22, 60000)` 布尔数组，约 29 MB；
- 总内存主要由整幅堆栈决定，而不是块大小；
- 1 万 × 1 万像元、22 年 float64 堆栈约 17.6 GB，需改用 `float32` 或逐行/分块读栅格。

当前脚本为清晰起见先把整幅堆栈读进内存（案例 961×1380 很小）。
超大研究区建议改成：按窗口读取 → 计算 → 写回，或使用 `rasterio` 的窗口迭代。

## 9. 自测

```bash
python tests/test_algorithms.py
```

6 项测试：

1. ties 修正 = 教科书公式（含真实 0）；
2. Sen 斜率 = 显式配对中位数；
3. S = 暴力双重循环；
4. |Z|=1.96 → p=0.05；
5. 五级分类边界；
6. 常量序列 → Z=0, p=1, 类 3。
