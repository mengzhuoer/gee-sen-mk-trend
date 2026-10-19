# 山西 GPP 案例的分阶段 GEE 脚本

这 5 个脚本是山西 2001—2022 GPP 案例在 GEE 侧实际跑通的版本。
它们不是通用库，而是**可核查的取数与核验记录**；趋势计算在本地完成
（`code/python/sen_mk_local.py`）。

## 运行顺序

| 脚本 | 作用 | 关键输出 |
|---|---|---|
| `01_check_annual_gpp.js` | 边界、影像数、每年 46 时段、波段与比例因子检查 | Console 诊断值；仅显示 2022 |
| `02_export_annual_gpp_to_assets.js` | 单年导出到 GEE Assets，验证导出链路 | `Shanxi_GAUL_GPP_raw_annual_2022` |
| `03_download_remaining_years.js` | 其余年份逐年生成下载链接 | 21 条签名 URL |
| `04_fingerprint_years.js` | 在与下载文件一致的网格上打印年度 8 点指纹 | `YYYY => [8 个值]` |
| `05_fill_missing_years.js` | 补下缺失年份（案例为 2017、2022） | 指纹 + 单年链接 |

## 前置条件

- GEE 工作区已绑定 Cloud Project（否则资产面板会报认证错误）；
- 使用 GAUL 边界：`FAO/GAUL/2015/level1`，太原点 `[112.549, 37.871]`
  对应 `Shanxi Sheng`；
- 数据：`MODIS/061/MOD17A2HGF`，8 天，`Gpp` 波段；
- 年度指标：日历年 46 个时段全部有效 → 全年累加，Int32 原始整数；
- 物理值：原始整数 × 0.1 = g C/m²/year。

## 注意

- `02` 的 `ASSET_ROOT = 'projects/songyuan/assets'` 需要换成你自己的项目路径；
- `03` / `04` / `05` 里的网格常量（EPSG:4326 + transform）是本案下载文件的网格，
  换研究区必须重新确定，不能照抄；
- `getDownloadURL` 的 `filename` 参数会被忽略，下载文件是哈希名，
  必须用指纹认领年份；
- 大区域一次 Run 容易触发 `Earth Engine capacity exceeded`，保持“一脚本一件事”。
