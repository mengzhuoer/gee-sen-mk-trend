# -*- coding: utf-8 -*-
"""
本地版 Sen's slope + Mann–Kendall 趋势分析。

适用场景：GEE 侧一次请求做完整趋势计算会触发
`Earth Engine capacity exceeded`，或需要完全离线复现。
数据仍由 GEE 导出，本脚本只负责计算与出图。

输入：同一网格的年度栅格序列，文件名含年份占位符。
默认示例（山西省 GPP 案例）：
    data/annual/GPP_raw_annual_2001.tif ... GPP_raw_annual_2022.tif

用法：
    python code/python/sen_mk_local.py \
        --data-dir data/annual \
        --pattern "GPP_raw_annual_{year}.tif" \
        --start-year 2001 --end-year 2022 \
        --scale 0.1 \
        --boundary data/shanxi.geojson \
        --out results

`--scale` 是“原始整数 -> 物理值”的比例因子；GPP 案例中
原始整数 × 0.1 = g C/m²/year。NDVI 案例若输入已是物理值，用 1.0。
"""

from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import rasterio
from rasterio.transform import array_bounds

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sen_mk_core import CLASS_NAMES, analyze  # noqa: E402

RE = 6371007.181  # WGS84 授权球半径（m），球面面积近似误差 < 0.5%


def parse_args():
    p = argparse.ArgumentParser(description="Sen + MK 本地趋势分析")
    p.add_argument("--data-dir", default="data/annual")
    p.add_argument("--pattern", default="GPP_raw_annual_{year}.tif",
                   help="文件名模板，必须含 {year}")
    p.add_argument("--start-year", type=int, default=2001)
    p.add_argument("--end-year", type=int, default=2022)
    p.add_argument("--scale", type=float, default=0.1,
                   help="原始值 -> 物理值比例因子")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--chunk", type=int, default=60_000,
                   help="每块处理的像元数")
    p.add_argument("--out", default="results")
    p.add_argument("--boundary", default=None,
                   help="可选：GeoJSON 边界，用于出图叠加")
    p.add_argument("--prefix", default="",
                   help="输出文件名前缀（多变量时避免覆盖）")
    p.add_argument("--no-plots", action="store_true")
    return p.parse_args()


def load_stack(args):
    years = list(range(args.start_year, args.end_year + 1))
    paths = []
    for y in years:
        p = os.path.join(args.data_dir, args.pattern.format(year=y))
        if not os.path.exists(p):
            raise SystemExit(f"缺少年份文件: {p}")
        paths.append(p)

    ref = None
    stack = []
    profile = None
    for p in paths:
        with rasterio.open(p) as ds:
            key = (ds.width, ds.height, str(ds.crs), tuple(ds.transform),
                   ds.count, ds.dtypes[0], ds.nodata)
            if ref is None:
                ref, profile = key, ds.profile.copy()
            elif key != ref:
                raise SystemExit(f"网格不一致，不能叠栈: {os.path.basename(p)}")
            arr = ds.read(1)
            if ds.nodata is not None:
                arr = np.where(arr == ds.nodata, np.nan, arr)
            stack.append(arr.astype(np.float64))

    print(f"[1/5] 读取 {len(stack)} 年：{profile['width']}×{profile['height']}，"
          f"网格一致 ✓")
    return np.stack(stack), profile, years


def pixel_area(profile):
    """返回每行像元面积（m²）。经纬度网格用球面公式，投影网格用仿射行列式。"""
    T = profile["transform"]
    h, w = profile["height"], profile["width"]
    crs = profile["crs"]
    is_geographic = bool(crs and crs.is_geographic)
    if not is_geographic:
        area = abs(T.a * T.e - T.b * T.d)
        print(f"      注意：CRS={crs} 非经纬度，面积按仿射行列式 {area:.1f} m² 估算")
        return np.full(h, area)

    dlon = np.deg2rad(abs(T.a))
    lat1 = T.f + np.arange(h + 1) * T.e
    lat_hi = np.maximum(lat1[:-1], lat1[1:])
    lat_lo = np.minimum(lat1[:-1], lat1[1:])
    return RE ** 2 * dlon * (np.sin(np.deg2rad(lat_hi)) -
                             np.sin(np.deg2rad(lat_lo)))


def compute(stack, args, years):
    n_years, h, w = stack.shape
    total = h * w
    flat = stack.reshape(n_years, total)
    complete = ~np.isnan(flat).any(axis=0)
    cols = np.flatnonzero(complete)
    print(f"[2/5] 完整序列像元：{cols.size:,} / {total:,} "
          f"({100 * cols.size / total:.1f}%)")

    slope = np.full(total, np.nan)
    z_out = np.full(total, np.nan)
    p_out = np.full(total, np.nan)
    cls = np.zeros(total, dtype=np.uint8)

    t0 = time.time()
    for s in range(0, cols.size, args.chunk):
        idx = cols[s:s + args.chunk]
        res = analyze(flat[:, idx], years=years, alpha=args.alpha)
        slope[idx] = res["slope"]
        z_out[idx] = res["Z"]
        p_out[idx] = res["p"]
        cls[idx] = res["class"]
        if (s // args.chunk) % 5 == 0:
            print(f"      块 {s // args.chunk + 1}，"
                  f"{s + idx.size:,}/{cols.size:,}")

    print(f"[3/5] 计算完成，用时 {time.time() - t0:.1f} 秒")
    return (slope.reshape(h, w), z_out.reshape(h, w), p_out.reshape(h, w),
            cls.reshape(h, w), complete.reshape(h, w), flat, years)


def write_rasters(profile, out_dir, prefix, slope, z, p, cls):
    def put(name, arr, nodata, dtype):
        prof = profile.copy()
        prof.update(driver="GTiff", count=1, dtype=dtype, nodata=nodata,
                    compress="deflate", tiled=False)
        path = os.path.join(out_dir, prefix + name)
        with rasterio.open(path, "w", **prof) as dst:
            dst.write(arr.astype(dtype), 1)
        print(f"      {path}")

    put("SenSlope.tif", slope, np.nan, "float32")
    put("MK_Z.tif", z, np.nan, "float32")
    put("MK_P.tif", p, np.nan, "float32")
    put("TrendClass.tif", cls, 0, "uint8")
    print("[4/5] 结果栅格已写出")


def area_stats(cls, area_row, out_dir, prefix, scale):
    h, w = cls.shape
    rows = []
    valid_km2 = 0.0
    per_class = {}
    for c in range(1, 6):
        cnt = (cls == c).sum(axis=1).astype(np.float64)
        km2 = float((cnt * area_row).sum() / 1e6)
        per_class[c] = km2
        valid_km2 += km2

    for c in range(1, 6):
        pct = 100 * per_class[c] / valid_km2 if valid_km2 else 0.0
        rows.append((c, CLASS_NAMES[c], per_class[c], pct))

    path = os.path.join(out_dir, prefix + "trend_area_stats.csv")
    with open(path, "w", encoding="utf-8-sig") as f:
        f.write("class_id,trend,area_km2,percent_of_valid\n")
        for c, name, km2, pct in rows:
            f.write(f"{c},{name},{km2:.2f},{pct:.4f}\n")

    print("\n===== 五级趋势面积 =====")
    for c, name, km2, pct in rows:
        print(f"  {c} {name:<6s} {km2:12,.1f} km²  {pct:6.2f}%")
    print(f"  有效分析面积合计 {valid_km2:,.1f} km² -> {path}")
    return rows


def annual_series(flat, complete, years, args, out_dir, prefix):
    x = flat[:, complete.ravel()]
    means = np.nanmean(x, axis=1) * args.scale
    path = os.path.join(out_dir, prefix + "annual_mean.csv")
    with open(path, "w", encoding="utf-8-sig") as f:
        f.write("year,mean\n")
        for y, v in zip(years, means):
            f.write(f"{y},{v:.4f}\n")
    print("\n===== 研究区年均值 =====")
    for y, v in zip(years, means):
        print(f"  {y}: {v:10.2f}")
    return means, path


def make_figures(slope, cls, means, years, args, out_dir, prefix):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import ListedColormap, BoundaryNorm
    except ImportError:
        print("跳过出图：未安装 matplotlib")
        return

    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(years, means, "o-", color="#006837", lw=2, ms=5)
    k, b = np.polyfit(years, means, 1)
    ax.plot(years, np.poly1d([k, b])(years), "--", color="#A50026",
            label=f"线性参考 {k:.2f}/yr")
    ax.set_xlabel("Year")
    ax.set_ylabel("Annual mean")
    ax.set_title(f"{years[0]}–{years[-1]} annual mean")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, prefix + "fig_annual_curve.png"), dpi=200)
    plt.close(fig)

    boundary = None
    if args.boundary and os.path.exists(args.boundary):
        try:
            import geopandas as gpd
            boundary = gpd.read_file(args.boundary)
        except Exception as e:  # noqa: BLE001
            print(f"边界读取失败，跳过叠加：{e}")

    cmap = ListedColormap(["#A50026", "#F46D43", "#D9D9D9",
                           "#66BD63", "#006837"])
    norm = BoundaryNorm([0.5, 1.5, 2.5, 3.5, 4.5, 5.5], cmap.N)
    fig, ax = plt.subplots(figsize=(7, 8))
    im = ax.imshow(np.where(cls == 0, np.nan, cls), cmap=cmap, norm=norm)
    if boundary is not None:
        boundary.boundary.plot(ax=ax, color="k", lw=0.8)
    ax.set_title(f"Sen + MK trend class ({years[0]}–{years[-1]}, "
                 f"alpha={args.alpha})")
    ax.set_axis_off()
    cb = fig.colorbar(im, ax=ax, shrink=0.6, ticks=[1, 2, 3, 4, 5])
    cb.ax.set_yticklabels([f"{k} {v}" for k, v in CLASS_NAMES.items()],
                          fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, prefix + "fig_trend_map.png"), dpi=200)
    plt.close(fig)

    vals = slope[np.isfinite(slope)]
    lo, hi = np.percentile(vals, [2, 98])
    fig, ax = plt.subplots(figsize=(7, 8))
    im = ax.imshow(slope, cmap="RdBu_r", vmin=lo, vmax=hi)
    if boundary is not None:
        boundary.boundary.plot(ax=ax, color="k", lw=0.8)
    ax.set_title("Sen slope (2%–98% stretch)")
    ax.set_axis_off()
    fig.colorbar(im, ax=ax, shrink=0.6)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, prefix + "fig_sen_slope.png"), dpi=200)
    plt.close(fig)
    print("      三张图已输出")


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    print("=" * 62)
    print("Sen's slope + Mann–Kendall 本地趋势分析")
    print("=" * 62)

    stack, profile, years = load_stack(args)
    slope_raw, z, p, cls, complete, flat, years = compute(stack, args, years)

    v = slope_raw[np.isfinite(slope_raw)]
    q = np.percentile(v, [2, 25, 50, 75, 98])
    print("\nSen 斜率分位数（原始单位/年）：")
    print("  " + "  ".join(f"p{a}={b:.4f}" for a, b in zip([2, 25, 50, 75, 98], q)))
    print("  物理单位/年：")
    print("  " + "  ".join(f"p{a}={b * args.scale:.4f}"
                           for a, b in zip([2, 25, 50, 75, 98], q)))
    sig = np.isfinite(p) & (p < args.alpha)
    print(f"显著像元占比（p<{args.alpha}）："
          f"{100 * sig.sum() / max(complete.sum(), 1):.1f}%")

    write_rasters(profile, args.out, args.prefix,
                  slope_raw * args.scale, z, p, cls)
    area_row = pixel_area(profile)
    area_stats(cls, area_row, args.out, args.prefix, args.scale)
    means, _ = annual_series(flat, complete, years, args, args.out, args.prefix)

    if not args.no_plots:
        make_figures(slope_raw * args.scale, cls, means, years, args,
                     args.out, args.prefix)

    print("\n[5/5] 全部完成。输出目录：", os.path.abspath(args.out))


if __name__ == "__main__":
    main()
