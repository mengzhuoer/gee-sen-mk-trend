# -*- coding: utf-8 -*-
"""
Sen's slope + Mann–Kendall 核心算法（纯 numpy，无 I/O）。

输入约定：x 形状为 (n_years, n_pixels) 或 (n_years,)。
缺失值请在此之前剔除（本仓库采用“完整序列”策略，不插值、不补 0）。

方法要点
--------
Sen 斜率
    β = median{ (x_j - x_i) / (j - i) : 1 <= i < j <= n }
    对异常值稳健，不要求正态分布。

Mann–Kendall
    S = Σ_{i<j} sign(x_j - x_i)
    Var(S) = [n(n-1)(2n+5) - Σ_p t_p(t_p-1)(2t_p+5)] / 18
    其中 t_p 是第 p 组重复值的个数（ties 修正）。
    连续性修正：S>0 时 Z=(S-1)/√Var；S<0 时 Z=(S+1)/√Var；S=0 时 Z=0。
    双侧 p = 2 * (1 - Φ(|Z|))。

分类
    1 显著下降 / 2 不显著下降 / 3 基本稳定 / 4 不显著上升 / 5 显著上升。
"""

from __future__ import annotations

import numpy as np
from scipy.stats import norm

CLASS_NAMES = {
    1: "显著下降",
    2: "不显著下降",
    3: "基本稳定",
    4: "不显著上升",
    5: "显著上升",
}


def _as_2d(x) -> np.ndarray:
    """把输入统一成 (n_years, n_pixels)；1 维输入视为单像元。"""
    arr = np.asarray(x, dtype=np.float64)
    if arr.ndim == 1:
        return arr[:, None]
    if arr.ndim != 2:
        raise ValueError(f"x 必须是 1 维或 2 维，收到 {arr.ndim} 维")
    return arr


def sen_slope(x, years=None) -> np.ndarray:
    """计算 Sen's slope。

    Parameters
    ----------
    x : array_like, shape (n,) or (n, k)
        时间序列，按年份升序排列。
    years : array_like, optional
        对应年份；缺省用 0..n-1，此时斜率单位为“每年”。

    Returns
    -------
    ndarray, shape (k,)
        每个像元的斜率，单位 = x 的单位 / 年。
    """
    arr = _as_2d(x)
    n = arr.shape[0]
    if years is None:
        t = np.arange(n, dtype=np.float64)
    else:
        t = np.asarray(years, dtype=np.float64)
        if t.shape != (n,):
            raise ValueError("years 长度必须等于时间序列长度")

    i, j = np.triu_indices(n, 1)
    dt = (t[j] - t[i])[:, None]          # (n_pairs, 1)
    slopes = (arr[j] - arr[i]) / dt      # (n_pairs, k)
    return np.median(slopes, axis=0)


def mk_test(x, continuity: bool = True):
    """Mann–Kendall 趋势检验（含 ties 修正）。

    Parameters
    ----------
    x : array_like, shape (n,) or (n, k)
    continuity : bool
        是否使用连续性修正（默认 True，与主流实现一致）。

    Returns
    -------
    S : ndarray, shape (k,)
    var_s : ndarray, shape (k,)
        修正后的方差。
    z : ndarray, shape (k,)
    p : ndarray, shape (k,)
        双侧 p 值。Var(S)=0（常量序列）时定义 Z=0, p=1。
    """
    arr = _as_2d(x)
    n = arr.shape[0]
    i, j = np.triu_indices(n, 1)

    d = arr[j] - arr[i]
    s = np.sign(d).sum(axis=0)

    # 等价于教科书式：Σ_p t_p(t_p-1)(2t_p+5)
    # 证明：大小为 t 的一组在“逐年重复计数”中出现 t 次，
    # 每次贡献 (t-1)(2t+5)，合计 t(t-1)(2t+5)。
    eq = arr[:, None, :] == arr[None, :, :]   # (n, n, k)
    m = eq.sum(axis=1)                        # 每一年取值在该像元的重复数
    tie_sum = ((m - 1) * (2 * m + 5)).sum(axis=0)

    base_var = n * (n - 1) * (2 * n + 5) / 18.0
    var_s = np.maximum(base_var - tie_sum / 18.0, 0.0)

    if continuity:
        numerator = np.where(s > 0, s - 1.0,
                             np.where(s < 0, s + 1.0, 0.0))
    else:
        numerator = s.astype(np.float64)

    safe = np.where(var_s > 0, np.sqrt(var_s), 1.0)
    z = np.where(var_s > 0, numerator / safe, 0.0)
    p = np.where(var_s > 0, 2.0 * norm.sf(np.abs(z)), 1.0)
    return s, var_s, z, p


def classify(slope, p, alpha: float = 0.05) -> np.ndarray:
    """五级趋势分类。

    返回 uint8 数组：0=无数据/无效，1..5 见 CLASS_NAMES。
    """
    slope = np.asarray(slope, dtype=np.float64)
    p = np.asarray(p, dtype=np.float64)
    slope, p = np.broadcast_arrays(slope, p)

    out = np.zeros(slope.shape, dtype=np.uint8)
    valid = np.isfinite(slope) & np.isfinite(p)
    sig = p < alpha

    out[valid & (slope < 0) & sig] = 1
    out[valid & (slope < 0) & ~sig] = 2
    out[valid & (slope == 0)] = 3
    out[valid & (slope > 0) & ~sig] = 4
    out[valid & (slope > 0) & sig] = 5
    return out


def analyze(x, years=None, alpha: float = 0.05):
    """一次算完 Sen + MK + 分类，返回字典。"""
    slope = sen_slope(x, years=years)
    s, var_s, z, p = mk_test(x)
    cls = classify(slope, p, alpha=alpha)
    return {"slope": slope, "S": s, "var_s": var_s, "Z": z, "p": p,
            "class": cls}
