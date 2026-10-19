# -*- coding: utf-8 -*-
"""
算法自测（无需 GEE、无需栅格数据）：

    python tests/test_algorithms.py

覆盖：
1. ties 修正与教科书公式逐例一致；
2. Sen 斜率与暴力双重循环一致；
3. MK 的 S 与暴力符号和一致；
4. |Z| = 1.96 时双侧 p ≈ 0.05；
5. 五级分类边界；
6. 常量序列（VarS=0）→ Z=0, p=1, 类 3。
"""

import os
import sys

import numpy as np
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "code", "python"))

from sen_mk_core import analyze, classify, mk_test, sen_slope  # noqa: E402


def textbook_tie_sum(x):
    """教科书式 Σ t(t-1)(2t+5)。"""
    _, counts = np.unique(x, return_counts=True)
    return sum(t * (t - 1) * (2 * t + 5) for t in counts)


def test_tie_correction():
    rng = np.random.default_rng(42)
    for n in (5, 10, 22, 24):
        for _ in range(300):
            # 小整数 → 必然产生 ties，且包含真实 0
            x = rng.integers(0, 4, size=n).astype(float)
            _, var_s, _, _ = mk_test(x)
            base = n * (n - 1) * (2 * n + 5) / 18.0
            expected = base - textbook_tie_sum(x) / 18.0
            assert abs(var_s[0] - max(expected, 0.0)) < 1e-9, (x, var_s, expected)
    print("PASS  tie correction == textbook formula")


def test_sen_slope():
    rng = np.random.default_rng(7)
    for n in (5, 12, 22):
        x = rng.normal(size=n)
        got = sen_slope(x)[0]
        slopes = [(x[j] - x[i]) / (j - i)
                  for i in range(n) for j in range(i + 1, n)]
        assert abs(got - np.median(slopes)) < 1e-12
    # 2 维批量
    x2 = rng.normal(size=(22, 9))
    got2 = sen_slope(x2)
    for k in range(9):
        assert abs(got2[k] - sen_slope(x2[:, k])[0]) < 1e-12
    print("PASS  Sen slope == explicit pairwise median")


def test_mk_s():
    rng = np.random.default_rng(11)
    for n in (5, 10, 22):
        x = rng.integers(0, 6, size=n).astype(float)
        s, _, _, _ = mk_test(x)
        brute = sum(np.sign(x[j] - x[i]) for i in range(n) for j in range(i + 1, n))
        assert int(s[0]) == int(brute)
    print("PASS  MK S == brute-force sign sum")


def test_p_value():
    # 反解：p=0.05 对应 |Z|=1.959964...
    z = norm.ppf(1 - 0.05 / 2)
    p = 2 * norm.sf(z)
    assert abs(p - 0.05) < 1e-12
    print(f"PASS  |Z|={z:.6f} -> p={p:.6f}")


def test_classification():
    slope = np.array([-2.0, -2.0, 0.0, 2.0, 2.0, np.nan])
    p = np.array([0.01, 0.50, 0.50, 0.50, 0.01, 0.01])
    got = classify(slope, p)
    assert got.tolist() == [1, 2, 3, 4, 5, 0], got
    print("PASS  five-class assignment")


def test_constant_series():
    x = np.full(22, 7.0)
    res = analyze(x)
    assert res["S"][0] == 0
    assert res["var_s"][0] == 0
    assert res["Z"][0] == 0
    assert res["p"][0] == 1
    assert res["class"][0] == 3
    print("PASS  constant series -> Z=0, p=1, class=3")


if __name__ == "__main__":
    test_tie_correction()
    test_sen_slope()
    test_mk_s()
    test_p_value()
    test_classification()
    test_constant_series()
    print("\nALL TESTS PASSED")
