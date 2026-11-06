"""Routh-Hurwitz 稳定性判据（Routh-Hurwitz Criterion）。

用途
----
不需求解特征方程的根，仅通过系数构造 Routh 表，即可判断系统闭环极点是否全部
位于复平面左半平面（LHP），从而判定稳定性。这是经典控制中最常用的代数判据。

数学原理
--------
对特征多项式 p(s) = a0·sⁿ + a1·sⁿ⁻¹ + ... + an，构造 Routh 表：前两行分别
放偶数与奇数下标系数，之后每行由前两行按行列式规则生成：

    b1 = (a1·a2 - a0·a3) / a1,  b2 = (a1·a4 - a0·a5) / a1, ...

**稳定性判据**：Routh 表第一列元素符号改变的次数 = 位于右半平面的根的个数。
因此系统稳定 ⇔ 第一列所有元素同号（且均非零）。

特殊情况（第一列出现 0）可用 ε 小量替代后再判断，本实现用极小正数近似处理。
"""


def routh_first_column(coeffs):
    """返回 Routh 表第一列元素。coeffs 按降幂排列（首项为 sⁿ 系数）。"""
    n = len(coeffs) - 1
    rows = [list(coeffs[0::2]), list(coeffs[1::2])]     # 前两行
    L = len(rows[0])
    for r in rows:
        while len(r) < L:
            r.append(0.0)

    for _ in range(2, n + 1):
        prev2, prev = rows[-2], rows[-1]
        if abs(prev[0]) < 1e-12:
            prev[0] = 1e-9                                # ε 法处理零元素
        new = []
        for j in range(L - 1):
            num = prev[0] * prev2[j + 1] - prev2[0] * prev[j + 1]
            new.append(num / prev[0])
        new.append(0.0)
        rows.append(new)

    return [rows[i][0] for i in range(len(rows))]


def count_unstable_poles(coeffs):
    """返回右半平面根的个数（第一列符号变化次数）。"""
    col = routh_first_column(coeffs)
    changes = 0
    for i in range(1, len(col)):
        if col[i - 1] * col[i] < 0:
            changes += 1
    return changes


def main():
    # 三个三阶多项式：稳定、临界稳定(K=6)、不稳定(K=7)
    examples = [
        ("s³+4s²+5s+2 (稳定)", [1.0, 4.0, 5.0, 2.0]),
        ("s³+3s²+2s+6 (临界)", [1.0, 3.0, 2.0, 6.0]),
        ("s³+3s²+2s+7 (不稳定)", [1.0, 3.0, 2.0, 7.0]),
    ]
    for name, coeffs in examples:
        k = count_unstable_poles(coeffs)
        print("%s -> 右半平面根数 = %d, %s" % (name, k, "稳定" if k == 0 else "不稳定"))


if __name__ == "__main__":
    main()
