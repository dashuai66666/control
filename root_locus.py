"""根轨迹（Root Locus）。

用途
----
当反馈增益 K 从 0 变化到 ∞ 时，闭环极点（特征方程 1 + K·G(s) = 0 的根）在
复平面上的运动轨迹。根轨迹直观地展示了增益对系统稳定性与动态性能的影响：
极点一旦越过虚轴进入右半平面，系统就不稳定。

数学原理
--------
闭环特征方程可写为：

    den(s) + K·num(s) = 0

对每个 K 求该多项式的根即得一组闭环极点。把不同 K 下的根连起来就是根轨迹。
本实现用 Durand-Kerner 方法（见 ``_control.poly_roots``）对若干 K 值求根。

典型规律
--------
- K→0 时根趋近开环极点；
- 根轨迹分支始于开环极点、终于开环零点（或无穷远）；
- 分支数等于 max(极点个数, 零点个数)。
"""

from _control import poly_roots, poly_add, poly_mul


def root_locus(num, den, Ks):
    """对每个增益 K 求闭环极点，返回 [(K, roots)]。"""
    result = []
    for K in Ks:
        # 特征多项式 = den + K·num
        char = poly_add(den, [K * c for c in num])
        result.append((K, poly_roots(char)))
    return result


def main():
    # G(s) = 1/(s(s+1)(s+2))：三条根轨迹，其中一条随 K 增大穿过虚轴（不稳定）
    num, den = [1.0], [1.0, 3.0, 2.0, 0.0]
    Ks = [0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0]

    print("G(s) = 1/(s(s+1)(s+2)) 的根轨迹（随 K 增大）:")
    for K, roots in root_locus(num, den, Ks):
        roots = sorted(roots, key=lambda z: (z.real, z.imag))
        shown = "  ".join("(%+.2f%+.2fj)" % (z.real, z.imag) for z in roots)
        print("  K=%-4.1f : %s" % (K, shown))


if __name__ == "__main__":
    main()
