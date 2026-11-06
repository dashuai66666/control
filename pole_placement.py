"""极点配置（Pole Placement / State Feedback）。

用途
----
给定期望的闭环极点位置（决定系统的响应速度与阻尼），设计全状态反馈律
u = -K·x，使闭环系统 A - B·K 的特征值恰好落在期望位置。这是现代控制里
最直接的设计方法，前提是系统能控。

数学原理（Ackermann 公式，单输入）
----------------------------------
设期望特征多项式为 φ(s) = sⁿ + aₙ₋₁sⁿ⁻¹ + ... + a₀，则反馈增益

    K = [0 ... 0 1] · Wc⁻¹ · φ(A)

其中 Wc 为能控性矩阵，φ(A) = Aⁿ + aₙ₋₁Aⁿ⁻¹ + ... + a₀I（把 A 代入特征多项式）。

本实现据此计算 K，并仿真闭环验证极点确实被移动。
"""

from _control import (
    identity, zeros, matmul, mat_add, mat_sub, mat_scale, inverse,
    matrix_power, controllability_matrix, poly_mul, poly_roots, ss_simulate,
)


def ackermann(A, B, poles):
    """Ackermann 公式求状态反馈增益 K（1×n）。"""
    n = len(A)
    # 期望特征多项式（首一，降幂）
    char = [1.0]
    for p in poles:
        char = poly_mul(char, [1.0, -p])

    # φ(A) = Aⁿ + a_{n-1}A^{n-1} + ... + a0·I
    phi = matrix_power(A, n)
    for k, c in enumerate(char[1:]):            # a_{n-1}, ..., a_0
        power = n - 1 - k
        term = identity(n) if power == 0 else matrix_power(A, power)
        phi = mat_add(phi, mat_scale(term, c))

    Wc_inv = inverse(controllability_matrix(A, B))
    last = Wc_inv[-1]                            # e_nᵀ · Wc⁻¹ 即 Wc⁻¹ 的最后一行
    K = [sum(last[j] * phi[j][i] for j in range(n)) for i in range(n)]
    return [K]


def main():
    # 原系统：A 的特征值为 0 与 1（不稳定）
    A = [[0.0, 1.0], [0.0, 1.0]]
    B = [[0.0], [1.0]]
    C = [[1.0, 0.0]]
    D = [[0.0]]

    # 期望闭环极点：-2 与 -3（稳定且更快）
    K = ackermann(A, B, [-2.0, -3.0])
    Acl = mat_sub(A, matmul(B, K))

    # 2×2 矩阵的特征多项式：λ² - tr·λ + det
    tr = Acl[0][0] + Acl[1][1]
    det = Acl[0][0] * Acl[1][1] - Acl[0][1] * Acl[1][0]
    actual = poly_roots([1.0, -tr, det])

    print("状态反馈增益 K =", ["%.3f" % v for v in K[0]])
    print("配置后闭环极点:", ["%.3f" % z.real for z in actual])
    print("（期望极点 -2, -3）")


if __name__ == "__main__":
    main()
