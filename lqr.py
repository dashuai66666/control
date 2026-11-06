"""线性二次型最优控制（Linear Quadratic Regulator, LQR）。

用途
----
在状态空间框架下，求使二次型代价函数最小的状态反馈律 u = -K·x：

    J = Σ (xᵀQx + uᵀRu)

其中 Q 加权状态偏差（通常半正定），R 加权控制能量（正定）。LQR 兼顾"控制性能"
与"控制代价"，是最经典的最优控制方法，广泛用于航天、机器人等。

数学原理（离散时间）
--------------------
对离散系统 x[k+1] = A x[k] + B u[k]，最优增益可由离散代数 Riccati 方程求得。
本实现用**值迭代**直接迭代 Riccati 差分方程直至收敛：

    P ← AᵀPA - AᵀPB (R + BᵀPB)⁻¹ BᵀPA + Q
    K = (R + BᵀPB)⁻¹ BᵀPA

收敛后 K 即最优反馈增益。
"""

from _control import matmul, mat_add, mat_sub, transpose, inverse, mat_scale, identity


def lqr_discrete(A, B, Q, R, iters=2000, tol=1e-10):
    """离散 LQR：返回最优增益 K（1×n）。"""
    P = [row[:] for row in Q]                 # 初始化 P = Q
    At = transpose(A)
    Bt = transpose(B)
    for _ in range(iters):
        BPB = matmul(Bt, matmul(P, B))        # R + BᵀPB
        S = [[R[i][j] + BPB[i][j] for j in range(len(R))] for i in range(len(R))]
        S_inv = inverse(S)
        AtPB = matmul(At, matmul(P, B))
        K = matmul(S_inv, transpose(AtPB))    # K = (R+BᵀPB)⁻¹ BᵀPA
        P_new = mat_add(matmul(At, matmul(P, A)), Q)
        P_new = mat_sub(P_new, matmul(AtPB, K))
        # 收敛判断
        diff = max(abs(P_new[i][j] - P[i][j])
                   for i in range(len(P)) for j in range(len(P)))
        P = P_new
        if diff < tol:
            break
    K = matmul(inverse(mat_add(R, matmul(Bt, matmul(P, B)))),
               matmul(Bt, matmul(P, A)))
    return K


def main():
    # 离散双积分器：x[k+1] = A x + B u
    A = [[1.0, 1.0], [0.0, 1.0]]
    B = [[0.0], [1.0]]
    Q = identity(2)          # 惩罚状态偏差
    R = [[0.1]]              # 惩罚控制能量

    K = lqr_discrete(A, B, Q, R)
    print("离散 LQR 最优增益 K =", ["%.3f" % v for v in K[0]])

    # 仿真闭环 x[k+1] = (A - B K) x，验证状态收敛
    Acl = mat_sub(A, matmul(B, K))
    x = [1.0, 0.5]           # 初始状态
    for k in range(20):
        x = [sum(Acl[i][j] * x[j] for j in range(2)) for i in range(2)]
    print("20 步后状态 x =", ["%.4f" % v for v in x], "（应趋于 0）")


if __name__ == "__main__":
    main()
