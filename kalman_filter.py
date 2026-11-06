"""卡尔曼滤波（Kalman Filter）。

用途
----
从含噪声的观测中**最优估计**线性系统的内部状态，是状态估计与数据融合的基石，
广泛用于导航、目标跟踪、传感器融合等。

数学模型
--------
系统方程：x[k+1] = A·x[k] + B·u[k] + w   （w 为过程噪声，协方差 Q）
观测方程：z[k]   = C·x[k] + v            （v 为观测噪声，协方差 R）

卡尔曼滤波交替执行两步：

1. **预测**：x̂ = A·x̂ + B·u，P = A·P·Aᵀ + Q
2. **更新**：K = P·Cᵀ·(C·P·Cᵀ + R)⁻¹，x̂ ← x̂ + K·(z - C·x̂)，P ← (I - K·C)·P

K 为卡尔曼增益，它按"过程噪声/观测噪声"的相对大小自适应地权衡预测与观测。
"""

import random
import math

from _control import matmul, mat_add, mat_sub, transpose, inverse, mat_scale, identity


def kalman_filter(A, B, C, Q, R, zs, us=None, x0=None, P0=None):
    """离散卡尔曼滤波，返回 (估计状态轨迹, 估计协方差轨迹)。"""
    n = len(A)
    x = [0.0] * n if x0 is None else list(x0)
    P = identity(n) if P0 is None else [row[:] for row in P0]
    xs = []
    for k, z in enumerate(zs):
        u = us[k] if us else 0.0
        # 预测
        x = [sum(A[i][j] * x[j] for j in range(n)) + B[i][0] * u for i in range(n)]
        P = mat_add(matmul(A, matmul(P, transpose(A))), Q)
        # 更新
        Cx = [sum(C[0][j] * x[j] for j in range(n))]
        S = mat_add(matmul(C, matmul(P, transpose(C))), R)
        K = matmul(matmul(P, transpose(C)), inverse(S))      # n×1
        innov = z - Cx[0]
        x = [x[i] + K[i][0] * innov for i in range(n)]
        P = mat_sub(P, matmul(K, matmul(C, P)))
        xs.append(list(x))
    return xs


def main():
    rng = random.Random(0)
    # 匀速运动模型：状态 = [位置, 速度]
    dt = 1.0
    A = [[1.0, dt], [0.0, 1.0]]
    B = [[0.0], [0.0]]
    C = [[1.0, 0.0]]                       # 只观测位置
    Q = [[0.1, 0.0], [0.0, 0.1]]           # 过程噪声
    R = [[1.0]]                            # 观测噪声（较大）

    # 生成真实轨迹与含噪观测
    true = [0.0, 1.0]                      # 初始位置 0，速度 1
    zs = []
    for _ in range(40):
        true = [sum(A[i][j] * true[j] for j in range(2)) + rng.gauss(0, 0.3) for i in range(2)]
        zs.append(true[0] + rng.gauss(0, 1.0))

    xs = kalman_filter(A, B, C, Q, R, zs)
    # 比较最后一步：观测噪声误差 vs 卡尔曼估计误差（与真实值）
    err_obs = abs(zs[-1] - true[0])
    err_kf = abs(xs[-1][0] - true[0])
    print("卡尔曼滤波（位置估计）:")
    print("  观测误差 %.3f -> 估计误差 %.3f" % (err_obs, err_kf))


if __name__ == "__main__":
    main()
