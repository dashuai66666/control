"""状态空间表示与仿真（State-Space Representation）。

用途
----
现代控制理论用一阶微分（差分）方程组描述系统：

    ẋ = A·x + B·u
    y = C·x + D·u

状态空间表示能刻画多输入多输出系统与系统内部状态，比传递函数更一般，是
极点配置、LQR、卡尔曼滤波等方法的基础。

本实现用前向欧拉法把连续系统离散化后逐步仿真，展示阶跃响应与稳态增益。
"""

from _control import matmul, inverse, mat_scale, ss_simulate


def step_response(A, B, C, D, dt=0.01, T=10.0):
    """单位阶跃输入下的输出序列。"""
    steps = int(T / dt)
    u = [1.0] * steps
    _, ys = ss_simulate(A, B, C, D, u, dt)
    return ys


def steady_state_gain(A, B, C, D):
    """DC 增益 = -C·A⁻¹·B + D（连续系统阶跃稳态值）。"""
    return -matmul(C, matmul(inverse(A), B))[0][0] + D[0][0]


def main():
    # 二阶系统 s²+3s+2 的能控标准型
    A = [[0.0, 1.0], [-2.0, -3.0]]
    B = [[0.0], [1.0]]
    C = [[1.0, 0.0]]
    D = [[0.0]]

    ys = step_response(A, B, C, D, dt=0.01, T=10.0)
    print("二阶系统状态空间仿真:")
    print("  稳态值 = %.4f（理论 DC 增益 = %.4f）" % (ys[-1], steady_state_gain(A, B, C, D)))
    print("  前几个采样点:", ["%.3f" % y for y in ys[:5]])


if __name__ == "__main__":
    main()
