"""Ziegler-Nichols 临界比例度整定法（Ziegler-Nichols Ultimate Gain Method）。

用途
----
PID 参数整定的经典工程方法：不用建立精确模型，只通过实验找到**临界增益**与
**临界振荡周期**，再按经验表直接得到 P / PI / PID 参数。

步骤
----
1. 只加比例控制（Ki=Kd=0），逐步增大 Kp；
2. 直到输出出现**等幅振荡**（临界稳定），记此时的增益为临界增益 Ku，
   振荡周期为临界周期 Tu；
3. 按 Z-N 表计算参数：

    P   : Kp = 0.5·Ku
    PI  : Kp = 0.45·Ku,  Ti = 0.83·Tu
    PID : Kp = 0.6·Ku,   Ti = 0.5·Tu,  Td = 0.125·Tu   （Ki=Kp/Ti, Kd=Kp·Td）

本实现用三阶对象 1/(s(s+1)(s+2)) 做闭环仿真，通过二分搜索逼近 Ku，再从
临界振荡中测出 Tu，演示整个整定流程。
"""

from _control import tf_to_ss


def _p_loop_response(plant, K, setpoint, dt, steps):
    """比例控制闭环仿真，返回输出序列。"""
    A, B, C, D = plant
    n = len(A)
    x = [0.0] * n
    ys = []
    for _ in range(steps):
        y = sum(C[0][j] * x[j] for j in range(n))
        u = K * (setpoint - y)
        Ax = [sum(A[i][j] * x[j] for j in range(n)) for i in range(n)]
        Bu = [B[i][0] * u for i in range(n)]
        x = [x[i] + dt * (Ax[i] + Bu[i]) for i in range(n)]
        ys.append(y)
    return ys


def _peaks(ys):
    """提取局部极大值（采样点下标, 峰值）。"""
    out = []
    for i in range(1, len(ys) - 1):
        if ys[i] > ys[i - 1] and ys[i] >= ys[i + 1]:
            out.append((i, ys[i]))
    return out


def find_ultimate(plant, dt=0.01, T=80.0, setpoint=1.0):
    """二分搜索临界增益 Ku，并测出临界振荡周期 Tu。"""
    steps = int(T / dt)

    def growing(K):
        ys = _p_loop_response(plant, K, setpoint, dt, steps)
        pk = _peaks(ys)
        if len(pk) < 6:
            return False
        first = sum(p[1] for p in pk[:3]) / 3
        last = sum(p[1] for p in pk[-3:]) / 3
        return last > first * 1.02            # 幅值持续增大 → 不稳定

    lo, hi = 0.0, 1.0
    while not growing(hi) and hi < 1e6:
        hi *= 2.0                             # 先扩大 hi 直到不稳定
    for _ in range(35):
        mid = (lo + hi) / 2
        if growing(mid):
            hi = mid
        else:
            lo = mid
    Ku = (lo + hi) / 2

    ys = _p_loop_response(plant, Ku, setpoint, dt, steps)
    pk = _peaks(ys)
    Tu = (pk[-1][0] - pk[0][0]) * dt / max(1, len(pk) - 1)   # 平均峰间隔
    return Ku, Tu


def ziegler_nichols(Ku, Tu):
    """返回 (P, PI, PID) 三组参数。每组为 (Kp, Ki, Kd)。"""
    p = (0.5 * Ku, 0.0, 0.0)
    pi = (0.45 * Ku, 0.45 * Ku / (0.83 * Tu), 0.0)
    pid = (0.6 * Ku, 0.6 * Ku / (0.5 * Tu), 0.6 * Ku * (0.125 * Tu))
    return p, pi, pid


def main():
    plant = tf_to_ss([1.0], [1.0, 3.0, 2.0, 0.0])   # 1/(s(s+1)(s+2))
    Ku, Tu = find_ultimate(plant)
    print("临界增益 Ku = %.3f, 临界周期 Tu = %.3f s" % (Ku, Tu))
    p, pi, pid = ziegler_nichols(Ku, Tu)
    print("Z-N 整定参数 (Kp, Ki, Kd):")
    print("  P   : %.3f, %.3f, %.3f" % p)
    print("  PI  : %.3f, %.3f, %.3f" % pi)
    print("  PID : %.3f, %.3f, %.3f" % pid)


if __name__ == "__main__":
    main()
