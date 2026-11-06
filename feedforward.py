"""前馈控制（Feedforward Control）。

用途
----
反馈控制是"看到偏差再纠正"，天然滞后；前馈控制则根据参考输入（或可测扰动）
**直接**按对象模型算出所需的控制量，做到"未偏差先动作"。两者结合（前馈 +
反馈）既能快速跟踪设定值，又能靠反馈抑制扰动与模型误差。

数学原理
--------
若对象稳态增益为 K，要使其输出 y 达到设定值 r，稳态所需输入即 u = r/K。
把这一项作为前馈直接叠加到反馈量上：

    u = u_ff + u_fb = (r/K) + Kp·(r - y)

反馈项只负责修正误差，前馈项负责"主要驱动"，从而显著加快响应并消除稳态误差。
"""


def simulate_first_order(tau, u_seq, dt):
    """一阶对象 y' = (u - y)/τ 的开环仿真。"""
    y = 0.0
    ys = []
    for u in u_seq:
        y += dt * (u - y) / tau
        ys.append(y)
    return ys


def main():
    tau, dt, r, Kp = 1.0, 0.01, 1.0, 5.0
    steps = 800

    # 方案1：纯 P 反馈
    y1, u = 0.0, 0.0
    out1 = []
    for _ in range(steps):
        u = Kp * (r - y1)
        y1 += dt * (u - y1) / tau
        out1.append(y1)

    # 方案2：前馈 + P 反馈（前馈用稳态逆模型 u_ff = r/K = r）
    y2 = 0.0
    out2 = []
    for _ in range(steps):
        u = r + Kp * (r - y2)                 # u_ff = r（K=1）
        y2 += dt * (u - y2) / tau
        out2.append(y2)

    print("阶跃设定值 r=1，纯 P 反馈 vs 前馈+P 反馈:")
    print("  纯 P 反馈   : 稳态值 %.3f（存在稳态误差）" % out1[-1])
    print("  前馈+P 反馈 : 稳态值 %.3f（无稳态误差）" % out2[-1])


if __name__ == "__main__":
    main()
