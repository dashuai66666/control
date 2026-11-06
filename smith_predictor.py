"""史密斯预估器（Smith Predictor）。

用途
----
控制含**纯滞后**的对象（如化工、传输过程）时，滞后会给反馈回路引入大相位
滞后，限制控制器增益、易引起振荡。史密斯预估器用一个对象模型"预测"未滞后的
输出，用预测值做反馈，从而把滞后"移出"闭环，使控制器可按无滞后对象设计。

数学原理
--------
对象模型 = G0(s)·e^(-θs)（G0 为无滞后部分，θ 为纯滞后时间）。控制器看到的是
模型的**未滞后输出** ŷ，因此闭环特征方程不含 e^(-θs)，滞后不再影响稳定性：

    u = K·(r - ŷ)          （ŷ 为模型的无滞后预测输出）

真实输出 y 只是 ŷ 的延迟版本，最终同样跟踪设定值，只是慢 θ 秒。

本实现仿真一个一阶 + 纯滞后对象，对比"直接 P 控制"与"史密斯预估器"在相同
高增益下的表现。
"""


def simulate(tau, delay_steps, Kp, setpoint, dt, steps, use_smith):
    """一阶 + 纯滞后对象的闭环仿真，返回输出序列。

    use_smith=True 时用无滞后模型输出做反馈（史密斯预估器）。
    """
    ym = 0.0                      # 无滞后模型输出（也是"预测"输出）
    buffer = [0.0] * (delay_steps + 1)   # 纯滞后缓冲区
    ys = []
    for _ in range(steps):
        if use_smith:
            u = Kp * (setpoint - ym)     # 反馈预测的未滞后输出
        else:
            u = Kp * (setpoint - buffer[0])   # 反馈真实（滞后）输出
        ym += dt * (u - ym) / tau        # 模型动态
        buffer = buffer[1:] + [ym]       # 延迟 θ
        ys.append(buffer[0])
    return ys


def main():
    tau, theta = 1.0, 2.0        # 时间常数 1s，纯滞后 2s
    dt = 0.01
    delay_steps = int(theta / dt)
    Kp = 3.0                     # 高增益（直接控制会因滞后而振荡）
    steps = 3000

    y_direct = simulate(tau, delay_steps, Kp, 1.0, dt, steps, use_smith=False)
    y_smith = simulate(tau, delay_steps, Kp, 1.0, dt, steps, use_smith=True)

    # 用输出波动幅度衡量稳定性
    def swing(ys):
        tail = ys[-500:]
        return max(tail) - min(tail)

    print("一阶 + 纯滞后对象（τ=1s, θ=2s），增益 Kp=%.1f:" % Kp)
    print("  直接 P 控制  : 稳态值 %.3f, 尾段波动 %.3f" % (y_direct[-1], swing(y_direct)))
    print("  史密斯预估器  : 稳态值 %.3f, 尾段波动 %.3f" % (y_smith[-1], swing(y_smith)))
    print("  （纯 P 控制存在 1/(1+Kp) 稳态误差；史密斯预估器的关键是消除滞后引起的振荡）")


if __name__ == "__main__":
    main()
