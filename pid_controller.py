"""PID 控制器（Proportional-Integral-Derivative）。

用途
----
PID 是工业控制中使用最广泛的反馈控制器：根据误差 e = 设定值 - 输出，
按比例、积分、微分三项加权产生控制量：

    u(t) = Kp·e + Ki·∫e dt + Kd·(de/dt)

- 比例项 P：对当前误差立即反应，增大可加快响应但易超调/震荡；
- 积分项 I：消除稳态误差（残差），但可能引起积分饱和与慢振荡；
- 微分项 D：预测误差趋势，提供阻尼、抑制超调，但对噪声敏感。

本实现用离散化形式（位置式）仿真一阶对象在阶跃设定值下的闭环跟踪。
"""


class PID:
    """位置式离散 PID。"""

    def __init__(self, kp, ki, kd, dt):
        self.kp, self.ki, self.kd, self.dt = kp, ki, kd, dt
        self.integral = 0.0
        self.prev_e = 0.0

    def update(self, e):
        self.integral += e * self.dt               # 积分累加
        deriv = (e - self.prev_e) / self.dt        # 微分（差分近似）
        self.prev_e = e
        return self.kp * e + self.ki * self.integral + self.kd * deriv


def simulate_first_order_plant(tau, setpoint, dt, steps, pid):
    """仿真一阶对象 y' = (u - y)/τ 在 PID 控制下的闭环响应。"""
    y = 0.0
    ys = []
    for _ in range(steps):
        e = setpoint - y
        u = pid.update(e)
        y += dt * (u - y) / tau                    # 欧拉积分对象
        ys.append(y)
    return ys


def main():
    tau = 1.0
    dt = 0.01
    setpoint = 1.0

    # 三组参数演示 P / PI / PID 的差别
    configs = [
        ("纯 P (Kp=2)", dict(kp=2.0, ki=0.0, kd=0.0)),
        ("PI  (Kp=2,Ki=0.8)", dict(kp=2.0, ki=0.8, kd=0.0)),
        ("PID (Kp=2,Ki=0.8,Kd=0.1)", dict(kp=2.0, ki=0.8, kd=0.1)),
    ]
    for name, c in configs:
        pid = PID(dt=dt, **c)
        ys = simulate_first_order_plant(tau, setpoint, dt, steps=1000, pid=pid)
        final = ys[-1]
        overshoot = max(ys) - setpoint
        print("%-22s 稳态值=%.3f 超调=%.3f" % (name, final, max(0.0, overshoot)))


if __name__ == "__main__":
    main()
