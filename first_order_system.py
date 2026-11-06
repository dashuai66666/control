"""一阶系统时域响应（First-Order System）。

用途
----
大量实际对象（RC 电路、温度过程、水箱液位等）可近似为一阶惯性环节。分析其
阶跃响应是理解系统动态特性的起点。

数学模型
--------
传递函数：G(s) = K / (τ·s + 1)

- K：稳态增益，决定最终稳态值；
- τ：时间常数，决定响应快慢。

单位阶跃输入下的响应为：

    y(t) = K·(1 - e^(-t/τ))

关键指标
--------
- t = τ 时达到稳态值的 63.2%；
- t = 4τ 时达到约 98.2%（常作为"进入稳态"的经验判据）；
- 上升时间（10%→90%）约 2.2τ。
"""

import math


def first_order_step(K, tau, t):
    """单位阶跃输入下一阶系统的输出 y(t) = K(1 - e^{-t/τ})。"""
    return K * (1.0 - math.exp(-t / tau))


def main():
    K, tau = 2.0, 1.5
    print("一阶系统 G(s) = %.1f/(%.1f·s + 1) 的阶跃响应:" % (K, tau))
    for t in [0.0, tau, 2 * tau, 3 * tau, 4 * tau, 10 * tau]:
        print("  t=%5.2fs  y=%.4f  (稳态的 %.1f%%)" % (t, first_order_step(K, tau, t),
                                                    100 * first_order_step(K, tau, t) / K))


if __name__ == "__main__":
    main()
