"""二阶系统时域响应（Second-Order System）。

用途
----
二阶系统是最重要的动态模型：弹簧-质量-阻尼、RLC 电路、伺服机构等都可化为
标准二阶形式。通过阻尼比 ζ 可判断系统是欠阻尼（振荡）、临界阻尼还是过阻尼。

数学模型
--------
传递函数：G(s) = ωn² / (s² + 2·ζ·ωn·s + ωn²)

- ωn：无阻尼自然频率，决定响应的"快慢"；
- ζ：阻尼比，决定"振荡程度"（0<ζ<1 欠阻尼振荡；ζ≥1 无振荡）。

欠阻尼（0<ζ<1）阶跃响应的关键性能指标
-------------------------------------
- 超调量   Mp = exp(-ζπ/√(1-ζ²))；
- 峰值时间 tp = π / (ωn·√(1-ζ²))；
- 调节时间 ts ≈ 4/(ζ·ωn)（2% 误差带）；
- 振荡频率 ωd = ωn·√(1-ζ²)。
"""

import math


def second_order_metrics(wn, zeta):
    """返回欠阻尼二阶系统的 (超调量, 峰值时间, 调节时间, 振荡频率)。"""
    if not (0.0 < zeta < 1.0):
        raise ValueError("该公式适用于欠阻尼 0 < ζ < 1")
    wd = wn * math.sqrt(1 - zeta * zeta)          # 阻尼振荡频率
    overshoot = math.exp(-zeta * math.pi / math.sqrt(1 - zeta * zeta))
    tp = math.pi / wd
    ts = 4.0 / (zeta * wn)                        # 2% 判据
    return overshoot, tp, ts, wd


def step_response(wn, zeta, t):
    """欠阻尼二阶系统单位阶跃响应解析解。"""
    wd = wn * math.sqrt(1 - zeta * zeta)
    phi = math.acos(zeta)
    return 1.0 - math.exp(-zeta * wn * t) / math.sqrt(1 - zeta * zeta) * math.sin(wd * t + phi)


def main():
    wn, zeta = 2.0, 0.5
    Mp, tp, ts, wd = second_order_metrics(wn, zeta)
    print("二阶系统 ωn=%.1f, ζ=%.2f:" % (wn, zeta))
    print("  超调量 Mp = %.1f%%" % (100 * Mp))
    print("  峰值时间 tp = %.3f s" % tp)
    print("  调节时间 ts = %.3f s（2%% 判据）" % ts)
    print("  振荡频率 ωd = %.3f rad/s" % wd)
    print("  阶跃响应采样: ", ["t=%.1f:%.3f" % (t, step_response(wn, zeta, t))
                              for t in (0.5, tp, 1.5, 2.5)])


if __name__ == "__main__":
    main()
