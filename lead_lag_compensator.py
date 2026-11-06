"""超前/滞后校正（Lead/Lag Compensator）。

用途
----
在频域里"塑造"开环频率响应以改善稳定裕度：

- **超前校正（lead）**：z < p，在某频率提供相位超前，用于增加相位裕度、
  加快响应（相当于微分作用）；
- **滞后校正（lag）**：z > p，压低高频增益、基本不引入相位滞后（或少量滞后），
  用于提高低频增益、减小稳态误差而不损害稳定性。

数学原理（超前校正）
--------------------
对一阶校正器 Gc(s) = (s+z)/(s+p)，其最大相位超前 φ_max 出现在几何平均频率
ωm = √(z·p)，且满足：

    sin(φ_max) = (p - z) / (p + z)

给定所需 φ_max 与作用频率 ωm，可反解出：

    α = p/z = (1 + sin φ_max) / (1 - sin φ_max)
    z = ωm / √α,   p = ωm · √α

本实现据此设计校正器，并归一化使其低频增益为 1。
"""

import math

from _control import tf_response


def lead_compensator(phi_deg, wm):
    """在频率 wm 处提供最大相位超前 phi_deg 的超前校正器 (num, den)。"""
    phi = math.radians(phi_deg)
    alpha = (1.0 + math.sin(phi)) / (1.0 - math.sin(phi))
    z = wm / math.sqrt(alpha)
    p = wm * math.sqrt(alpha)
    # (p/z)·(s+z)/(s+p)，低频增益归一化为 1
    return [p / z, p / z * z], [1.0, p]


def lag_compensator(z, p):
    """滞后校正器 (num, den)，z > p，低频增益 1。"""
    return [p / z, p / z * z], [1.0, p]


def main():
    # 在 ωm=2 处提供 45° 相位超前的超前校正器
    num, den = lead_compensator(45.0, 2.0)
    print("超前校正器 Gc(s) = (%.3f s + %.3f)/(s + %.3f)" % (num[0], num[1], den[1]))
    g = tf_response(num, den, 2.0)
    print("  在 ω=2 处提供相位 %.1f°（期望约 45°）" % math.degrees(math.atan2(g.imag, g.real)))


if __name__ == "__main__":
    main()
