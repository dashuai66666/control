"""增益裕度与相位裕度（Gain Margin & Phase Margin）。

用途
----
量化一个反馈系统离"不稳定"还有多远，是频域设计中最常用的稳定裕度指标：

- **增益裕度 GM**：相位为 -180° 时，增益还能再增大多少分贝才到临界；
- **相位裕度 PM**：增益为 0dB 时，相位离 -180° 还差多少度。

经验上，PM ≥ 30°~60°、GM ≥ 6dB 通常认为有足够鲁棒性。

数学原理
--------
- 相位穿越频率 ω_pc：∠G(jω_pc) = -180°，则 GM = -20·log10|G(jω_pc)|；
- 增益穿越频率 ω_gc：|G(jω_gc)| = 1（0dB），则 PM = 180° + ∠G(jω_gc)。

本实现通过对数扫频 + 线性插值寻找穿越频率。
"""

import math

from _control import tf_response
from bode_plot import logspace


def find_crossing(ws, values, target):
    """在采样序列中找到 value 穿过 target 的位置（线性插值）。"""
    for i in range(len(ws) - 1):
        if (values[i] - target) * (values[i + 1] - target) < 0:
            t = (target - values[i]) / (values[i + 1] - values[i])
            return ws[i] + t * (ws[i + 1] - ws[i])
    return None


def margins(num, den):
    """返回 (gain_margin_dB, phase_margin_deg)。

    增益裕度用"G(jω) 的虚部过零且实部为负"来定位相位穿越频率（相位恰为 -180°），
    避免 atan2 相位在 -180° 处的跳变问题。
    """
    ws = logspace(-3, 3, 4000)
    mags_db = []
    for w in ws:
        g = tf_response(num, den, w)
        mags_db.append(20 * math.log10(abs(g) + 1e-12))

    # 增益裕度：找 Im(G)=0 且 Re<0 的穿越频率（相位穿越）
    w_pc = None
    prev = tf_response(num, den, ws[0])
    for i in range(1, len(ws)):
        cur = tf_response(num, den, ws[i])
        if prev.imag * cur.imag < 0 and cur.real < 0:
            t = prev.imag / (prev.imag - cur.imag)
            w_pc = ws[i - 1] + t * (ws[i] - ws[i - 1])
            break
        prev = cur
    gm = None
    if w_pc is not None:
        g = tf_response(num, den, w_pc)
        gm = -20 * math.log10(abs(g) + 1e-12)

    # 相位裕度：找 0dB 增益穿越频率
    w_gc = find_crossing(ws, mags_db, 0.0)
    pm = None
    if w_gc is not None:
        g = tf_response(num, den, w_gc)
        pm = 180.0 + math.degrees(math.atan2(g.imag, g.real))

    return gm, pm


def main():
    # 三阶对象 5/(s(s+1)(s+2))
    num, den = [5.0], [1.0, 3.0, 2.0, 0.0]
    gm, pm = margins(num, den)
    print("G(s) = 5/(s(s+1)(s+2)) 的稳定裕度:")
    print("  增益裕度 GM = %.2f dB" % (gm if gm is not None else float("inf")))
    print("  相位裕度 PM = %.2f°" % (pm if pm is not None else float("inf")))


if __name__ == "__main__":
    main()
