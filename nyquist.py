"""Nyquist 图与稳定性判据（Nyquist Plot & Criterion）。

用途
----
在复平面上画出开环频率响应 G(jω)（ω 从 -∞ 到 +∞），根据曲线对 (-1, 0) 点的
**包围情况**判断闭环系统稳定性。Nyquist 判据比 Bode/Routh 更一般，可直接处理
含纯滞后的系统。

数学原理（简化版，开环稳定时）
------------------------------
若开环 G(s) 稳定（极点都在左半平面），则单位负反馈闭环稳定 ⇔ Nyquist 曲线
**不包围** (-1, 0) 点。包围圈数 = 曲线逆时针穿过"负实轴上 -1 左侧"的净次数。

本实现通过扫描 ω ∈ [0, ωmax]，利用实系数系统的对称性（G(-jω)=conj(G(jω))），
统计曲线与 (-∞, -1) 段的净穿越次数来估计包围圈数。
"""

import math

from _control import tf_response
from bode_plot import logspace


def nyquist_crossings(num, den, wmin=-2, wmax=2, n=2000):
    """返回 Nyquist 曲线对 (-∞,-1) 段的净穿越次数（近似包围圈数）。"""
    ws = logspace(wmin, wmax, n)
    crossings = 0
    prev = tf_response(num, den, ws[0])
    for w in ws[1:]:
        cur = tf_response(num, den, w)
        # 虚部变号 → 穿过实轴；且实部 < -1 → 在 (-1,0) 左侧
        if prev.imag * cur.imag < 0 and prev.real < -1.0 and cur.real < -1.0:
            crossings += 1 if (prev.imag < 0 and cur.imag > 0) else -1
        prev = cur
    # 实系数系统曲线关于实轴对称，[0,∞) 的穿越圈数即总圈数
    return crossings


def main():
    # 例1：G=1/((s+1)(s+2))，Nyquist 不包围 -1 → 闭环稳定
    g1 = ([1.0], [1.0, 3.0, 2.0])
    # 例2：G=10/(s+1)^3（临界增益 K=8），包围 -1 → 闭环不稳定
    g2 = ([10.0], [1.0, 3.0, 3.0, 1.0])

    for name, (num, den) in [("1/((s+1)(s+2))", g1), ("10/(s+1)^3", g2)]:
        n = nyquist_crossings(num, den)
        print("%s: 包围(-1,0)圈数 = %d -> 闭环%s"
              % (name, n, "稳定" if n == 0 else "不稳定"))


if __name__ == "__main__":
    main()
