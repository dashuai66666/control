"""Bode 图（Bode Plot）。

用途
----
用两张对数坐标图描述系统的**频率响应**：幅频图（增益 dB 对频率）与相频图
（相位对频率）。它是频域分析与控制器设计的核心工具，可直观读出带宽、
谐振峰、高频滚降等特性。

数学原理
--------
令 s = jω 代入传递函数 G(s)，得到频率响应：

    幅值  |G(jω)|，用分贝表示：20·log10|G(jω)|
    相位  ∠G(jω) = atan2(Im, Re)

频率轴取对数，因此"每十倍频程下降 -20dB"对应一个积分环节等规律一眼可见。

本实现输出文本表格；生产环境可用 matplotlib 绘制平滑曲线。
"""

import math

from _control import tf_response


def logspace(a, b, n):
    """10^a 到 10^b 之间取 n 个对数均匀点。"""
    return [10 ** (a + (b - a) * i / (n - 1)) for i in range(n)]


def bode(num, den, ws):
    """返回 (magnitudes_dB, phases_deg) 列表。"""
    mags, phases = [], []
    for w in ws:
        g = tf_response(num, den, w)
        mags.append(20 * math.log10(abs(g) + 1e-12))
        phases.append(math.degrees(math.atan2(g.imag, g.real)))
    return mags, phases


def main():
    # 一阶惯性环节 1/(s+1)：转折频率 ω=1，此处增益 -3dB、相位 -45°
    num, den = [1.0], [1.0, 1.0]
    ws = logspace(-2, 2, 9)
    mags, phases = bode(num, den, ws)

    print("G(s) = 1/(s+1) 的 Bode 数据:")
    print("%8s %10s %10s" % ("ω", "增益(dB)", "相位(°)"))
    for w, m, p in zip(ws, mags, phases):
        print("%8.3f %10.2f %10.1f" % (w, m, p))


if __name__ == "__main__":
    main()
