# ControlTheory

经典控制（含现代控制基础）算法集合（纯 Python，仅标准库）。与 ``DataAnalysis``、
``SignalAnalysis`` 同风格：每个文件实现一个典型算法，含详细中文注释与可运行
示例，覆盖时域分析、PID 及其整定、稳定性判据、频率法、根轨迹、状态空间、
最优控制、状态估计、校正器与离散化。

## 目录结构

```
ControlTheory/
├── _control.py                    # 共用工具：线性代数/多项式求根/传函/状态空间
├── first_order_system.py          # 一阶系统时域响应
├── second_order_system.py         # 二阶系统时域响应（超调/峰值/调节时间）
├── pid_controller.py              # PID 控制器
├── ziegler_nichols.py             # Ziegler-Nichols 临界比例度整定
├── routh_hurwitz.py               # Routh-Hurwitz 稳定性判据
├── bode_plot.py                   # Bode 图
├── nyquist.py                     # Nyquist 图与稳定性判据
├── root_locus.py                  # 根轨迹
├── gain_phase_margin.py           # 增益/相位裕度
├── state_space.py                 # 状态空间表示与仿真
├── controllability_observability.py # 能控性与能观性
├── pole_placement.py              # 极点配置（Ackermann 公式）
├── lqr.py                         # 线性二次型最优控制 LQR
├── kalman_filter.py               # 卡尔曼滤波（状态估计）
├── lead_lag_compensator.py        # 超前/滞后校正
├── discretization.py              # 连续系统离散化（Euler/Tustin/ZOH）
├── feedforward.py                 # 前馈控制
└── smith_predictor.py             # 史密斯预估器（纯滞后补偿）
```

## 算法分类

| 类别 | 文件 |
|------|------|
| 时域分析 | `first_order_system`、`second_order_system` |
| PID 与整定 | `pid_controller`、`ziegler_nichols` |
| 稳定性/频域 | `routh_hurwitz`、`bode_plot`、`nyquist`、`gain_phase_margin`、`root_locus` |
| 现代控制 | `state_space`、`controllability_observability`、`pole_placement`、`lqr`、`kalman_filter` |
| 校正/离散 | `lead_lag_compensator`、`discretization`、`feedforward`、`smith_predictor` |

## 运行

所有文件相互独立、可直接运行（在 `ControlTheory` 目录下）：

```bash
cd ControlTheory
python3 pid_controller.py           # 或 python3 -m pid_controller
python3 bode_plot.py
# ... 每个文件同理
```

每个文件末尾的 `main()` 会用典型对象演示算法效果。

## 说明

- 频率响应用 Python 内置 `complex` 计算；多项式求根用 Durand-Kerner 迭代；
- 矩阵运算、传递函数→状态空间、状态空间仿真集中在 `_control.py` 复用；
- 面向教学与中小规模系统；生产环境请使用 `scipy.signal`、`control`、`python-control`
  等成熟库。
