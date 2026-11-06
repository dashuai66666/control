"""能控性与能观性（Controllability & Observability）。

用途
----
- **能控性**：能否通过输入 u 把系统状态从任意初始值驱动到任意目标值；
- **能观性**：能否仅通过输出 y 的观测反推出系统的全部内部状态。

它们是状态反馈控制与状态观测器设计的前提条件。

数学判据（秩判据）
------------------
- 能控性矩阵 Wc = [B, AB, A²B, ..., Aⁿ⁻¹B]，系统能控 ⇔ rank(Wc) = n；
- 能观性矩阵 Wo = [C; CA; CA²; ...; CAⁿ⁻¹]，系统能观 ⇔ rank(Wo) = n。

本实现用高斯消元求秩（见 ``_control.rank``）。
"""

from _control import rank, controllability_matrix, observability_matrix


def main():
    # 例1：能控且能观
    A1 = [[0.0, 1.0], [-2.0, -3.0]]
    B1 = [[0.0], [1.0]]
    C1 = [[1.0, 0.0]]

    # 例2：不能控（第二个模态不可达，B 只作用在第一个状态）
    A2 = [[1.0, 0.0], [0.0, 2.0]]
    B2 = [[1.0], [0.0]]
    C2 = [[1.0, 1.0]]

    for name, A, B, C in [("例1", A1, B1, C1), ("例2", A2, B2, C2)]:
        Wc = controllability_matrix(A, B)
        Wo = observability_matrix(A, C)
        print("%s: rank(Wc)=%d/%d %s, rank(Wo)=%d/%d %s"
              % (name, rank(Wc), len(A), "能控" if rank(Wc) == len(A) else "不能控",
                 rank(Wo), len(A), "能观" if rank(Wo) == len(A) else "不能观"))


if __name__ == "__main__":
    main()
