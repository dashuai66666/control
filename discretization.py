"""连续系统的离散化（Discretization）。

用途
----
数字控制器运行在离散时间上，需要把连续传递函数 G(s) 变成离散传递函数 G(z)。
常见的映射方法有前向欧拉、后向欧拉、双线性（Tustin）变换，以及精确的
零阶保持（ZOH）离散化。

数学原理
--------
把 G(s) 中的 s 用 z 的有理函数替换：

- 前向欧拉：  s = (z - 1) / T
- 后向欧拉：  s = (z - 1) / (T·z)
- Tustin（双线性）： s = (2/T)·(z - 1)/(z + 1)   ← 稳定域映射到单位圆内

ZOH 则精确描述"采样保持"：先把 s 域系统转成状态空间，再算
Ad = e^(A·T)、Bd = (Ad - I)·A⁻¹·B，得到与连续系统在采样点处完全一致的离散模型。

本实现用多项式替换实现三种代数方法，用矩阵指数实现 ZOH。
"""

from _control import poly_mul, poly_add


def _poly_pow(p, k):
    out = [1.0]
    for _ in range(k):
        out = poly_mul(out, p)
    return out


def _substitute(num, den, s_num, s_den):
    """G(s)=num(s)/den(s)，令 s = s_num(z)/s_den(z)，返回离散 (N, D) 多项式。

    分子分母各自替换后，还需通分：分子乘 s_den^m、分母乘 s_den^n。
    """
    n, m = len(num) - 1, len(den) - 1
    N = [0.0]
    for i in range(n + 1):
        term = poly_mul([num[n - i]],
                        poly_mul(_poly_pow(s_num, i), _poly_pow(s_den, n - i)))
        N = poly_add(N, term)
    D = [0.0]
    for i in range(m + 1):
        term = poly_mul([den[m - i]],
                        poly_mul(_poly_pow(s_num, i), _poly_pow(s_den, m - i)))
        D = poly_add(D, term)
    # 通分：G(z) = [N·s_den^m] / [D·s_den^n]
    N = poly_mul(N, _poly_pow(s_den, m))
    D = poly_mul(D, _poly_pow(s_den, n))
    return N, D


def _normalize(N, D):
    """把分母归一化为首一，分子同除，返回降幂系数。"""
    lead = D[0]
    return [c / lead for c in N], [c / lead for c in D]


def forward_euler(num, den, T):
    s_num, s_den = [1.0 / T, -1.0 / T], [1.0]
    return _normalize(*_substitute(num, den, s_num, s_den))


def backward_euler(num, den, T):
    s_num, s_den = [1.0, -1.0], [T, 0.0]
    return _normalize(*_substitute(num, den, s_num, s_den))


def tustin(num, den, T):
    s_num, s_den = [2.0 / T, -2.0 / T], [1.0, 1.0]
    return _normalize(*_substitute(num, den, s_num, s_den))


def zoh(num, den, T):
    """ZOH 精确离散化，返回状态空间 (Ad, Bd, Cd, Dd)。"""
    from _control import tf_to_ss, matrix_exp, inverse, identity, matmul, mat_sub
    A, B, C, D = tf_to_ss(num, den)
    n = len(A)
    Ad = matrix_exp(A, T)
    Ai = inverse(A)
    Bd = matmul(mat_sub(Ad, identity(n)), matmul(Ai, B))
    return Ad, Bd, C, D


def _fmt(poly):
    """把降幂系数列表格式化成字符串。"""
    terms = []
    n = len(poly) - 1
    for i, c in enumerate(poly):
        if abs(c) < 1e-12:
            continue
        power = n - i
        if power == 0:
            terms.append("%.4f" % c)
        elif power == 1:
            terms.append("%.4f z" % c)
        else:
            terms.append("%.4f z^%d" % (c, power))
    return " + ".join(terms) if terms else "0"


def main():
    num, den = [1.0], [1.0, 1.0]           # 1/(s+1)，时间常数 1s
    T = 0.1

    for name, (N, D) in [("前向欧拉", forward_euler(num, den, T)),
                         ("后向欧拉", backward_euler(num, den, T)),
                         ("Tustin", tustin(num, den, T))]:
        print("%s: G(z) = (%s) / (%s)" % (name, _fmt(N), _fmt(D)))

    Ad, Bd, Cd, Dd = zoh(num, den, T)
    print("ZOH 离散状态空间: Ad =", ["%.4f" % v for v in Ad[0]],
          "Bd =", ["%.4f" % v for v in Bd[0]])


if __name__ == "__main__":
    main()
