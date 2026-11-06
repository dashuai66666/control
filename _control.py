"""经典控制算法共用的纯 Python 工具：线性代数、多项式、传递函数与状态空间。

这些都是"基础件"而非算法本身。控制算法文件复用这里的矩阵运算、多项式求根、
频率响应与状态空间仿真，避免重复实现。全部基于嵌套 list，不依赖 numpy。
"""

import math


# ---------------------------------------------------------------------------
# 线性代数
# ---------------------------------------------------------------------------

def zeros(r, c):
    return [[0.0] * c for _ in range(r)]


def identity(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]


def matmul(A, B):
    """矩阵乘法 A (m,k) × B (k,n) -> (m,n)。"""
    m, k = len(A), len(A[0])
    n = len(B[0])
    out = zeros(m, n)
    for i in range(m):
        for j in range(n):
            out[i][j] = sum(A[i][t] * B[t][j] for t in range(k))
    return out


def matvec(A, x):
    return [sum(a * b for a, b in zip(row, x)) for row in A]


def transpose(A):
    return [list(col) for col in zip(*A)]


def mat_add(A, B):
    return [[a + b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def mat_sub(A, B):
    return [[a - b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def mat_scale(A, k):
    return [[a * k for a in row] for row in A]


def inverse(A):
    """高斯-约当消元求逆（方阵）。"""
    n = len(A)
    M = [list(row) + [1.0 if i == j else 0.0 for j in range(n)]
         for i, row in enumerate(A)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[pivot][col]) < 1e-12:
            raise ValueError("矩阵奇异，无法求逆")
        M[col], M[pivot] = M[pivot], M[col]
        pv = M[col][col]
        M[col] = [v / pv for v in M[col]]
        for r in range(n):
            if r != col:
                f = M[r][col]
                M[r] = [a - f * b for a, b in zip(M[r], M[col])]
    return [row[n:] for row in M]


def rank(A, tol=1e-9):
    """矩阵的秩（高斯消元）。"""
    M = [row[:] for row in A]
    m, n = len(M), len(M[0])
    r = 0
    for col in range(n):
        pivot = next((i for i in range(r, m) if abs(M[i][col]) > tol), None)
        if pivot is None:
            continue
        M[r], M[pivot] = M[pivot], M[r]
        pv = M[r][col]
        M[r] = [v / pv for v in M[r]]
        for i in range(m):
            if i != r and abs(M[i][col]) > tol:
                f = M[i][col]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        r += 1
    return r


def matrix_exp(A, dt, terms=80):
    """矩阵指数 e^(A·dt)，用级数展开（dt 较小时收敛很快）。"""
    n = len(A)
    P = mat_scale(A, dt)
    S = identity(n)
    term = identity(n)
    for k in range(1, terms):
        term = mat_scale(matmul(term, P), 1.0 / k)
        S = mat_add(S, term)
    return S


def matrix_power(A, k):
    """方阵的 k 次幂（k=0 返回单位阵）。"""
    n = len(A)
    out = identity(n)
    for _ in range(k):
        out = matmul(out, A)
    return out


def controllability_matrix(A, B):
    """能控性矩阵 [B, AB, A²B, ..., A^{n-1}B]（B 为 n×1）。"""
    n = len(A)
    cols = []
    cur = [row[:] for row in B]
    for _ in range(n):
        cols.append([cur[i][0] for i in range(n)])
        cur = matmul(A, cur)
    return transpose(cols)            # n×n，第 k 列为 A^k B


def observability_matrix(A, C):
    """能观性矩阵 [C; CA; CA²; ...; CA^{n-1}]（C 为 1×n）。"""
    n = len(A)
    rows = []
    cur = [row[:] for row in C]
    for _ in range(n):
        rows.append(cur[0][:])
        cur = matmul(cur, A)
    return rows                       # n×n，第 k 行为 C A^k


# ---------------------------------------------------------------------------
# 多项式（系数按降幂排列）
# ---------------------------------------------------------------------------

def poly_eval(coeffs, s):
    """霍纳法则求多项式值，s 可为实数或复数。"""
    val = coeffs[0]
    for c in coeffs[1:]:
        val = val * s + c
    return val


def poly_mul(a, b):
    """多项式乘法（即系数卷积）。"""
    out = [0.0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out


def poly_add(a, b):
    """多项式加法（按降幂对齐）。"""
    la, lb = len(a), len(b)
    n = max(la, lb)
    pa = [0.0] * (n - la) + list(a)
    pb = [0.0] * (n - lb) + list(b)
    return [x + y for x, y in zip(pa, pb)]


def poly_roots(coeffs, iters=500):
    """用 Durand-Kerner 方法同时求多项式的所有复根。

    适用于根轨迹等需要多次求根的场景，无需推导闭式解。
    """
    coeffs = [c for c in coeffs]
    while len(coeffs) > 1 and abs(coeffs[0]) < 1e-12:
        coeffs = coeffs[1:]                      # 去掉前导零
    n = len(coeffs) - 1
    if n <= 0:
        return []
    a = coeffs[0]
    # 初始猜测散布在复平面，避免对称退化
    roots = []
    for k in range(n):
        ang = 2 * math.pi * k / n + 0.7
        rad = 0.4 + 0.9 * ((k * 0.61803398875) % 1.0)
        roots.append(complex(rad * math.cos(ang), rad * math.sin(ang)))
    for _ in range(iters):
        for i in range(n):
            p = poly_eval(coeffs, roots[i])
            denom = a
            for j in range(n):
                if j != i:
                    denom *= (roots[i] - roots[j])
            if abs(denom) < 1e-12:
                roots[i] += complex(1e-3, 1e-3)
                continue
            roots[i] -= p / denom
    return roots


# ---------------------------------------------------------------------------
# 传递函数 / 状态空间
# ---------------------------------------------------------------------------

def tf_response(num, den, w):
    """传递函数在 s=jω 处的频率响应（复数）。num/den 降幂。"""
    return poly_eval(num, 1j * w) / poly_eval(den, 1j * w)


def tf_to_ss(num, den):
    """传递函数 → 能控标准型状态空间 (A, B, C, D)。

    要求分母归一化为首一多项式，num 次数不超过 den。
    """
    lead = den[0]
    den = [c / lead for c in den]
    num = [c / lead for c in num]
    n = len(den) - 1
    num_pad = [0.0] * (n - len(num)) + list(num)      # 补齐到长度 n

    A = zeros(n, n)
    for j in range(n):
        A[0][j] = -den[j + 1]
    for i in range(1, n):
        A[i][i - 1] = 1.0
    B = [[1.0] if i == 0 else [0.0] for i in range(n)]   # n×1
    C = [num_pad]                                        # 1×n
    D = [[0.0]]
    return A, B, C, D


def ss_simulate(A, B, C, D, u, dt, x0=None):
    """状态空间前向欧拉仿真，返回 (状态轨迹, 输出轨迹)。u 为输入序列。"""
    n = len(A)
    x = [0.0] * n if x0 is None else list(x0)
    ys, xs = [], []
    for uk in u:
        Ax = [sum(A[i][j] * x[j] for j in range(n)) for i in range(n)]
        Bu = [B[i][0] * uk for i in range(n)]
        x = [x[i] + dt * (Ax[i] + Bu[i]) for i in range(n)]
        y = sum(C[0][j] * x[j] for j in range(n)) + D[0][0] * uk
        ys.append(y)
        xs.append(list(x))
    return xs, ys
