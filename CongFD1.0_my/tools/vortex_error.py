#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
三维等熵涡对流基准：数值解（CGNS）与解析解逐点对比 → L∞/L² 误差与网格收敛阶
（任务 8.5-外验；仅依赖系统 h5dump，无需 h5py/numpy）

用法:
  python3 tools/vortex_error.py <t> <N1> <f1.cgns> [<N2> <f2.cgns> ...]
  例: python3 tools/vortex_error.py 2.0 16 a.cgns 32 b.cgns 64 c.cgns

约定（与 src/src/initializer.cpp 的 nCase=2 初值一致）:
  周期域 [0,2pi]^3；涡心=(pi,pi,pi)；涡轴=平流方向=(1,1,1)/sqrt(3)；
  eps=5(涡强度)；基态 rho=p=1；解析解 = 初值沿 U*t 刚体平移: q(x,t)=q0(x-Ut)
  格点 = 格心 (i+0.5)h，i=0..N-1，h=2pi/N；CGNS 数组序 idx=i+j*N+k*N*N
"""
import subprocess, re, sys, math

GAMMA = 1.4
EPS   = 5.0
PI    = math.pi
N3    = 1.0 / math.sqrt(3.0)
UB    = (N3, N3, N3)

DS_PATHS = {
    "rho": "/Base/Zone  1/flowSolution/Density/ data",
    "u":   "/Base/Zone  1/flowSolution/XVelocity/ data",
    "v":   "/Base/Zone  1/flowSolution/YVelocity/ data",
    "w":   "/Base/Zone  1/flowSolution/ZVelocity/ data",
    "p":   "/Base/Zone  1/flowSolution/Pressure/ data",
}
KEYS = ["rho", "u", "v", "w", "p"]
NUMRE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")

def read_array(path, dspath):
    # -m %.17g：double 全精度输出（默认 %g 仅 6 位有效数字，会成为误差地板）
    out = subprocess.run(["h5dump", "-y", "-m", "%.17g", "-d", dspath, path],
                         capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError("h5dump 失败(%s): %s" % (dspath, out.stderr.strip()))
    body = out.stdout.split("DATA {", 1)[1].rsplit("}", 1)[0]
    return [float(v) for v in NUMRE.findall(body)]

def exact(x, y, z, t):
    """解析解：初值场沿 UB*t 平移"""
    dx = x - PI - UB[0] * t
    dy = y - PI - UB[1] * t
    dz = z - PI - UB[2] * t
    dn = (dx + dy + dz) * N3
    r2 = dx * dx + dy * dy + dz * dz - dn * dn
    g  = math.exp((1.0 - r2) * 0.5)
    amp = EPS / (2.0 * PI) * g
    u = UB[0] + amp * N3 * (dz - dy)
    v = UB[1] + amp * N3 * (dx - dz)
    w = UB[2] + amp * N3 * (dy - dx)
    T = 1.0 - (GAMMA - 1.0) * EPS * EPS / (8.0 * GAMMA * PI * PI) * math.exp(1.0 - r2)
    rho = T ** (1.0 / (GAMMA - 1.0))
    p   = T ** (GAMMA / (GAMMA - 1.0))
    return rho, u, v, w, p

def errors(path, N, t, order2=False):
    A = {k: read_array(path, v) for k, v in DS_PATHS.items()}
    nN = N ** 3
    if len(A["rho"]) != nN:
        raise RuntimeError("点数不符: %d vs %d" % (len(A["rho"]), nN))
    h = 2.0 * PI / N
    linf = {k: 0.0 for k in KEYS}
    l2   = {k: 0.0 for k in KEYS}
    for k in range(N):
        for j in range(N):
            for i in range(N):
                idx = i + j * N + k * N * N
                if order2:
                    idx = k + j * N + i * N * N
                x = (i + 0.5) * h
                y = (j + 0.5) * h
                z = (k + 0.5) * h
                ex = exact(x, y, z, t)
                for m, key in enumerate(KEYS):
                    d = abs(A[key][idx] - ex[m])
                    if d > linf[key]:
                        linf[key] = d
                    l2[key] += d * d
    for key in KEYS:
        l2[key] = math.sqrt(l2[key] / nN)
    return linf, l2

def main():
    args = sys.argv[1:]
    order2 = "--order2" in args
    args = [a for a in args if a != "--order2"]
    t = float(args[0])
    cases = [(int(args[i]), args[i + 1]) for i in range(1, len(args), 2)]
    prev = None
    print("== 三维等熵涡对流：数值 vs 解析（t=%.4g，order%s）" % (t, "=k-fast" if order2 else "=i-fast"))
    for N, path in cases:
        linf, l2 = errors(path, N, t, order2)
        print("N=%d  h=%.4e" % (N, 2.0 * PI / N))
        for key in KEYS:
            print("  %-3s Linf=%.3e  L2=%.3e" % (key, linf[key], l2[key]))
        if prev is not None:
            pN, plinf, pl2 = prev
            print("  收敛阶 (N=%d->%d): " % (pN, N)
                  + "  ".join("%s: Linf %.2f / L2 %.2f" %
                              (key, math.log2(plinf[key] / linf[key]),
                               math.log2(pl2[key] / l2[key])) for key in KEYS))
        prev = (N, linf, l2)
    print("== 完成")

if __name__ == "__main__":
    main()
