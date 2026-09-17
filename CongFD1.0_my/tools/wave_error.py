#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
三维熵波基准（nCase=3）：数值解（CGNS）vs 解析解逐点对比 → Linf/L2 误差与收敛阶
（任务 8 外部真值验证的查证脚本；只读 runs/wav2_n* 下已生成的 CGNS）

用法:
  cd <repo> && python3 tools/wave_error.py [t(默认2.0)]

设置（与 src/src/initializer.cpp nCase=3 一致）:
  rho = 1 + 0.2 sin(x+y+z)，u=v=w=1/sqrt(3)，p=1（波矢 (1,1,1)，2pi 域严格周期）
  解析解: rho(x,t) = 1 + 0.2 sin(x+y+z - sqrt(3) t)（欧拉方程严格精确解）
"""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vortex_error import read_array, DS_PATHS, PI

U = 1.0 / math.sqrt(3.0)          # 平流速度分量；n·U = 3/sqrt(3) = sqrt(3)
AMP = 0.2

def errs(iMax, t):
    dirn = "runs/wav2_n%d" % iMax
    f = "%s/3D_wave - TENO-AS-Ff_5_10 - %dx%dx%d - t=%.4f.cgns" % (dirn, iMax, iMax, iMax, t)
    rho = read_array(f, DS_PATHS["rho"])
    N = iMax - 1
    h = 2.0 * PI / N
    linf = 0.0; l2 = 0.0
    for k in range(N):
        for j in range(N):
            for i in range(N):
                x = (i + .5) * h; y = (j + .5) * h; z = (k + .5) * h
                ex = 1.0 + AMP * math.sin(x + y + z - 3.0 * U * t)
                d = abs(rho[i + j * N + k * N * N] - ex)
                if d > linf:
                    linf = d
                l2 += d * d
    return linf, math.sqrt(l2 / N ** 3)

def main():
    t = float(sys.argv[1]) if len(sys.argv) > 1 else 2.0
    print("== 三维熵波 vs 解析解（t=%.4g）" % t)
    prev = None
    for iMax in (17, 33, 65):
        if not os.path.isdir("runs/wav2_n%d" % iMax):
            print("  (缺 runs/wav2_n%d，跳过)" % iMax); continue
        linf, l2 = errs(iMax, t)
        print("  iMax=%3d (格 %2d): Linf=%.3e  L2=%.3e" % (iMax, iMax - 1, linf, l2))
        if prev is not None:
            print("    收敛阶 (%d->%d 格): Linf %.2f / L2 %.2f"
                  % (prev[0] - 1, iMax - 1, math.log2(prev[1] / linf), math.log2(prev[2] / l2)))
        prev = (iMax, linf, l2)

if __name__ == "__main__":
    main()
