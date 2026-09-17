#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
任务 8.6：dim3CrossTest 输出比对——2D 全场 vs 3D 的 z 第一层（逐点），
并检验 3D 的 z 方向均匀性（k=0 与 k=1 层应完全一致）。

用法（在 dim3CrossTest 运行目录内）:
  cd runs/cross && python3 ../../tools/cross_compare.py
"""
import sys, os, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vortex_error import read_array

NX = 32          # 格数（iMax=33）
DS = {           # 变量名 -> CGNS 数据集路径（2D/3D 通用；2D 无 ZVelocity）
    "rho": "/Base/Zone  1/flowSolution/Density/ data",
    "u":   "/Base/Zone  1/flowSolution/XVelocity/ data",
    "v":   "/Base/Zone  1/flowSolution/YVelocity/ data",
    "p":   "/Base/Zone  1/flowSolution/Pressure/ data",
}

def main():
    f2 = glob.glob("*2D_TGV* - t=*.cgns")
    f3 = glob.glob("*3D_TGV2D* - t=*.cgns")
    if not f2 or not f3:
        print("未找到输出文件（先在 runs/cross 下运行 ./dim3CrossTest）"); sys.exit(2)
    f2, f3 = sorted(f2)[-1], sorted(f3)[-1]
    print("2D 文件:", f2)
    print("3D 文件:", f3)
    worst = 0.0
    for key, dspath in DS.items():
        a = read_array(f2, dspath)                 # 2D 场：idx = i + j*NX
        b = read_array(f3, dspath)                 # 3D 场：idx = i + j*NX + k*NX*NX
        n2 = NX * NX
        dmax = 0.0; amax = 0.0
        cnt = 0
        for k in range(1):                          # z 第一层
            for j in range(NX):
                for i in range(NX):
                    va = a[i + j * NX]
                    vb = b[i + j * NX + k * n2]
                    d = abs(va - vb)
                    if d > dmax: dmax = d
                    aa = abs(va)
                    if aa > amax: amax = aa
        rel = dmax / amax if amax > 0 else dmax
        print("  %-3s: max|2D-3D| = %.3e   相对 = %.3e" % (key, dmax, rel))
        worst = max(worst, rel)
        # z 均匀性：k=0 vs k=1 层
        dz = 0.0
        for j in range(NX):
            for i in range(NX):
                dz = max(dz, abs(b[i + j * NX] - b[i + j * NX + n2]))
        print("       z 均匀性（层0 vs 层1）max|diff| = %.3e" % dz)
    print("== 最大相对差 = %.3e  → %s（判据 1e-12）" %
          (worst, "PASS" if worst < 1e-12 else "FAIL"))

if __name__ == "__main__":
    main()
