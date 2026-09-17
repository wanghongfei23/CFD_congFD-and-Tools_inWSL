#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TGV 镜像对称检验（任务 8 系统对称性证据的查证脚本）
读取 runs/tgv8_n65 的 TGV 场（64^3, t=1），检验四组精确镜像关系的破缺量：
  x-反射: u 反对称、v 对称；y-反射: u 对称、v 反对称
（TGV 初值严格满足；格式若保持对称性，破缺应为机器零量级）

用法:
  cd <repo> && python3 tools/tgv_symmetry.py [cgns文件（默认 runs/tgv8_n65 的 t=1 档）]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vortex_error import read_array, DS_PATHS

DEFAULT = "runs/tgv8_n65/3D_TGV - TENO-AS-Ff_5_10 - 65x65x65 - t=1.0000.cgns"

def main():
    f = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    if not os.path.exists(f):
        print("文件不存在: %s" % f); sys.exit(2)
    N = 64
    u = read_array(f, DS_PATHS["u"])
    v = read_array(f, DS_PATHS["v"])
    idx = lambda i, j, k: i + j * N + k * N * N
    eu = ev = eu2 = ev2 = 0.0
    for k in range(N):
        for j in range(N):
            for i in range(N):
                ii = N - 1 - i; jj = N - 1 - j
                d1 = abs(u[idx(i, j, k)] + u[idx(ii, j, k)])
                d2 = abs(v[idx(i, j, k)] - v[idx(ii, j, k)])
                d3 = abs(u[idx(i, j, k)] - u[idx(i, jj, k)])
                d4 = abs(v[idx(i, j, k)] + v[idx(i, jj, k)])
                if d1 > eu: eu = d1
                if d2 > ev: ev = d2
                if d3 > eu2: eu2 = d3
                if d4 > ev2: ev2 = d4
    print("TGV 镜像对称破缺（%s）:" % f)
    print("  x-反射: u 反对称 %.3e   v 对称 %.3e" % (eu, ev))
    print("  y-反射: u 对称 %.3e   v 反对称 %.3e" % (eu2, ev2))
    ok = max(eu, ev, eu2, ev2) < 1e-12
    print("  判定（阈值 1e-12）: %s" % ("PASS" if ok else "FAIL"))
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
