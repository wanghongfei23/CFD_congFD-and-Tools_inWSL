#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
事件复现工具：统计"边界污染"事件中的超阈值偏差占比（核验记录 §7 数字的复核命令）

用法（在仓库根目录）:
  python3 tools/incident_stats.py vortex  <cgns> [阈值=0.1]    # 等熵涡：|u_num-u_exact| 超阈值占比
  python3 tools/incident_stats.py wave_old <cgns> [阈值=0.01]  # 初版熵波(波矢 (1,2,3)/sqrt(14))：|rho_num-rho_exact|
  python3 tools/incident_stats.py wave_new <cgns> [阈值=0.01]  # 现行熵波(波矢 (1,1,1))
说明：等熵涡的解析式与 tools/vortex_error.py 相同；旧波初值源码未入库（§7），
      但其失败运行的原始数据保留于 runs/wav_n33/，可用本工具对其复算占比。
"""
import sys, os, math, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vortex_error import read_array, DS_PATHS, PI, exact as vortex_exact

def wave_old(x, y, z, t):
    u = 1.0 / math.sqrt(3.0); nrm = 1.0 / math.sqrt(14.0)
    return 1.0 + 0.2 * math.sin((x + 2.0 * y + 3.0 * z - 6.0 * u * t) * nrm)

def wave_new(x, y, z, t):
    u = 1.0 / math.sqrt(3.0)
    return 1.0 + 0.2 * math.sin(x + y + z - 3.0 * u * t)

def main():
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(2)
    mode, path = sys.argv[1], sys.argv[2]
    thr = float(sys.argv[3]) if len(sys.argv) > 3 else (0.1 if mode == "vortex" else 0.01)
    m = re.search(r"t=([0-9.]+)\.cgns", path)
    t = float(m.group(1)) if m else 0.0

    if mode == "vortex":
        num = read_array(path, DS_PATHS["u"])
        def ex(x, y, z): return vortex_exact(x, y, z, t)[1]     # u 分量
    elif mode == "wave_old":
        num = read_array(path, DS_PATHS["rho"])
        def ex(x, y, z): return wave_old(x, y, z, t)
    elif mode == "wave_new":
        num = read_array(path, DS_PATHS["rho"])
        def ex(x, y, z): return wave_new(x, y, z, t)
    else:
        print("mode 须为 vortex/wave_old/wave_new"); sys.exit(2)

    N = round(len(num) ** (1.0 / 3.0))
    assert N ** 3 == len(num), "点数非立方: %d" % len(num)
    h = 2.0 * PI / N
    linf = 0.0; cnt = 0
    for k in range(N):
        for j in range(N):
            for i in range(N):
                d = abs(num[i + j * N + k * N * N] - ex((i+.5)*h, (j+.5)*h, (k+.5)*h))
                if d > linf: linf = d
                if d > thr: cnt += 1
    print("%s  t=%.4g  N=%d格  阈值=%.3g  →  Linf=%.4e  超阈值格点占比=%.2f%%"
          % (mode, t, N, thr, linf, 100.0 * cnt / N ** 3))

if __name__ == "__main__":
    main()
