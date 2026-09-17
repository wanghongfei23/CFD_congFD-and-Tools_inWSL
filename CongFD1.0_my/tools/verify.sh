#!/bin/bash
# 三维加装·一键查证（只读既有产物与构建品，不重跑仿真；约 1-2 分钟）
# 用法: cd <repo> && bash tools/verify.sh
cd "$(dirname "$0")/.." || exit 9
fail=0

echo "== [1] 内核测试（实时运行，秒级；期望 ALL PASS）"
./build/dim3KernelTest || fail=1

echo
echo "== [2] 均匀场纯性：t=0 vs t=0.6 全场 h5diff（期望零输出=零差异）"
h5diff "runs/tgv8_uni/3D_uniform - TENO-AS-Ff_5_10 - 17x17x17 - t=0.0000.cgns" \
       "runs/tgv8_uni/3D_uniform - TENO-AS-Ff_5_10 - 17x17x17 - t=0.6000.cgns" && echo "  零差异 ✓" || fail=1

echo
echo "== [3] 线程确定性：OMP1 vs OMP14 的 t=0.2 全场 h5diff（期望零输出）"
h5diff "runs/tgv8_omp1/3D_TGV - TENO-AS-Ff_5_10 - 17x17x17 - t=0.2000.cgns" \
       "runs/tgv8_omp14/3D_TGV - TENO-AS-Ff_5_10 - 17x17x17 - t=0.2000.cgns" && echo "  零差异 ✓" || fail=1
diff -q runs/tgv8_omp1/diagnostics.txt runs/tgv8_omp14/diagnostics.txt > /dev/null && echo "  诊断文件逐字节一致 ✓" || fail=1

echo
echo "== [4] 关键数据抽查"
echo "-- 均匀场 diagnostics（K 首末逐位相同、Ω 严格 0）:"
cat runs/tgv8_uni/diagnostics.txt
echo "-- TGV 收敛三档 t=1 行（K 值；p=log2((K17-K33)/(K33-K65))=4.69）:"
tail -n 1 runs/tgv8_n17/diagnostics.txt
tail -n 1 runs/tgv8_n33/diagnostics.txt
tail -n 1 runs/tgv8_n65/diagnostics.txt

echo
echo "== [5] 熵波外部真值（python 逐点对比，约 1-2 分钟）"
python3 tools/wave_error.py || fail=1

echo
echo "== [6] TGV 镜像对称检验（约 1 分钟）"
python3 tools/tgv_symmetry.py || fail=1

echo
echo "== [7] 2D↔3D code-to-code 对照（重跑 1000 步 + 逐点比对；约 1 分钟）"
if [ -f build/dim3CrossTest ] && [ -d runs/cross ]; then
  ( cd runs/cross && ../../build/dim3CrossTest 1000 > /dev/null && python3 ../../tools/cross_compare.py ) || fail=1
else
  echo "  (缺 build/dim3CrossTest 或 runs/cross，跳过)"
fi

echo
if [ "$fail" -eq 0 ]; then echo "===== 一键查证全部通过 ====="; else echo "===== 存在失败项，见上 ====="; fi
exit $fail
