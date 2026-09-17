/**
 * @file fixedDtWaveTest.cpp
 * @brief 固定 dt 复核程序（独立审查建议）：以固定小步长重算熵波三档，
 *        排除"CFL 固定（dt∝h）导致时间误差掺入"的可能，验证空间收敛阶。
 *
 * 用法: ./fixedDtWaveTest <iMax> [步数=4000] [dt=5e-4]
 * 输出: 当前目录下 3D_wave 的 CGNS（t = 步数×dt），由 tools/wave_error.py 对比。
 */
#include "blockSolver.hpp"
#include <iostream>

int main(int argc, char** argv)
{
    const int  iMax  = (argc > 1) ? std::atoi(argv[1]) : 33;   // 网格点数（格数 = iMax-1）
    const int  nStep = (argc > 2) ? std::atoi(argv[2]) : 4000; // 默认 4000 步
    const real dt    = (argc > 3) ? std::atof(argv[3]) : 5e-4; // 默认固定 dt=5e-4 → t=2
    const real pi    = 3.14159265358979323846;

    Info* i3 = new Info;
    i3->dim    = 3;
    i3->nCase  = 3;                                  // 熵波（与 runs/wav2_* 同初值）
    i3->calZone = { 0, 2*pi, 0, 2*pi, 0, 2*pi };
    i3->iMax   = { iMax, iMax, iMax };
    i3->interMethod = (InterMethod)31;
    BlockSolver b3(i3);

    for (int s = 0; s < nStep; s++) b3.solve(dt);

    std::cout << "fixedDtWave: iMax=" << iMax << " nStep=" << nStep
              << " dt=" << dt << "  t=" << i3->t << "\n";
    b3.outputGrid();
    b3.outputPrim();
    std::cout << "done\n";
    return 0;
}
