/**
 * @file dim3CrossTest.cpp
 * @brief 任务 8.6：2D↔3D code-to-code 对照（固定 dt 同步推进、逐点比对）
 *
 * 同一初值（2D TGV）分别在 dim=2（nCase=8）与 dim=3（nCase=4，z 均匀嵌入）下，
 * 用【非 CFL、固定 dt】的 BlockSolver::solve() 推进相同步数——两版共用 RK3 时间
 * 积分与逐维框架，唯一差异是空间离散的 2D/3D 分派链（面心重构、Roe 通量、
 * 差分、周期装配），因此两条数值轨迹应在舍入级一致。
 * 输出：各自 outputPrim() 的 CGNS（2D 全场 / 3D 各 z 层），比对脚本逐点核验。
 *
 * 用法: ./dim3CrossTest [步数（默认1000，dt=0.001 → t=1）]
 */
#include "blockSolver.hpp"
#include <iostream>

int main(int argc, char** argv)
{
    const real dt    = 0.001;                                        // 固定步长（远小于显式稳定限）
    const int  nStep = (argc > 1) ? std::atoi(argv[1]) : 1000;       // 默认 1000 步 → t=1.0
    const real pi    = 3.14159265358979323846;

    // ---- 2D：dim=2、nCase=8（2D TGV，四面周期）；末维惯例 2（即 1 格）----
    Info* i2 = new Info;
    i2->dim    = 2;
    i2->nCase  = 8;
    i2->calZone = { 0, 2*pi, 0, 2*pi, 0, 0 };
    i2->iMax   = { 33, 33, 2 };
    i2->interMethod = (InterMethod)31;                               // 与主算例一致：TENO-AS-Ff_5_10
    BlockSolver b2(i2);

    // ---- 3D：dim=3、nCase=4（2D TGV 的 z 均匀嵌入，六面周期；z 方向 8 格）----
    // 注：单向格数须 ≥ 幽灵层数（当前配置 nGhost=4），否则周期回卷需绕域多圈、
    // 超出 OneDBnd 单次映射能力（实测 z=2 格时出现幽灵层越界→局部负压；z=8 格正常）。
    Info* i3 = new Info;
    i3->dim    = 3;
    i3->nCase  = 4;
    i3->calZone = { 0, 2*pi, 0, 2*pi, 0, 2*pi };
    i3->iMax   = { 33, 33, 9 };
    i3->interMethod = (InterMethod)31;
    BlockSolver b3(i3);

    // ---- 固定 dt 同步推进相同步数 ----
    for (int s = 0; s < nStep; s++)
    {
        b2.solve(dt);
        b3.solve(dt);
    }
    std::cout << "cross test: nStep=" << nStep << " dt=" << dt
              << "  t(2D)=" << i2->t << "  t(3D)=" << i3->t << "\n";

    // ---- 输出（先写网格 Base/Zone，再追加解；与主程序 output 流程一致）----
    b2.outputGrid();
    b2.outputPrim();
    b3.outputGrid();
    b3.outputPrim();
    std::cout << "cross test done\n";
    return 0;
}
