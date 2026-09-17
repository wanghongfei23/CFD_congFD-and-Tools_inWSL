/**
 * @file diagnostics.hpp
 * @brief 三维诊断输出：动能 K 与涡度（enstrophy）Ω 的体积分（任务 7 新增）
 */

#pragma once
#include "data.hpp"
#include "macro.hpp"

/**
 * @brief 计算三维动能 K 与涡度 Ω（周期域体积分）
 * @param prim 原始变量场（列序：rho, u, v, w, p）
 * @param icMax 各方向内部格数（格心网格）
 * @param h 均匀网格间距（三向相同）
 * @param K 输出：K = Σ ρ(u²+v²+w²)/2 · h³
 * @param Omega 输出：Ω = Σ (ωx²+ωy²+ωz²)/2 · h³（四阶中心差分 + 周期回卷）
 */
void calcDiagnostics(Data* prim, const std::array<int,3>& icMax, real h, real& K, real& Omega);
