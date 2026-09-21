/**
 * @file diagnostics.hpp
 * @brief 三维诊断输出：K/Ω 与 HIT 扩展量 U2/Tvar/Thvar 的体积分
 */

#pragma once
#include "data.hpp"
#include "macro.hpp"

/**
 * @brief 计算三维诊断量（周期域体积分）
 * @param prim 原始变量场（列序：rho, u, v, w, p）
 * @param icMax 各方向内部格数（格心网格）
 * @param h 均匀网格间距（三向相同）
 * @param K 输出：K = Σ ρ(u²+v²+w²)/2 · h³
 * @param Omega 输出：Ω = Σ (ωx²+ωy²+ωz²)/2 · h³（四阶中心差分 + 周期回卷）
 * @param U2 输出：U2 = Σ uᵢ² · h³（未加权均方速度；Johnsen et al. 2010 图 11 口径）
 * @param Tvar 输出：Tvar = Σ (T−T̄)² · h³，T̄ 为体积平均温度（T=p/ρ）
 * @param Thvar 输出：Thvar = Σ (∇·u)² · h³（胀量平方积分）
 */
void calcDiagnostics(Data* prim, const std::array<int,3>& icMax, real h, real& K, real& Omega,
                     real& U2, real& Tvar, real& Thvar);
