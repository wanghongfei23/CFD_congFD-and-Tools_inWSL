/**
 * @file macro.hpp
 * @brief 定义项目中使用的宏、枚举类型和通用函数
 */

#pragma once

#define real double              ///< 实数类型别名
#define ind int                 ///< 整数类型别名

// #include <cmath>
// #include <stdio.h>
#include <array>
#include <cstring>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <format>
#include <map>

#include <chrono>
inline long timepp = 0;         ///< 时间统计变量
inline long timesss = 0;        ///< 时间统计变量

// 【王鸿飞】注意gamma
#define GAMMA 1.4               ///< 比热比常数（其他算例）
// #define GAMMA 5.0/3.0        ///< 比热比常数（RTI算例）

/**
 * @brief 边界类型枚举
 */
enum BndType {
    TYPENULL,           ///< 空边界类型
    PERIODIC1D,         ///< 一维周期边界
    SYMMETRY1D,         ///< 一维对称边界
    DIRICLET,           ///< Dirichlet边界
    DIRICLET_SODL,      ///< Sod问题左边界
    DIRICLET_SODR,      ///< Sod问题右边界
    FLUXGHOST,          ///< 通量Ghost边界
    SUPERSONICOUTLET,   ///< 超声速出口边界
    SYMMETRYX,          ///< X方向对称边界（仅用于2D）
    SYMMETRYY,          ///< Y方向对称边界（仅用于2D）
    DoubleMachUp,       ///< 双马赫反射上边界
    ARDBND,             ///< ARD边界
    ShockShearIn        ///< 激波/剪切层算例左入流（时变扰动）
};

/**
 * @brief 插值方法枚举
 */
enum InterMethod {
    // 0
    FIRSTORDER,      // 0
    // 1-10
    MUSCL,           // 1
    WCNS5,           // 2
    WCNSZ5,          // 3
    WCNS5Char,       // 4
    WCNSZ5Char,      // 5
    WCNS5CONG,       // 6
    TCNSCongA,       // 7
    WCNS5CONGZ,      // 8
    WCNS5CONGZCT4,   // 9
    WCNS5CONGZCT7,   // 10
    // 11-20
    TCNS5,           // 11
    TCNS5CT4,        // 12
    TCNS5CT7,        // 13
    LINEAR5,         // 14
    NICEST5,         // 15
    INTERMAX,        // 16
    WHFTCNSA,        // 17
    // 21-30
    // 31-40
    WHFTCNSASFf3_5_9_time_improve,     // 35
    WHFTCNSLADSFf2_5_10,     // 37
    WHFTCNSLAD,     // 38
    WHFTCNSLADS_g1,     // 39
};

/**
 * @brief 差分方法枚举
 */
enum DiffMethod {
    HDS6,                       ///< HDS6差分方法
    TRAD6,                      ///< TRAD6差分方法
    TRAD2,                      ///< TRAD2差分方法
    MND6                        ///< MND6差分方法
};

/**
 * @brief 方程类型枚举
 */
enum EquationType {
    LINEARCONV1D,               ///< 一维线性对流方程
    BURGERS1D,                  ///< 一维Burgers方程
    EULER,                      ///< 欧拉方程
    ACCURACYTEST                ///< 精度测试方程
};

/**
 * @brief 通量方法枚举
 */
enum FluxMethod {
    HLLC1D,                     ///< 一维HLLC通量
    ROE1D,                      ///< 一维Roe通量
    HLLC2D                      ///< 二维HLLC通量
};

/**
 * @brief 时间积分方法枚举
 */
enum TimeMethod {
    RK3SSP,                     ///< 三阶SSP Runge-Kutta方法
    EulerFront                  ///< 欧拉向前积分方法
};

/**
 * @brief 源项类型枚举
 */
enum SourceType {
    SOURCENULL,                 ///< 无源项
    GRAVITY,                    ///< 重力源项
    VISCOUSMMS                  ///< 粘性 MMS 造解源项（任务 4；测试用，见 04 实施分解 S9）
};

/**
 * @brief 计算三维数组的线性索引
 * @param i x方向索引
 * @param j y方向索引
 * @param k z方向索引
 * @param iMax 各方向最大索引数组
 * @return 线性索引
 */
constexpr int index(int i, int j, int k, std::array<int, 3> iMax)
{
    return i + j * iMax[0] + k * iMax[0] * iMax[1];
}

/**
 * @brief 计算偏移量
 * @param idim 维度
 * @param i 第一索引
 * @param j 第二索引
 * @param iMax 各方向最大索引数组
 * @return 包含起始索引和偏移量的数组
 */
constexpr std::array<int, 2> calOffset(int idim, int i, int j, std::array<int, 3> iMax)
{
    std::array<int, 3> offsets { 1, iMax[0], iMax[0] * iMax[1] };
    std::array<int, 2> res;
    if (idim == 1) {
        res[0] = i * offsets[1] + j * offsets[2]; // i0
        res[1] = offsets[0]; // offset
    } else if (idim == 2) {
        res[0] = i * offsets[0] + j * offsets[2]; // i0
        res[1] = offsets[1]; // offset
    } else if (idim == 3) {
        res[0] = i * offsets[0] + j * offsets[1]; // i0
        res[1] = offsets[2]; // offset
    }
    return res;
}

/**
 * @brief 计算逆向偏移量
 * @param idim 维度
 * @param i 第一索引
 * @param j 第二索引
 * @param iMax 各方向最大索引数组
 * @return 包含起始索引和偏移量的数组
 */
constexpr std::array<int, 2> calOffsetInverse(int idim, int i, int j, std::array<int, 3> iMax)
{
    std::array<int, 3> offsets { 1, iMax[0], iMax[0] * iMax[1] };
    std::array<int, 2> res;
    if (idim == 1) {
        res[0] = i * offsets[1] + j * offsets[2] + (iMax[0] - 1) * offsets[0]; // i0
        res[1] = -offsets[0]; // offset
    } else if (idim == 2) {
        res[0] = i * offsets[0] + j * offsets[2] + (iMax[1] - 1) * offsets[1]; // i0
        res[1] = -offsets[1]; // offset
    } else if (idim == 3) {
        res[0] = i * offsets[0] + j * offsets[1] + (iMax[2] - 1) * offsets[2]; // i0
        res[1] = -offsets[2]; // offset
    }
    return res;
}