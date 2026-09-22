/**
 * @file sp_distributor.hpp
 * @brief 定义空间分布器类，用于分发空间离散任务
 */

#pragma once
#include "Bnds.hpp"
#include "SpaceDis.hpp"
#include <omp.h>

/**
 * @brief 空间分布器类
 * 
 * 该类负责将空间离散任务分发到各个计算维度，支持1D/2D/3D问题的并行计算。
 */
class SpDistributor
{
    public:
    /**
     * @brief 求解右手端
     * 
     * 对各个维度分别进行空间离散计算，得到右端项。
     */
    void rhsSolve();

    private:
    friend class Initializer;

    int nCons;              ///< 守恒变量数
    int nPrim;              ///< 原始变量数
    int dim;                ///< 问题维度
    std::array<int,3> iMax; ///< 各方向最大索引
    Data* prim;             ///< 原始变量数据指针
    Data* cons;             ///< 守恒变量数据指针
    Data* rhs;              ///< 右手端数据指针
    Bnds* bnds;             ///< 边界条件管理器指针
    Info* info;             ///< Info对象指针
};