/**
 * @file Bnds.hpp
 * @brief 定义边界条件管理类
 */

#pragma once
#include "info.hpp"
#include "oneDBnd.hpp"
#include "data.hpp"

/**
 * @brief 边界条件管理类
 * 
 * 该类用于管理计算域的边界条件，包括获取一维边界和更新边界数据。
 */
class Bnds
{
    public:
    /**
     * @brief 获取指定方向的一维边界
     * @param idim 方向索引
     * @param i 第一索引
     * @param j 第二索引
     * @return 包含左右边界的数组
     */
    std::array<std::shared_ptr<OneDBnd>,2> getOneDBnd(int,int,int);
    
    /**
     * @brief 更新所有边界条件
     */
    void update();

    private:
    friend class Initializer;
    int dim;                                    ///< 问题维度
    std::array<int,3> iMax;                     ///< 各方向最大索引
    //idim*
    std::vector<std::shared_ptr<OneDBnd>> oneDBnds; ///< 一维边界条件向量
};