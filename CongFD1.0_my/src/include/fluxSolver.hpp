/**
 * @file fluxSolver.hpp
 * @brief 定义通量求解器基类，用于计算数值通量
 */

#pragma once

#include "data.hpp"
#include "info.hpp"

/**
 * @brief 通量求解器基类
 * 
 * 该类是所有通量求解器的基类，定义了通量计算的基本接口。
 */
class fluxSolver
{
    public:
    /**
     * @brief 构造函数
     * @param valsR_ 重构后数据指针
     * @param fluxes_ 通量数据指针
     * @param info_ Info对象指针
     * @param idim_ 维度索引
     */
    fluxSolver(std::shared_ptr<Data> valsR_,std::shared_ptr<Data> fluxes,Info* info_,int idim_);
    
    /**
     * @brief 通量求解虚函数，需要在派生类中实现
     */
    virtual void fluxSolve()=0;

    protected:
    std::shared_ptr<Data> valsR;    ///< 重构后数据指针
    std::shared_ptr<Data> fluxes;   ///< 通量数据指针
    Info* info;                     ///< Info对象指针
    int idim;                       ///< 维度索引
};