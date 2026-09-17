/**
 * @file equation.hpp
 * @brief 定义方程类，用于处理守恒变量与原始变量之间的转换
 */

#pragma once

#include"data.hpp"

/**
 * @brief 方程类，用于处理不同类型方程的守恒变量与原始变量之间的转换
 * 
 * 该类提供了守恒变量到原始变量的转换功能，支持多种类型的方程，
 * 包括线性对流方程、Burgers方程和欧拉方程等。
 */
class Equation
{
    public:
    /**
     * @brief 默认构造函数
     */
    Equation(){};
    
    /**
     * @brief 将守恒变量转换为原始变量
     */
    void consToPrim();
    
    /**
     * @brief 获取守恒变量数据指针
     * @return 守恒变量数据指针
     */
    Data* getCons();
    
    /**
     * @brief 获取原始变量数据指针
     * @return 原始变量数据指针
     */
    Data* getPrim();
    
    /**
     * @brief 获取右手端数据指针
     * @return 右手端数据指针
     */
    Data* getRhs();

    protected:
    friend class Initializer;

    /**
     * @brief 一维欧拉方程的守恒变量到原始变量转换
     */
    void consToPrimEuler1D();
    
    /**
     * @brief 二维欧拉方程的守恒变量到原始变量转换
     */
    void consToPrimEuler2D();

    /**
     * @brief 三维欧拉方程的守恒变量到原始变量转换
     *
     * 将守恒变量（rho, rho*u, rho*v, rho*w, rho*E）转换为原始变量（rho, u, v, w, p）
     */
    void consToPrimEuler3D();

    bool inited=false;                ///< 初始化标志
    int n;                      ///< 网格点数
    int nCons;                  ///< 守恒变量数
    int nPrim;                  ///< 原始变量数
    int dim;                    ///< 问题维度
    Data* cons;                 ///< 守恒变量数据指针
    Data* rhs;                  ///< 右手端数据指针
    Data* prim;                 ///< 原始变量数据指针
    EquationType type;          ///< 方程类型
};