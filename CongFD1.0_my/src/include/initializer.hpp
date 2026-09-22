/**
 * @file initializer.hpp
 * @brief 定义初始化器类，用于初始化网格、方程和边界条件
 */

#pragma once
#include"block.hpp"
#include"info.hpp"
#include"equation.hpp"
#include"Bnds.hpp"
#include"sp_distributor.hpp"

/**
 * @brief 初始化器类
 * 
 * 该类负责初始化网格、方程、边界条件和空间分布器等组件。
 */
class Initializer
{
    public:
    /**
     * @brief 默认构造函数
     */
    Initializer();
    
    /**
     * @brief 带参数的构造函数
     * @param 第二个参数 Info对象指针
     */
    Initializer(Info*);
    
    /**
     * @brief 解初始化
     * @param 第一个参数 网格块指针
     * @param 第二个参数 解数据指针
     */
    void solInit(Block*,Data*);
    
    /**
     * @brief 初始化均匀网格
     * @param 第一个参数 网格块指针
     */
    void initUniformBlock(Block*);
    
    /**
     * @brief 初始化方程
     * @param 第一个参数 方程对象指针
     * @param 第二个参数 网格块指针
     */
    void initEqution(Equation*,Block*);

    /**
     * @brief 初始化边界条件
     * @param bnds 边界条件管理器指针
     * @param 第二个参数 方程对象指针
     * @param iMax 网格最大索引数组
     * @param block 网格块指针
     */
    void initBnds(Bnds* bnds,Equation*,std::array<int,3> iMax,Block* block);
    
    /**
     * @brief 初始化双马赫反射问题边界条件
     * @param bnds 边界条件管理器指针
     * @param eqn 方程对象指针
     * @param iMax 网格最大索引数组
     * @param block 网格块指针
     */
    void initDoubleMachBnds(Bnds* bnds,Equation* eqn,std::array<int,3> iMax,Block* block);

    /**
     * @brief 初始化激波/剪切层相互作用问题边界条件
     * @param bnds 边界条件管理器指针
     * @param eqn 方程对象指针
     * @param iMax 网格最大索引数组
     * @param block 网格块指针
     */
    void initShockShearBnds(Bnds* bnds,Equation* eqn,std::array<int,3> iMax,Block* block);

    /**
     * @brief 初始化空间分布器
     * @param 第一个参数 空间分布器指针
     * @param 第二个参数 方程对象指针
     * @param 第三个参数 网格块指针
     * @param 第四个参数 边界条件管理器指针
     */
    void initSpDistributor(SpDistributor*,Equation*,Block*,Bnds*);

    private:
    Info* info;     ///< Info对象指针
};