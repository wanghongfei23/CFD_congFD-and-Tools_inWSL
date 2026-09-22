/**
 * @file reconstructor.hpp
 * @brief 定义重构器类，用于变量重构
 */

#pragma once
#include "macro.hpp"
#include "data.hpp"
#include "info.hpp"
#include "oneDBnd.hpp"

/**
 * @brief 重构器基类
 * 
 * 该类是所有重构器的基类，定义了重构的基本接口和通用功能。
 */
class Reconstructor
{
    public:
    /**
     * @brief 构造函数
     * @param n_ 网格点数
     * @param vals_ 原始数据指针
     * @param valsR_ 重构后数据指针
     * @param bndL_ 左边界数据指针
     * @param bndR_ 右边界数据指针
     * @param info_ Info对象指针
     * @param i0_ 起始索引
     * @param offset_ 偏移量
     */
    Reconstructor(int n_,Data* vals_,std::shared_ptr<Data> valsR_,std::shared_ptr<OneDBnd> bndL_,std::shared_ptr<OneDBnd> bndR_,Info* info_,int i0_, int offset_);
    
    /**
     * @brief 重构虚函数，需要在派生类中实现
     */
    virtual void recon()=0;
    
    std::shared_ptr<Data> valsR;    ///< 重构后数据指针

    protected:
    Data* vals;                     ///< 原始数据指针
    /**
     * @brief 访问指定位置的数据
     * @param i 索引
     * @param ivar 变量索引
     * @return 数据引用
     */
    real& at(int i,int ivar);
    int i0;                         ///< 起始索引
    int offset;                     ///< 偏移量
    std::shared_ptr<OneDBnd> bndL;  ///< 左边界数据指针
    std::shared_ptr<OneDBnd> bndR;  ///< 右边界数据指针
    Info* info;                     ///< Info对象指针
    int n;                          ///< 网格点数
    int nval;                       ///< 变量数
};

/**
 * @brief 一维欧拉方程特征重构器
 * 
 * 该类实现了一维欧拉方程的特征重构方法。
 */
class ReconEigen1DEuler : public Reconstructor
{
    public:
    /**
     * @brief 实现特征重构方法
     */
    void recon() override;
};