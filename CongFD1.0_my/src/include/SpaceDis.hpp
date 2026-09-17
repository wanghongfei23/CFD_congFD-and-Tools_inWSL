/**
 * @file SpaceDis.hpp
 * @brief 空间离散类定义文件
 * @author [作者姓名]
 * @date [创建日期]
 */

#pragma once
#include "block.hpp"
#include "oneDBnd.hpp"
#include "info.hpp"
#include <boost/circular_buffer.hpp>
#include "eigenSystem.hpp"

/**
 * @brief 空间离散类，用于处理计算流体力学中的空间离散化
 * 
 * 该类负责执行空间离散化操作，包括通量计算、重构方法以及不同的离散格式。
 * 支持多种方程类型和差分方法。
 */
class SpaceDis
{
    public:
    /**
     * @brief 构造函数，初始化空间离散对象
     * @param n_ 网格点数量
     * @param data_ 数据指针
     * @param rhs_ 右-hand side数据指针
     * @param bndL_ 左边界条件
     * @param bndR_ 右边界条件
     * @param info_ 信息对象指针
     */
    SpaceDis(int n_,Data* data_,Data* rhs_
            ,std::shared_ptr<OneDBnd> bndL_,std::shared_ptr<OneDBnd> bndR_,Info* info);
    
    /**
     * @brief 默认构造函数
     */
    SpaceDis();
    
    /**
     * @brief 执行差分计算
     */
    void difference();
    
    /**
     * @brief 设置偏移量
     * @param i0 起始索引
     * @param offset 偏移值
     */
    void setOffset(int,int);
    
    /**
     * @brief 设置方程类型和差分方法
     * @param eqType 方程类型
     * @param diffMethod 差分方法
     */
    void setMethod(EquationType ,DiffMethod);
    
    /**
     * @brief 设置空间维度
     * @param idim 空间维度索引
     */
    void setIDim(int);

    /**
     * @brief 设置法向量
     * @param norm 法向量数组
     */
    void setConstNorm(std::array<real,3>&&);

    /// 计时相关变量
    long long timep=0;
    
    private:
    
    Data* data; ///< 数据指针
    Data* rhs; ///< 右-hand side数据指针
    Info* info; ///< 信息对象指针
    std::shared_ptr<OneDBnd> bndL; ///< 左边界条件
    std::shared_ptr<OneDBnd> bndR; ///< 右边界条件

    /**
     * @brief 计算通量
     */
    void calFlux();
    
    /**
     * @brief 计算对流通量
     * @param i 网格点索引
     */
    void calFluxConv(int);
    
    /**
     * @brief 计算Burgers方程通量
     * @param i 网格点索引
     */
    void calFluxBurgers(int);
    
    /**
     * @brief 计算一维欧拉方程通量
     * @param i 网格点索引
     */
    void calFluxEuler1D(int);
    
    /**
     * @brief 计算二维欧拉方程通量
     * @param i 网格点索引
     */
    void calFluxEuler2D(int);

    /**
     * @brief 计算三维欧拉方程通量（任务 5 新增）
     * @param i 网格点索引
     */
    void calFluxEuler3D(int);

    /**
     * @brief 计算精度测试通量
     * @param i 网格点索引
     */
    void calFluxAccuracyTest(int i);
    
    /// 函数指针，指向特定类型的通量计算函数
    void (SpaceDis::*calTypeFlux)(int);

    /**
     * @brief 访问数据元素
     * @param i 第一个索引
     * @param j 第二个索引
     * @return 对应位置的数据引用
     */
    real& at(int,int);
    
    /**
     * @brief 访问通量数据元素
     * @param i 第一个索引
     * @param j 第二个索引
     * @return 对应位置的通量数据引用
     */
    real& fluxAt(int,int);
    
    /// 函数指针，指向左界面重构方法
    std::vector<real> (SpaceDis::*reconLMethod)(int i);
    
    /// 函数指针，指向右界面重构方法
    std::vector<real> (SpaceDis::*reconRMethod)(int i);
    
    /**
     * @brief 左界面重构方法
     * @param i 网格点索引
     * @return 重构后的变量值向量
     */
    std::vector<real> reconL(int);
    
    /**
     * @brief 左界面原始变量重构方法
     * @param i 网格点索引
     * @return 重构后的原始变量值向量
     */
    std::vector<real> reconLprim(int);
    
    /**
     * @brief 一维左界面特征重构方法
     * @param i 网格点索引
     * @return 重构后的特征变量值向量
     */
    std::vector<real> reconLChar1D(int i);
    
    /**
     * @brief 二维左界面特征重构方法
     * @param i 网格点索引
     * @return 重构后的特征变量值向量
     */
    std::vector<real> reconLChar2D(int i);
    
    /**
     * @brief 右界面重构方法
     * @param i 网格点索引
     * @return 重构后的变量值向量
     */
    std::vector<real> reconR(int);
    
    /**
     * @brief 右界面原始变量重构方法
     * @param i 网格点索引
     * @return 重构后的原始变量值向量
     */
    std::vector<real> reconRprim(int);
    
    /**
     * @brief 一维右界面特征重构方法
     * @param i 网格点索引
     * @return 重构后的特征变量值向量
     */
    std::vector<real> reconRChar1D(int i);
    
    /**
     * @brief 二维右界面特征重构方法
     * @param i 网格点索引
     * @return 重构后的特征变量值向量
     */
    std::vector<real> reconRChar2D(int i);

    /**
     * @brief 一维界面中心重构方法
     * @param i 网格点索引
     * @return 重构后的变量值向量
     */
    std::vector<real> recon1DFaceCenter(int i);
    
    /**
     * @brief 二维界面中心重构方法
     * @param i 网格点索引
     * @return 重构后的变量值向量
     */
    std::vector<real> recon2DFaceCenter(int i);

    /**
     * @brief 三维界面中心重构方法（任务 5 新增）
     * @param i 网格点索引
     * @return 重构后的变量值向量（左右各 5 个，共 10 个）
     */
    std::vector<real> recon3DFaceCenter(int i);

    /// 五点插值函数指针
    real (*inter5) (std::array<real,5>)=nullptr;
    
    /// 正五点插值函数指针
    real (*inter5Positive) (std::array<real,5>)=nullptr;

    /**
     * @brief HCS差分方法
     */
    void difHCS();
    
    /**
     * @brief 传统六阶差分方法
     */
    void difTraditional6();
    
    /**
     * @brief 二阶差分方法
     */
    void dif2Order();
    
    /**
     * @brief MND6差分方法
     */
    void difMND6();
    
    /// 函数指针，指向具体的差分方法实现
    void (SpaceDis::*difMethod)();
    
    int n;      ///< 网格点数量
    int nVar;   ///< 变量数量
    int nPrim;  ///< 原始变量数量
    int idim;   ///< 空间维度索引
    int i0=0;     ///< 起始索引，默认为0
    int offset=1; ///< 偏移值，默认为1
    int center=0; ///< 中心索引，默认为0
    int centerOffset=1; ///< 中心偏移值，默认为1
    int nHalf;  ///< 半网格点数

    std::array<real,3> norm; ///< 法向量

    std::shared_ptr<Data> flux_d; ///< 通量数据
    std::shared_ptr<OneDBnd> fBndL; ///< 左边界通量
    std::shared_ptr<OneDBnd> fBndR; ///< 右边界通量
    EquationType fluxType;    ///< 通量类型
    DiffMethod diffMethod=TRAD2;    ///< 差分方法，默认为TRAD2
    InterMethod interMethod=FIRSTORDER;  ///< 插值方法，默认为FIRSTORDER

};