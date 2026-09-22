/**
 * @file info.hpp
 * @brief 定义求解器相关信息的类
 */

#pragma once
#include "macro.hpp"

/**
 * @class Info
 * @brief 存储和管理求解器的各种配置参数和状态信息
 * 
 * Info 类包含了求解器运行所需的所有配置参数，包括方程类型、空间离散方法、
 * 时间步长控制、网格信息、边界条件等。同时也跟踪求解过程的状态信息，
 * 如当前时间、步数等。
 */
class Info
{
    public:
    InterMethod spMethod=WCNS5;                    ///< 空间离散方法
    EquationType eqType=EULER;                      ///< 控制方程类型
    int nCase=2;                                   ///< 算例编号
    real t=0;                                      ///< 当前物理时间
    real dt=0.0001;                                ///< 时间步长
    int step=0;                                    ///< 当前步数
    int endStep=25000;                             ///< 结束步数
    int dim=1;                                     ///< 问题维度

    int outputInterval=100;                        ///< 输出间隔步数
    real outputDt=0.01;                            ///< 输出时间间隔
    real outputT=0;                                ///< 上次输出时间

    real CFL=0.5;                                  ///< CFL数
    bool fixedtimeSteps=true;                      ///< 是否使用固定时间步长

    //粘性扩展（任务 4；info.txt 第 17 至 19 项，缺省粘性关，旧配置零改动）
    bool viscous=false;                            ///< 粘性开关（常粘 μ=1/Re）
    real Re=1600;                                  ///< 雷诺数（本档 1600 → μ=6.25e-4）
    real Pr=0.72;                                  ///< 普朗特数（κ=μc_p/Pr；c_p=γ/(γ−1)）

    //幂律粘性扩展（info.txt 第 20 至 21 项，缺省常粘，旧配置零改动；HIT 开发）
    real powerLawExp=0.0;                          ///< 幂律指数（μ=μref·(T/Tref)^n；≤0=常粘）
    real Tref=1.0;                                 ///< 幂律参考温度 T_ref（幂律开时有效）

    //对于隐式求解器
    real implicitCFL=0.01;                         ///< 隐式求解器CFL数
    real maxImplicitStep=100;                      ///< 最大隐式步数

    DiffMethod diffMethod=TRAD6;                   ///< 差分方法
    InterMethod interMethod=WCNS5;                 ///< 插值方法
    // 王鸿飞 默认设置处  RIT重力源项
    SourceType sourceType=SOURCENULL;              ///< 源项类型
    // SourceType sourceType=GRAVITY;                 ///< RTI的重力源项

    std::array<int,3> iMax{201,201,2};             ///< 网格在各方向上的节点数
    std::array<double,6> calZone{0,0.3,0,0.3,0,2}; ///< 计算域范围[xmin,xmax,ymin,ymax,zmin,zmax]

    /**
     * @brief 构造函数
     */
    Info();
    
    /**
     * @brief 获取虚拟点单元数
     * @return 虚拟点单元数量
     */
    int nGhostCell();
    
    /**
     * @brief 获取通量点数
     * @return 通量点数量
     */
    int nFluxPoint();
    
    /**
     * @brief 获取原始变量数
     * @return 原始变量数量
     */
    int nPrim();
    
    /**
     * @brief 获取守恒变量数
     * @return 守恒变量数量
     */
    int nCons();
    
    /**
     * @brief 获取问题维度
     * @return 问题维度
     */
    int getDim();

    bool constH=true;                              ///< 是否使用恒定网格间距
    real interval=0;                               ///< 网格间距
    /**
     * @brief 获取指定方向的网格间距
     * @param idim 方向索引(0-x方向, 1-y方向, 2-z方向)
     * @return 网格间距
     */
    real geth(int);

    /**
     * @brief 获取默认边界类型
     * @return 默认边界类型
     */
    BndType defaultBndType();
    
    /**
     * @brief 获取内部网格最大索引
     * @return 内部网格最大索引数组
     */
    std::array<int,3> icMax();
    
    /**
     * @brief 生成文件名
     * @return 生成的文件名字符串
     */
    std::string filename();

    /**
     * @brief 获取守恒变量名称列表
     * @return 守恒变量名称向量
     */
    std::vector<std::string> getVarNameListCons();
    
    /**
     * @brief 获取原始变量名称列表
     * @return 原始变量名称向量
     */
    std::vector<std::string> getVarNameListPrim();
    
    /**
     * @brief 获取RHS变量名称列表
     * @return RHS变量名称向量
     */
    std::vector<std::string> getVarNameListRhs();
};