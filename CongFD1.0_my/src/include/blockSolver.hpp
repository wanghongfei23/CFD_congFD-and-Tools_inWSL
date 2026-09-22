#pragma once
#include "initializer.hpp"
#include "cgnsio.hpp"
#include "SourceTerm.hpp"
#include "ViscousTerm.hpp"

/**
 * @brief 区块求解器类
 * 
 * 该类负责管理CFD计算的主要流程，包括初始化、时间推进、输出等核心功能。
 * 它集成了网格读写、方程求解、边界条件处理等功能模块。
 */
class BlockSolver
{
    public:
    /**
     * @brief 默认构造函数
     */
    BlockSolver();
    
    /**
     * @brief 带参数的构造函数
     * @param info 指向信息对象的指针，包含求解所需的各种参数
     */
    BlockSolver(Info*);
    
    /**
     * @brief 执行求解过程
     * @param t 时间步长或相关参数
     */
    void solve(real);
    
    /**
     * @brief 析构函数
     */
    ~BlockSolver();
    
    /// 时间步数计数器
    long timesteps=0;

    /**
     * @brief 输出网格信息
     */
    void outputGrid();
    
    /**
     * @brief 输出守恒变量
     */
    void outputCons();
    
    /**
     * @brief 输出原始变量
     */
    void outputPrim();

    /**
     * @brief 时间步循环
     */
    void stepsLoop();
    
    /**
     * @brief 基于CFL条件的时间步循环
     */
    void stepsLoopCFL();
    
    /**
     * @brief 基于固定时间步长的时间步循环
     */
    void stepsLoopDTS();
    
    /**
     * @brief 测试函数
     */
    void Test();

    /**
     * @brief 取守恒变量数据指针（测试驱动用；只读访问，无行为影响）
     */
    Data* getConsPtr(){ return cons; }

    /**
     * @brief 取原始变量数据指针（测试驱动用；只读访问，无行为影响）
     */
    Data* getPrimPtr(){ return eqn->getPrim(); }

    private:
    
    CgnsIO cgnsIO;              ///< CGNS输入输出处理对象
    Info* info;                 ///< 信息对象指针
    Block* block;               ///< 区块对象指针
    Initializer* initer;        ///< 初始化器对象指针
    Equation* eqn;              ///< 方程对象指针
    Bnds* bnds;                 ///< 边界条件对象指针
    SpDistributor* spDis;       ///< 空间分布器对象指针
    Data* cons,*rhs;            ///< 守恒变量和右端项数据指针
    SourceTerm* sourceTerm;     ///< 源项对象指针
    ViscousTerm* viscousTerm;   ///< 粘性项对象指针（任务 4）
    
    /**
     * @brief 三阶SSP Runge-Kutta时间积分方法
     * @param dt 时间步长
     */
    void RK3_SSP(real);
    
    /**
     * @brief 四阶SSP Runge-Kutta时间积分方法
     * @param dt 时间步长
     */
    void RK4_SSP(real);
    
    /**
     * @brief 显式欧拉时间积分方法
     * @param dt 时间步长
     */
    void DTS_Euler(real);
    
    /**
     * @brief 计算局部CFL数
     * @return 包含各点CFL数的向量
     */
    std::vector<real> calLocalCFL();
    
    /**
     * @brief 获取显式时间积分的时间间隔
     * @return 时间间隔值
     */
    real getTimeIntervalExplicit();

    /**
     * @brief 获取粘性档显式时间间隔（对流＋黏性谱半径；不动无粘路线，开发方案 §4.3）
     * @return 时间间隔值
     */
    real getTimeIntervalExplicitViscous();
    
    TimeMethod timeMethod=RK3SSP;  ///< 时间积分方法，默认为三阶SSP Runge-Kutta方法
};