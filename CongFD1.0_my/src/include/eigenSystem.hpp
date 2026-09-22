/**
 * @file eigenSystem.hpp
 * @brief 定义特征系统相关的类，用于欧拉方程的特征分解
 */

#pragma once
#include<array>
#include<macro.hpp>

enum{X,Y};

/**
 * @brief 二维欧拉方程的特征系统类
 * 
 * 该类实现了二维欧拉方程的特征分解，包括将原始变量转换为特征变量，
 * 以及将特征变量转换回原始变量的功能。
 */
class eigensystemEuler2D
{
    public:
    /**
     * @brief 默认构造函数
     */
    eigensystemEuler2D(){};
    
    /**
     * @brief 构造函数，使用原始变量和法向量初始化
     * @param prim 原始变量数组 [rho, u, v, p]
     * @param norm_ 法向量数组
     */
    eigensystemEuler2D(const std::array<real,4> & prim,const std::array<real,3> & norm_);
    
    /**
     * @brief 构造函数，使用左右两侧原始变量和法向量初始化
     * @param priml 左侧原始变量数组 [rho, u, v, p]
     * @param primr 右侧原始变量数组 [rho, u, v, p]
     * @param norm_ 法向量数组
     */
    eigensystemEuler2D(const std::array<real,4> &priml,const std::array<real,4> &primr,const std::array<real,3> & norm_);
    
    /**
     * @brief 将原始变量转换为特征变量
     * @param prim 原始变量数组 [rho, u, v, p]
     * @return 特征变量数组
     */
    std::array<real,4> primToChar(const std::array<real,4> & prim);
    
    /**
     * @brief 将特征变量转换为原始变量
     * @param chars 特征变量数组
     * @return 原始变量数组 [rho, u, v, p]
     */
    std::array<real,4> charToPrim(const std::array<real,4> & chars);
    

    private:
    real r;                     ///< 密度
    real u;                     ///< x方向速度
    real v;                     ///< y方向速度
    real p;                     ///< 压力
    real gamma;                 ///< 比热比
    real ek;                    ///< 动能
    real h;                     ///< 焓
    real c;                     ///< 声速
    real Vn;                    ///< 法向速度
    bool xOrY;                  ///< 方向标志
    std::array<real,3> norm;    ///< 法向量
    std::array<real,4*4> leftEig;   ///< 左特征矩阵
    std::array<real,4*4> rightEig;  ///< 右特征矩阵
};

/**
 * @brief 一维欧拉方程的特征系统类
 * 
 * 该类实现了一维欧拉方程的特征分解，包括将原始变量转换为特征变量，
 * 以及将特征变量转换回原始变量的功能。
 */
class eigensystemEuler1D
{
    public:
    /**
     * @brief 构造函数，使用原始变量初始化
     * @param prim 原始变量数组 [rho, u, p]
     */
    eigensystemEuler1D(const std::array<real,3> & prim);
    
    /**
     * @brief 构造函数，使用左右两侧原始变量初始化
     * @param priml 左侧原始变量数组 [rho, u, p]
     * @param primr 右侧原始变量数组 [rho, u, p]
     */
    eigensystemEuler1D(const std::array<real,3> &priml,const std::array<real,3> &primr);
    
    /**
     * @brief 将原始变量转换为特征变量
     * @param prim 原始变量数组 [rho, u, p]
     * @return 特征变量数组
     */
    std::array<real,3> primToChar(const std::array<real,3> & prim);
    
    /**
     * @brief 将特征变量转换为原始变量
     * @param chars 特征变量数组
     * @return 原始变量数组 [rho, u, p]
     */
    std::array<real,3> charToPrim(const std::array<real,3> & chars);
    

    private:
    real r;                     ///< 密度
    real u;                     ///< 速度
    real p;                     ///< 压力
    real gamma;                 ///< 比热比
    real ek;                    ///< 动能
    real h;                     ///< 焓
    real c;                     ///< 声速
    bool xOrY=false;                  ///< 方向标志
    std::array<real,3*3> leftEig;   ///< 左特征矩阵
    std::array<real,3*3> rightEig;  ///< 右特征矩阵
};

// std::array<real,4> characteristicDecomposition(std::array<real,4> prim)
// {
//     return{0,0,0,0};
// }
// std::array<real,4*4> leftEigen(const std::array<real,4>& prim,const std::array<real,3>& norm)
// {
//     std::array<real,4*4> res;
//     
//     real r=prim[0],u=prim[1],v=prim[2],p=prim[3];
//     real gamma=GAMMA,ek=(u*u+v*v)/2;
//     real h=p/r*gamma/(1-gamma);
//     real Vn=norm[0]*u+norm[1]*v;
//     real c=std::sqrt(gamma*p/r);
//     //first line
//     res[0]=-norm[Y]*u+norm[X]*v;
//     res[1]=norm[Y];
//     res[2]=-norm[X];
//     res[3]=0;
//
//     //second line
//     res[4]=h-ek;
//     res[5]=u;
//     res[6]=v;
//     res[7]=-1;
//
//     //third line
//     res[8 ]=(Vn/c+ek/h)/2;
//     res[9 ]=(-norm[X]/c-u/h)/2;
//     res[10]=(-norm[Y]/c-v/h)/2;
//     res[11]=1.0/(2*h);
//
//     //fourth line
//     res[12]=(-Vn/c+ek/h)/2;
//     res[13]=(norm[X]/c-u/h)/2;
//     res[14]=(norm[Y]/c-v/h)/2;
//     res[15]=1.0/(2*h);
//
//     return res;
// }
//
// std::array<real,4*4> rightEigen(const std::array<real,4>& prim,const std::array<real,3>& norm)
// {
//     std::array<real,4*4> res;
//     enum{X,Y};
//     real r=prim[0],u=prim[1],v=prim[2],p=prim[3];
//     real gamma=GAMMA,ek=(u*u+v*v)/2;
//     real h=p/r*gamma/(1-gamma);
//     real Vn=norm[0]*u+norm[1]*v;
//     real c=std::sqrt(gamma*p/r);
//     //first line
//     res[0]=0;
//     res[1]=1/h;
//     res[2]=1;
//     res[3]=1;
//
//     //second line
//     res[4]=norm[Y];
//     res[5]=u/h;
//     res[6]=u-norm[X]*c;
//     res[7]=u+norm[X]*c;
//
//     //third line
//     res[8 ]=-norm[X];
//     res[9 ]=v/h;
//     res[10]=v-norm[Y]*c;
//     res[11]=v+norm[Y]*c;
//
//     //fourth line
//     res[12]=norm[Y]*u-norm[X]*v;
//     res[13]=ek/h;
//     res[14]=h+ek-Vn*c;
//     res[15]=h+ek+Vn*c;
//     return res;
// }

/* ================= 三维（任务 4 新增；与上面 2D 类逐行对应） ================= */

/**
 * @brief 三维欧拉方程的特征系统类
 *
 * 变量序与 2D 类一致（输入/输出均为 [rho, u, v, w, p]），内部按守恒量做变换：
 *   primToChar：把原始量组成守恒量 U=(rho, rho*u, rho*v, rho*w, rho*E)，乘左特征矩阵；
 *   charToPrim：乘右特征矩阵还原守恒量，再拆回原始量。
 * leftEig/rightEig 为 5×5（行主序）互为逆矩阵，行/列构造与 2D 类逐行对应：
 *   剪切波×2（切向基 t1/t2）、熵波、声波±。
 * 说明：leftEig/rightEig 置于 public 区，供内核测试直接取用。
 */
class eigensystemEuler3D
{
    public:
    eigensystemEuler3D(){};

    /**
     * @brief 单态构造：以同一状态构造（内部委托两态构造，等价于不做 Roe 平均）
     */
    eigensystemEuler3D(const std::array<real,5> & prim,const std::array<real,3> & norm_);

    /**
     * @brief 两态构造：左右状态 Roe 平均得到参考态（recon3DFaceCenter 使用）
     */
    eigensystemEuler3D(const std::array<real,5> &priml,const std::array<real,5> &primr,const std::array<real,3> & norm_);

    /**
     * @brief 将原始变量转换为特征变量
     */
    std::array<real,5> primToChar(const std::array<real,5> & prim);

    /**
     * @brief 将特征变量转换为原始变量
     */
    std::array<real,5> charToPrim(const std::array<real,5> & chars);

    std::array<real,25> leftEig;    ///< 左特征矩阵（行主序，5×5）
    std::array<real,25> rightEig;   ///< 右特征矩阵（行主序，5×5）

    private:
    real r;                     ///< 密度
    real u;                     ///< x方向速度
    real v;                     ///< y方向速度
    real w;                     ///< z方向速度
    real p;                     ///< 压力
    real gamma;                 ///< 比热比
    real ek;                    ///< 动能
    real h;                     ///< 比焓
    real c;                     ///< 声速
    real Vn;                    ///< 法向速度
    std::array<real,3> norm;    ///< 法向量
};