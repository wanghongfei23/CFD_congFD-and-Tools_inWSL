/**
 * @file sp_flux.cpp
 * @brief 空间离散通量计算实现文件
 * 
 * 该文件实现了各种空间离散格式下的数值通量计算方法
 */

#include "SpaceDis.hpp"
#include "fluxScheme.hpp"
#include "interScheme.hpp"

/**
 * @brief 计算对流项通量
 * @param i 网格点索引
 * 
 * 适用于方程 u_t + a * u_x == 0，使用Lax-Friedrichs数值通量公式计算界面通量
 */
void SpaceDis::calFluxConv(int i)
{
    /*for u_t + a * u_x == 0*/
    real a = 1.0;  // 对流速度
    real ul, ur;
    // 使用重构方法获取左右界面值
    ul = (this->*reconLMethod)(i).at(0);
    ur = (this->*reconRMethod)(i).at(0);
    
    // 另一种简单的通量计算方式（根据特征线方向选择）
    /*
    ul=data(i-1,0);
    ur=data(i,0);*/
    /*if (a>0)
    {
        flux(i,0)=a*ul;
    }
    else
    {
        flux(i,0)=a*ur;
    }*/

    // 使用Lax-Friedrichs数值通量公式计算界面通量
    // 公式: 0.5*(a*ul + a*ur - |a|*(ur-ul))
    // 这是一种守恒型的数值通量格式，具有良好的稳定性
    fluxAt(i, 0) = 0.5 * (a * ul + a * ur - std::abs(a) * (ur - ul));
}

/**
 * @brief 计算精度测试问题的通量
 * @param i 网格点索引
 * 
 * 主要用于验证代码的精度和正确性，直接使用左界面重构值作为通量
 */
void SpaceDis::calFluxAccuracyTest(int i)
{
    /*for u_t + a * u_x == 0*/
    real a = 1.0;  // 波速
    real ul, ur;
    // 获取左界面重构值
    ul = (this->*reconLMethod)(i).at(0);
    // 精度测试中直接使用左界面值作为通量
    // 这种简单的通量计算方式便于分析数值方法的精度特性
    fluxAt(i, 0) = ul;
}

/**
 * @brief 计算Burgers方程的通量
 * @param i 网格点索引
 * 
 * 方程: u_t + (u^2/2)_x = 0，使用Lax-Friedrichs格式计算非线性对流项的数值通量
 */
void SpaceDis::calFluxBurgers(int i)
{
    /*for u_t + a * u_x == 0*/

    real ul, ur, aLF;
    // 获取数据中的最大元素值作为LF通量的速度参数
    // 这是为了保证数值格式的CFL条件稳定
    aLF = data->maxElement(0);

    // 获取左右界面的重构值
    ul = (this->*reconLMethod)(i).at(0);
    ur = (this->*reconRMethod)(i).at(0);

    // ul=data(i-1,0);
    // ur=data(i,0);
    // L-F flux
    real al, ar;
    // 计算左右界面的绝对值
    al = std::abs(ul);
    ar = std::abs(ur);
    // real a=std::max(ar,ar);
    // real a=(al+ar)/2;
    real a = aLF;  // 使用最大速度作为通量参数
    // 使用Lax-Friedrichs格式计算Burgers方程的界面通量
    // 公式: 0.5*(f(ul)+f(ur)-a*(ur-ul)), 其中f(u)=u^2/2
    fluxAt(i, 0) = 0.5 * (ul * ul / 2 + ur * ur / 2 - a * (ur - ul));
}

typedef std::array<real, 2> arr2;

// 计算一维欧拉方程的界面通量函数（旧版实现方式，已被注释）
// void SpaceDis::calFluxEuler1D(int i)
// {
//     /*for u_t + a * u_x == 0*/
//     auto WL=(this->*reconLMethod)(i);
//     auto WR=(this->*reconRMethod)(i);
//     enum{R,U,P};
//     if (WL[P]<0||WL[R]<0)
//     {
//         WL[R]=at(i-1,R);
//         WL[U]=at(i-1,U);
//         WL[P]=at(i-1,P);
//         std::cout<<"SpaceDis error: negative pressure/Density \n";
//     }
//     if (WR[P]<0||WR[R]<0)
//     {
//         WR[R]=at(i,R);
//         WR[U]=at(i,U);
//         WR[P]=at(i,P);
//         std::cout<<"SpaceDis error: negative pressure/Density \n";
//     }
//     std::vector<real> iflux=roeFlux1D2(WL[0],WR[0],WL[1],WR[1],WL[2],WR[2]);

//     //auto iflux=HLLCFlux1D(WL[0],WR[0],WL[1],WR[1],WL[2],WR[2]);
//     //std::vector<real> iflux2=roeFlux1D(r,u,p,H,RT);
//     for (int ivar = 0; ivar < 3; ivar++)
//     {
//         fluxAt(i,ivar)=iflux[ivar];
//     }
// }

/**
 * @brief 计算一维欧拉方程的界面通量
 * @param i 网格点索引
 */
void SpaceDis::calFluxEuler1D(int i)
{
    /*for u_t + a * u_x == 0*/
    // 使用一阶重构方法获取界面左右状态
    auto W = this->recon1DFaceCenter(i);
    // 定义变量索引常量
    enum { R,  // 密度
        U,     // 速度
        P      // 压力
    };
    
    // 检查左界面状态是否合理（密度和压力必须为正）
    if (W[0] < 0 || W[2] < 0) {
        W[0] = at(i - 1, R);
        W[1] = at(i - 1, U);
        W[2] = at(i - 1, P);
        std::cout << "SpaceDis error: negative pressure/Density \n";
    }
    
    // 检查右界面状态是否合理（密度和压力必须为正）
    if (W[3] < 0 || W[5] < 0) {
        W[3] = at(i, R);
        W[4] = at(i, U);
        W[5] = at(i, P);
        std::cout << "SpaceDis error: negative pressure/Density \n";
    }
    
    // 使用Roe平均方法计算欧拉方程的界面通量
    std::vector<real> iflux = roeFlux1D2(W[0], W[3], W[1], W[4], W[2], W[5]);

    // 另一种可选的HLLC通量计算方法（被注释）
    // auto iflux=HLLCFlux1D(W[0],W[3],W[1],W[4],W[2],W[5]);
    // std::vector<real> iflux2=roeFlux1D(r,u,p,H,RT);
    // 将计算得到的通量赋值给各守恒变量
    for (int ivar = 0; ivar < 3; ivar++) {
        fluxAt(i, ivar) = iflux[ivar];
    }
}

/**
 * @brief 计算二维欧拉方程的界面通量
 * @param i 网格点索引
 */
void SpaceDis::calFluxEuler2D(int i)
{
    /*for u_t + a * u_x == 0*/

    // auto WL=(this->*reconLMethod)(i);
    // auto WR=(this->*reconRMethod)(i);
    // 使用二阶重构方法获取界面左右状态
    auto W = this->recon2DFaceCenter(i);
    
    // 检查左界面状态的合理性（密度和压力必须为正）
    // 注：这里的检查包含了更复杂的调试信息输出代码，但现在被简化了
    if (W[3] < 0 || W[0] < 0) {
        // int ivar=(WL[3]<0||isnan(WL[3]))? 3:0;
        // std::array<real,5> q={at(i-3,ivar),at(i-2,ivar),at(i-1,ivar),at(i,ivar),at(i+1,ivar)};
        // std::array<real,3> beta;
        // unsigned flag=0;
        // beta[0]= 1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
        //     + 1.0/4.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);
        // beta[1]= 1.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
        //     + 1.0/4.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);

        // beta[2]= 1.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
        //     + 1.0/4.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);

        // auto minBeta=std::min_element(beta.begin(),beta.end())-beta.begin();
        // //if (beta[minBeta]<eps) return q[2];
        // real CT=CTi,sumGamma=0;
        // real tau=std::abs(beta[2]-beta[0]),KK=CT*(beta[minBeta]+tau);
        // real eps=1e-8;
        // if (minBeta==0||(beta[0]+tau)*beta[minBeta]>KK*beta[0]+eps)
        // {
        //     flag+=1;
        // }
        // if (minBeta==1||(beta[1]+tau)*beta[minBeta]>KK*beta[1]+eps)
        // {
        //     flag+=2;
        // }
        // if (minBeta==2||(beta[2]+tau)*beta[minBeta]>KK*beta[2]+eps)
        // {
        //     flag+=4;
        // }

        // std::cout<<std::format("SpaceDis error L: negative pressure i={} WLR={} WRR={} WLP={} WRP={} \
        //                         WSR={:.4f} {:.4f} {:.4f} {:.4f} {:.4f} \
        //                         WSU={:.4f} {:.4f} {:.4f} {:.4f} {:.4f} \
        //                         WSV={:.4f} {:.4f} {:.4f} {:.4f} {:.4f} \
        //                         WSP={:.4f} {:.4f} {:.4f} {:.4f} {:.4f} \
        //                         beta={:.4f} {:.4f} {:.4f} flag={}\n"
        //                       ,i,WL[0],WR[0],WL[3],WR[3]
        //                       ,at(i-3,0),at(i-2,0),at(i-1,0),at(i,0),at(i+1,0)
        //                       ,at(i-3,1),at(i-2,1),at(i-1,1),at(i,1),at(i+1,1)
        //                       ,at(i-3,2),at(i-2,2),at(i-1,2),at(i,2),at(i+1,2)
        //                       ,at(i-3,3),at(i-2,3),at(i-1,3),at(i,3),at(i+1,3)
        //                       ,beta[0],beta[1],beta[2],flag);
        std::cout << std::format("SpaceDis error L: negative pressure i={}\n", i);
        // 如果不合理，则使用原始数据点的值进行修正
        W = { at(i - 1, 0), at(i - 1, 1), at(i - 1, 2), at(i - 1, 3), W[4], W[5], W[6], W[7] };
        auto WX = (this->*reconLMethod)(i);
    }
    
    // 检查右界面状态的合理性（密度和压力必须为正）
    // 注：这里的检查也包含了复杂调试信息输出代码，但现在被简化了
    if (W[3 + 4] < 0 || W[0 + 4] < 0) //||isnan(WR[3])||isnan(WR[0]))
    {
        // std::array<real,5> q={at(i+2,3),at(i+1,3),at(i,3),at(i-1,3),at(i-2,3)};
        // std::array<real,3> beta;
        // beta[0]= 1.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
        //     + 1.0/4.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);
        // beta[1]= 1.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
        //     + 1.0/4.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);
        // beta[2]= 1.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
        //     + 1.0/4.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);
        // std::cout<<std::format("SpaceDis error R: negative pressure i={} WLR={} WRR={} WLP={} WRP={} WS={:.4f} {:.4f} {:.4f} {:.4f} {:.4f} beta={:.4f} {:.4f} {:.4f}\n"
        //                       ,i,WL[0],WR[0],WL[3],WR[3],at(i-2,3),at(i-1,3),at(i,3),at(i+1,3),at(i+2,3),beta[0],beta[1],beta[2]);
        std::cout << std::format("SpaceDis error R: negative pressure i={}\n", i);
        // 如果不合理，则使用原始数据点的值进行修正
        W = { W[0], W[1], W[2], W[3], at(i, 0), at(i, 1), at(i, 2), at(i, 3) };
        // auto WX=(this->*reconRMethod)(i);
    }

    // 根据算例类型选择通量方法
    // 只有当是二维Double Mach算例(dim==2且nCase==4)时使用HLLC通量，其他算例使用Roe通量
    std::array<real, 4> iflux;
    if (info->dim == 2 && info->nCase == 4) {
        // 使用HLLC通量计算方法（针对Double Mach算例）
        iflux = HLLCFlux2D2(W[0], W[4], W[1], W[5], W[2], W[6], W[3], W[7], norm);
    } else {
        // 使用Roe平均方法计算二维欧拉方程的界面通量（默认方法）
        iflux = roeFlux2DSym(W[0], W[4], W[1], W[5], W[2], W[6], W[3], W[7], norm);
    }
    
    // std::vector<real> iflux2=roeFlux1D(r,u,p,H,RT);
    // 将计算得到的通量赋值给各守恒变量
    for (int ivar = 0; ivar < 4; ivar++) {
        fluxAt(i, ivar) = iflux[ivar];
    }
}

/**
 * @brief 计算三维欧拉方程的界面通量（任务 5 新增；镜像 calFluxEuler2D）
 * @param i 网格点索引（面索引）
 *
 * W 前 5 值为左状态 [rho,u,v,w,p]、后 5 值为右状态；
 * 负压/负密度时回退为相邻格值（含格间简单回退，与 2D 版同型）。
 */
void SpaceDis::calFluxEuler3D(int i)
{
    auto W = this->recon3DFaceCenter(i);

    if (W[4] < 0 || W[0] < 0) {
        std::cout << "SpaceDis error L: negative pressure i=" << i << "\n";
        W = { at(i-1,0), at(i-1,1), at(i-1,2), at(i-1,3), at(i-1,4), W[5], W[6], W[7], W[8], W[9] };
    }
    if (W[9] < 0 || W[5] < 0) {
        std::cout << "SpaceDis error R: negative pressure i=" << i << "\n";
        W = { W[0], W[1], W[2], W[3], W[4], at(i,0), at(i,1), at(i,2), at(i,3), at(i,4) };
    }

    std::array<real, 5> iflux = roeFlux3DSym(W[0],W[5],W[1],W[6],W[2],W[7],W[3],W[8],W[4],W[9],norm);
    for (int ivar = 0; ivar < 5; ivar++) {
        fluxAt(i, ivar) = iflux[ivar];
    }
}
