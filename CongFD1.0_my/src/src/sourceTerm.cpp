/**
 * @file sourceTerm.cpp
 * @brief 源项类的实现文件
 */

#include <complex>
#include <cstdio>
#include <cstdlib>
#include "SourceTerm.hpp"
#include "viscousMMS.hpp"

// MMS 源项符号：rhs_total = −∂_t q*（主循环 cons=temp−dt·rhs 约定）
// ⇒ rhs 中的源项贡献 = −[∂_t q* + ∇·F_inv − ∇·F_vis]（推导见 04 实施分解 S9；由阶数测试裁决）
static const real MMS_SIGN=-1.0;

/**
 * @brief 计算重力源项
 * 
 * 该函数专门用于计算Rayleigh-Taylor不稳定性问题中的重力源项。
 */
void SourceTerm::calGravitySource()
{
    //here it is still only for R-T instability case
    if(nprim!=4) 
    {
        std::cout<<"SourceTerm error: nprim incorrect\n";
    }
    real r,u,v,p;
    for (int i = 0; i < n; i++)
    {
        r=(*prim)(i,0);
        u=(*prim)(i,1);
        v=(*prim)(i,2);
        p=(*prim)(i,3);
        (*rhs)(i,2)-=r;
        (*rhs)(i,3)-=r*v;
    }
    
}

/**
 * @brief 计算粘性 MMS 造解源项（任务 4 S9；测试用）
 *
 * S = ∂_t q* + ∇·F_inv(q*) − ∇·F_vis(q*,∇q*)（对制造解的 NS 残差），
 * 外层 x/y/z/t 导数全部由复步微分求取（机器精度、免手工展开；q* 与 ∇q* 闭式见 viscousMMS.hpp）。
 */
void SourceTerm::calViscousMMS()
{
    if(nprim!=5)
    {
        std::cout<<"SourceTerm error: MMS requires nprim=5 (3D Euler)\n";
        return;
    }
    auto iM=info->icMax();
    real h0=(info->calZone[1]-info->calZone[0])/(info->iMax[0]-1);
    real h1=(info->calZone[3]-info->calZone[2])/(info->iMax[1]-1);
    real h2=(info->calZone[5]-info->calZone[4])/(info->iMax[2]-1);
    real mu=1.0/info->Re;
    real kappa=mu*(GAMMA/(GAMMA-1.0))/info->Pr;
    const real eps=1e-28;

    #pragma omp parallel for
    for(int k=0;k<iM[2];k++)
    for(int j=0;j<iM[1];j++)
    for(int i=0;i<iM[0];i++)
    {
        int idx=i+j*iM[0]+k*iM[0]*iM[1];
        real x=info->calZone[0]+(i+0.5)*h0;
        real y=info->calZone[2]+(j+0.5)*h1;
        real z=info->calZone[4]+(k+0.5)*h2;
        real t=info->t;
        real S[5]={0,0,0,0,0};
        for(int d=0;d<3;d++)
        {
            std::complex<real> X(x,0.0),Y(y,0.0),Z(z,0.0);
            if(d==0) X=std::complex<real>(x,eps);
            if(d==1) Y=std::complex<real>(y,eps);
            if(d==2) Z=std::complex<real>(z,eps);
            std::complex<real> Fi[5],Fv[5];
            mms::fluxDiff< std::complex<real> >(d,X,Y,Z,std::complex<real>(t,0.0),Fi,Fv,mu,kappa);
            for(int c=0;c<5;c++) S[c]+=std::imag(Fi[c]-Fv[c])/eps;
        }
        {
            std::complex<real> q[5],dq[3][5],dtv[5];
            mms::state< std::complex<real> >(std::complex<real>(x,0.0),std::complex<real>(y,0.0),
                std::complex<real>(z,0.0),std::complex<real>(t,eps),q,dq,dtv);
            // ∂_t 须作用于守恒量：先由原始量 q=[ρ,u,v,w,p] 代换后取虚部
            // （修正记录：初版误取原始量时间导数；质量行恰等价而不显，动量/能量行系统性偏）
            std::complex<real> r=q[0],u=q[1],v=q[2],w=q[3],pr=q[4];
            std::complex<real> cc[5];
            cc[0]=r;
            cc[1]=r*u;
            cc[2]=r*v;
            cc[3]=r*w;
            cc[4]=pr/(GAMMA-1.0)+0.5*r*(u*u+v*v+w*w);
            for(int c=0;c<5;c++) S[c]+=std::imag(cc[c])/eps;
        }
        {
            static bool dbgDone=false;
            if(!dbgDone && std::getenv("MMS_DEBUG") && i==1 && j==1 && k==1)
            {
                dbgDone=true;
                std::printf("MMSDBG t=%.15e x=%.15e y=%.15e z=%.15e\n",t,x,y,z);
                std::printf("MMSDBG S:");
                for(int c=0;c<5;c++) std::printf(" %.12e",S[c]);
                std::printf("\n");
            }
        }
        for(int c=0;c<5;c++) (*rhs)(idx,c)+=MMS_SIGN*S[c];
    }
}

/**
 * @brief 空操作函数
 *
 * 当没有源项需要计算时调用此函数。
 */
void SourceTerm::nothingHappened()
{

}


/**
 * @brief 构造函数
 * @param prim_ 原始变量数据指针
 * @param rhs_ 右手端数据指针
 * @param info_ Info对象指针
 */
SourceTerm::SourceTerm(Data* prim_,Data* rhs_,Info* info_)
{
    prim=prim_;
    rhs=rhs_;
    info=info_;

    auto iMax=info->icMax();
    n=iMax[0]*iMax[1]*iMax[2];
    nprim=info->nPrim();
    nCons=info->nCons();
    switch (info->sourceType)
    {
    case GRAVITY:
        calSourceMethod=(&SourceTerm::calGravitySource);
        break;
    case VISCOUSMMS:
        calSourceMethod=(&SourceTerm::calViscousMMS);
        break;

    default:
        calSourceMethod=(&SourceTerm::nothingHappened);
        break;
    }
}

/**
 * @brief 计算源项
 * 
 * 根据初始化时设定的源项类型，调用相应的源项计算函数。
 */
void SourceTerm::calSource()
{
    (this->*calSourceMethod)();
}