/**
 * @file viscousMMS.hpp
 * @brief 粘性 MMS 制造解与其上 NS 通量差的模板实现（任务 4 S9；源项与阶数测试共用）
 *
 * 制造解（平滑、含时、整数波数、[0,2π]³ 周期相容）：
 *   ρ = 1 + 0.1 sin(x+t) cos y cos z
 *   u = 1 + 0.1 sin(x+t) cos y cos z
 *   v = 1 + 0.1 cos x sin(y+t) cos z
 *   w = 1 + 0.1 cos x cos y sin(z+t)
 *   p = 100 + 10 sin(x+t) cos(y+t) cos z
 * 一阶导全部闭式给出；外层对 x/y/z/t 的导数由调用方以复步微分求取。
 * 本文件只含解析函数与四则运算、无按值分支，保证复步（complex-step）精确。
 *
 * 注意：STL 头须先于 macro.hpp 包含（#define real 会破坏 <complex> 的成员名 real）。
 */
#include <complex>
#pragma once
#include "macro.hpp"
#include <cmath>

namespace mms
{
constexpr real Ar=0.1, Au=0.1, Av=0.1, Aw=0.1, Ap=10.0;

// q=[ρ,u,v,w,p]；dq[d][ivar]=∂q_ivar/∂x_d（d=0,1,2）；dt[ivar]=∂/∂t
template<class T>
inline void state(T x,T y,T z,T t,T q[5],T dq[3][5],T dt[5])
{
    T sx=std::sin(x+t), cx=std::cos(x+t);
    T sx0=std::sin(x),  cx0=std::cos(x);
    T sy=std::sin(y),   cy=std::cos(y);
    T sz=std::sin(z),   cz=std::cos(z);
    T syt=std::sin(y+t), cyt=std::cos(y+t);
    T szt=std::sin(z+t), czt=std::cos(z+t);

    q[0]=1.0+Ar*sx*cy*cz;
    q[1]=1.0+Au*sx*cy*cz;
    q[2]=1.0+Av*cx0*syt*cz;
    q[3]=1.0+Aw*cx0*cy*szt;
    q[4]=100.0+Ap*sx*cyt*cz;

    dq[0][0]=Ar*cx*cy*cz;    dq[0][1]=Au*cx*cy*cz;
    dq[0][2]=-Av*sx0*syt*cz; dq[0][3]=-Aw*sx0*cy*szt;
    dq[0][4]=Ap*cx*cyt*cz;

    dq[1][0]=-Ar*sx*sy*cz;   dq[1][1]=-Au*sx*sy*cz;
    dq[1][2]=Av*cx0*cyt*cz;  dq[1][3]=-Aw*cx0*sy*szt;
    dq[1][4]=-Ap*sx*syt*cz;

    dq[2][0]=-Ar*sx*cy*sz;   dq[2][1]=-Au*sx*cy*sz;
    dq[2][2]=-Av*cx0*syt*sz; dq[2][3]=Aw*cx0*cy*czt;
    dq[2][4]=-Ap*sx*cyt*sz;

    dt[0]=Ar*cx*cy*cz;       dt[1]=Au*cx*cy*cz;
    dt[2]=Av*cx0*cyt*cz;     dt[3]=Aw*cx0*cy*czt;
    dt[4]=Ap*(cx*cyt*cz-sx*syt*cz);
}

// d 方向的欧拉通量 Fi 与粘性通量 Fv（组装与 viscousTerm.cpp 同一套公式与符号约定）
template<class T>
inline void fluxDiff(int d,T x,T y,T z,T t,T Fi[5],T Fv[5],real mu,real kappa)
{
    T q[5],dq[3][5],dtv[5];
    state(x,y,z,t,q,dq,dtv);
    T r=q[0],u=q[1],v=q[2],w=q[3],p=q[4];
    T Tt=p/r;
    T Tx=(dq[0][4]-Tt*dq[0][0])/r;
    T Ty=(dq[1][4]-Tt*dq[1][0])/r;
    T Tz=(dq[2][4]-Tt*dq[2][0])/r;

    T E=p/(GAMMA-1.0)+0.5*r*(u*u+v*v+w*w);
    T vel[3]={u,v,w};
    T un=vel[d];
    Fi[0]=r*un;
    Fi[1]=r*u*un+(d==0?p:T(0.0));
    Fi[2]=r*v*un+(d==1?p:T(0.0));
    Fi[3]=r*w*un+(d==2?p:T(0.0));
    Fi[4]=(E+p)*un;

    T div=dq[0][1]+dq[1][2]+dq[2][3];
    T s11=mu*(2.0*dq[0][1]-(2.0/3.0)*div);
    T s22=mu*(2.0*dq[1][2]-(2.0/3.0)*div);
    T s33=mu*(2.0*dq[2][3]-(2.0/3.0)*div);
    T s12=mu*(dq[1][1]+dq[0][2]);
    T s13=mu*(dq[2][1]+dq[0][3]);
    T s23=mu*(dq[2][2]+dq[1][3]);
    T sd0,sd1,sd2,Td;
    if(d==0)      { sd0=s11; sd1=s12; sd2=s13; Td=Tx; }
    else if(d==1) { sd0=s12; sd1=s22; sd2=s23; Td=Ty; }
    else          { sd0=s13; sd1=s23; sd2=s33; Td=Tz; }
    Fv[0]=T(0.0);
    Fv[1]=sd0; Fv[2]=sd1; Fv[3]=sd2;
    Fv[4]=u*sd0+v*sd1+w*sd2+kappa*Td;
}
}
