/**
 * @file diagnostics.cpp
 * @brief 三维诊断输出实现（HIT 开发扩展：增 U2/Tvar/Thvar；方案见 case_for_comment8/HIT/HIT-实施方案.md §5）
 */

#include "diagnostics.hpp"

/**
 * @brief 周期回卷索引（把任意整数索引卷回 [0,n)）
 */
static inline int wrap(int a,int n){ a%=n; return a<0? a+n : a; }

/**
 * @brief 计算三维诊断量（周期域体积分）
 *
 * K：直接按 ρ(u²+v²+w²)/2 · h³ 逐格累加；
 * U2：未加权均方速度 Σuᵢ² · h³；
 * Ω 与散度由四阶中心差分（带周期回卷）求得；Tvar/Thvar 为温度方差（T=p/ρ）与胀量平方积分。
 */
void calcDiagnostics(Data* prim, const std::array<int,3>& icMax, real h, real& K, real& Omega,
                     real& U2, real& Tvar, real& Thvar)
{
    const int nx=icMax[0], ny=icMax[1], nz=icMax[2];
    const real vol=h*h*h;
    K=0.0; Omega=0.0; U2=0.0; Tvar=0.0; Thvar=0.0;

    // 动能 K、未加权均方速度 U2 与温度总和（T=p/ρ）
    real sumT=0.0;
    for(int k=0;k<nz;k++)
    for(int j=0;j<ny;j++)
    for(int i=0;i<nx;i++)
    {
        int id=i+j*nx+k*nx*ny;
        real r=(*prim)(id,0),u=(*prim)(id,1),v=(*prim)(id,2),w=(*prim)(id,3),p=(*prim)(id,4);
        K+=0.5*r*(u*u+v*v+w*w)*vol;
        U2+=(u*u+v*v+w*w)*vol;
        sumT+=p/r;
    }
    const real Tbar=sumT/(real)(nx*ny*nz);

    // 涡度、胀量与温度脉动：四阶中心差分 + 周期回卷
    const real inv=1.0/(12.0*h);
    for(int k=0;k<nz;k++)
    for(int j=0;j<ny;j++)
    for(int i=0;i<nx;i++)
    {
        // 各量在 (a,b,c) 处的取值（lambda 内做周期回卷）
        auto u=[&](int a,int b,int c){ return (*prim)(wrap(a,nx)+wrap(b,ny)*nx+wrap(c,nz)*nx*ny,1); };
        auto v=[&](int a,int b,int c){ return (*prim)(wrap(a,nx)+wrap(b,ny)*nx+wrap(c,nz)*nx*ny,2); };
        auto w=[&](int a,int b,int c){ return (*prim)(wrap(a,nx)+wrap(b,ny)*nx+wrap(c,nz)*nx*ny,3); };
        auto T=[&](int a,int b,int c){ int id=wrap(a,nx)+wrap(b,ny)*nx+wrap(c,nz)*nx*ny; return (*prim)(id,4)/(*prim)(id,0); };

        real dudx=(u(i-2,j,k)-8.0*u(i-1,j,k)+8.0*u(i+1,j,k)-u(i+2,j,k))*inv;
        real dvdy=(v(i,j-2,k)-8.0*v(i,j-1,k)+8.0*v(i,j+1,k)-v(i,j+2,k))*inv;
        real dwdz=(w(i,j,k-2)-8.0*w(i,j,k-1)+8.0*w(i,j,k+1)-w(i,j,k+2))*inv;
        real dwdy=(w(i,j-2,k)-8.0*w(i,j-1,k)+8.0*w(i,j+1,k)-w(i,j+2,k))*inv;
        real dvdz=(v(i,j,k-2)-8.0*v(i,j,k-1)+8.0*v(i,j,k+1)-v(i,j,k+2))*inv;
        real dudz=(u(i,j,k-2)-8.0*u(i,j,k-1)+8.0*u(i,j,k+1)-u(i,j,k+2))*inv;
        real dwdx=(w(i-2,j,k)-8.0*w(i-1,j,k)+8.0*w(i+1,j,k)-w(i+2,j,k))*inv;
        real dvdx=(v(i-2,j,k)-8.0*v(i-1,j,k)+8.0*v(i+1,j,k)-v(i+2,j,k))*inv;
        real dudy=(u(i,j-2,k)-8.0*u(i,j-1,k)+8.0*u(i,j+1,k)-u(i,j+2,k))*inv;

        real wx=dwdy-dvdz;
        real wy=dudz-dwdx;
        real wz=dvdx-dudy;
        Omega+=0.5*(wx*wx+wy*wy+wz*wz)*vol;

        real div=dudx+dvdy+dwdz;
        Thvar+=div*div*vol;

        real Tp=T(i,j,k)-Tbar;
        Tvar+=Tp*Tp*vol;
    }
}
