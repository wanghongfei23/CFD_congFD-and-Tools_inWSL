/**
 * @file diagnostics.cpp
 * @brief 三维诊断输出实现（任务 7 新增）
 */

#include "diagnostics.hpp"

/**
 * @brief 周期回卷索引（把任意整数索引卷回 [0,n)）
 */
static inline int wrap(int a,int n){ a%=n; return a<0? a+n : a; }

/**
 * @brief 计算三维动能 K 与涡度 Ω（周期域体积分）
 *
 * K：直接按 ρ(u²+v²+w²)/2 · h³ 逐格累加；
 * Ω：涡量由四阶中心差分（带周期回卷）求得，再逐格累加 |ω|²/2 · h³。
 * 用于与 Fu et al. 2018 图 24 口径的 K(t)/Ω(t) 时间序列比对。
 */
void calcDiagnostics(Data* prim, const std::array<int,3>& icMax, real h, real& K, real& Omega)
{
    const int nx=icMax[0], ny=icMax[1], nz=icMax[2];
    const real vol=h*h*h;
    K=0.0; Omega=0.0;

    // 动能：K = Σ ρ(u²+v²+w²)/2 · h³
    for(int k=0;k<nz;k++)
    for(int j=0;j<ny;j++)
    for(int i=0;i<nx;i++)
    {
        int id=i+j*nx+k*nx*ny;
        real r=(*prim)(id,0),u=(*prim)(id,1),v=(*prim)(id,2),w=(*prim)(id,3);
        K+=0.5*r*(u*u+v*v+w*w)*vol;
    }

    // 涡度：四阶中心差分 + 周期回卷；Ω = Σ (ωx²+ωy²+ωz²)/2 · h³
    const real inv=1.0/(12.0*h);
    for(int k=0;k<nz;k++)
    for(int j=0;j<ny;j++)
    for(int i=0;i<nx;i++)
    {
        // 三个速度分量在 (a,b,c) 处的取值（lambda 内做周期回卷）
        auto u=[&](int a,int b,int c){ return (*prim)(wrap(a,nx)+wrap(b,ny)*nx+wrap(c,nz)*nx*ny,1); };
        auto v=[&](int a,int b,int c){ return (*prim)(wrap(a,nx)+wrap(b,ny)*nx+wrap(c,nz)*nx*ny,2); };
        auto w=[&](int a,int b,int c){ return (*prim)(wrap(a,nx)+wrap(b,ny)*nx+wrap(c,nz)*nx*ny,3); };

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
    }
}
