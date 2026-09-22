/**
 * @file initializer.cpp
 * @brief 初始化器类的实现文件
 */

#include"initializer.hpp"
#include <random>
#include <cstdio>

// =====================================================================
// HIT 初值谱合成（Johnsen et al. 2010 附录 A 式 (19)-(21)）
// 口径与方案见 case_for_comment8/HIT/（HIT-算例定义.md §3.1、HIT-实施方案.md §3）。
// 说明：极简复数与基-2 FFT 自含于此（不复用 std::complex：macro.hpp 的
// real 宏会破坏其成员名）；生成串行且逐模式确定，结果与线程数无关。
// =====================================================================
namespace
{
struct Cpx { real re; real im; };
inline Cpx cpxAdd(Cpx a,Cpx b){ return Cpx{a.re+b.re,a.im+b.im}; }
inline Cpx cpxSub(Cpx a,Cpx b){ return Cpx{a.re-b.re,a.im-b.im}; }
inline Cpx cpxMul(Cpx a,Cpx b){ return Cpx{a.re*b.re-a.im*b.im,a.re*b.im+a.im*b.re}; }

// 一维基-2 FFT（原位）；inv=false 正变换 e^{-i2πjk/n}，inv=true 逆变换 e^{+i2πjk/n}；均不归一
inline void fft1d(std::vector<Cpx>& a,bool inv)
{
    const int n=(int)a.size();
    for(int i=1,j=0;i<n;i++)
    {
        int bit=n>>1;
        for(;j&bit;bit>>=1) j^=bit;
        j^=bit;
        if(i<j){ Cpx t=a[i]; a[i]=a[j]; a[j]=t; }
    }
    const real PI=3.14159265358979323846;
    for(int len=2;len<=n;len<<=1)
    {
        const real ang=(inv? 2.0:-2.0)*PI/(real)len;
        const real wr=std::cos(ang), wi=std::sin(ang);
        for(int s=0;s<n;s+=len)
        {
            real cr=1.0, ci=0.0;
            for(int j=0;j<len/2;j++)
            {
                const Cpx u=a[s+j];
                const Cpx v=cpxMul(Cpx{cr,ci},a[s+j+len/2]);
                a[s+j]=cpxAdd(u,v);
                a[s+j+len/2]=cpxSub(u,v);
                const real ncr=cr*wr-ci*wi, nci=cr*wi+ci*wr;
                cr=ncr; ci=nci;
            }
        }
    }
}

// 三维变换（逐向一维 FFT；数组按 idx=i+j*nx+k*nx*ny 存储）
inline void fft3dRun(std::vector<Cpx>& f,int nx,int ny,int nz,bool inv)
{
    std::vector<Cpx> line;
    line.resize(nx);
    for(int k=0;k<nz;k++)
    for(int j=0;j<ny;j++)
    {
        const size_t base=(size_t)j*nx+(size_t)k*nx*ny;
        for(int i=0;i<nx;i++) line[i]=f[base+i];
        fft1d(line,inv);
        for(int i=0;i<nx;i++) f[base+i]=line[i];
    }
    line.resize(ny);
    for(int k=0;k<nz;k++)
    for(int i=0;i<nx;i++)
    {
        const size_t base=(size_t)i+(size_t)k*nx*ny;
        for(int j=0;j<ny;j++) line[j]=f[base+(size_t)j*nx];
        fft1d(line,inv);
        for(int j=0;j<ny;j++) f[base+(size_t)j*nx]=line[j];
    }
    line.resize(nz);
    for(int j=0;j<ny;j++)
    for(int i=0;i<nx;i++)
    {
        const size_t base=(size_t)i+(size_t)j*nx;
        for(int k=0;k<nz;k++) line[k]=f[base+(size_t)k*nx*ny];
        fft1d(line,inv);
        for(int k=0;k<nz;k++) f[base+(size_t)k*nx*ny]=line[k];
    }
}

// HIT 初值：随机螺线管速度场谱合成；返回前打印自检（无散、厄米、往返、归一、壳层谱）
//   幅度：每模式 |û(k)|²=2E(k)/(4πk²)，E(k)=16√(2/π)·k⁴/k₀⁵·exp(−2k²/k₀²)
//   相位：φ1,φ2,φ3~U[0,2π)，每波数三元组重抽；固定种子（64 位原始位映射，确定性）
//   实场：半空间抽样＋厄米对称；Nyquist 面（|k_i|=N/2，振幅 ~e-128）与 (0,0,0) 置零
//   归一：逆变换后按离散场实测 u_rms 全局归一至目标值
inline void hitInitialVelocity(std::vector<real>& u,std::vector<real>& v,std::vector<real>& w,
                               int nx,int ny,int nz,real k0,real urms,unsigned long long seed)
{
    const size_t n=(size_t)nx*ny*nz;
    const real PI=3.14159265358979323846;

    auto isPow2=[](int m){ return m>0 && ((m&(m-1))==0); };
    if(!isPow2(nx)||!isPow2(ny)||!isPow2(nz))
    {
        std::printf("[HIT-IC] error: grid must be powers of two, got %d,%d,%d\n",nx,ny,nz);
        return;
    }

    std::vector<Cpx> uh[3];
    for(int c=0;c<3;c++) uh[c].assign(n,Cpx{0.0,0.0});

    std::mt19937_64 rng(seed);
    auto uni01=[&rng](){ return (real)(rng()>>11)*0x1.0p-53; };            // [0,1)
    auto Ek=[&](real kk){ const real k2=kk*kk; return 16.0*std::sqrt(2.0/PI)*k2*k2/(k0*k0*k0*k0*k0)*std::exp(-2.0*k2/(k0*k0)); };
    auto kidx=[](int i,int m){ return (i<=m/2)? i : i-m; };                // 0..m−1 → −m/2..m/2−1

    for(int kz=0;kz<nz;kz++)
    {
        const int KZ=kidx(kz,nz);
        if(KZ==nz/2||KZ==-nz/2) continue;
        for(int ky=0;ky<ny;ky++)
        {
            const int KY=kidx(ky,ny);
            if(KY==ny/2||KY==-ny/2) continue;
            for(int kx=0;kx<nx;kx++)
            {
                const int KX=kidx(kx,nx);
                if(KX==nx/2||KX==-nx/2) continue;
                if(KX==0&&KY==0&&KZ==0) continue;
                if(!(KZ>0||(KZ==0&&KY>0)||(KZ==0&&KY==0&&KX>0))) continue; // 半空间

                const real k2=(real)KX*KX+(real)KY*KY+(real)KZ*KZ;
                const real kk=std::sqrt(k2);
                const real k12=std::sqrt((real)KX*KX+(real)KY*KY);
                const real r=std::sqrt(2.0*Ek(kk)/(4.0*PI*k2));            // 幅度 √(2E/(4πk²))

                const real ph1=2.0*PI*uni01(), ph2=2.0*PI*uni01(), ph3=2.0*PI*uni01();
                const real ca=r*std::cos(ph3), sb=r*std::sin(ph3);         // a=r e^{iφ1}cosφ3，b=r e^{iφ2}sinφ3
                const real ar=ca*std::cos(ph1), ai=ca*std::sin(ph1);
                const real br=sb*std::cos(ph2), bi=sb*std::sin(ph2);

                real P1x,P1y,P2x,P2y,P2z;
                if(k12>0.0)
                {
                    P1x=(real)KY/k12; P1y=-(real)KX/k12;
                    P2x=(real)KX*(real)KZ/(k12*kk); P2y=(real)KY*(real)KZ/(k12*kk); P2z=-k12/kk;
                }
                else
                {
                    P1x=1.0; P1y=0.0;                                      // k₁₂=0 约定（Johnsen 附录 A）
                    P2x=0.0; P2y=(real)KZ/kk; P2z=0.0;
                }
                const size_t idx=(size_t)kx+(size_t)ky*nx+(size_t)kz*nx*ny;
                uh[0][idx]=Cpx{ar*P1x+br*P2x,ai*P1x+bi*P2x};               // û = a·P1 + b·P2
                uh[1][idx]=Cpx{ar*P1y+br*P2y,ai*P1y+bi*P2y};
                uh[2][idx]=Cpx{br*P2z,          bi*P2z};

                const int mx=(nx-kx)%nx, my=(ny-ky)%ny, mz=(nz-kz)%nz;     // 厄米对称
                const size_t mi=(size_t)mx+(size_t)my*nx+(size_t)mz*nx*ny;
                uh[0][mi]=Cpx{uh[0][idx].re,-uh[0][idx].im};
                uh[1][mi]=Cpx{uh[1][idx].re,-uh[1][idx].im};
                uh[2][mi]=Cpx{uh[2][idx].re,-uh[2][idx].im};
            }
        }
    }

    // 无散自检（k·û=0，构造保证）
    real kdot=0.0,knorm=0.0;
    for(int kz=0;kz<nz;kz++)
    for(int ky=0;ky<ny;ky++)
    for(int kx=0;kx<nx;kx++)
    {
        const size_t idx=(size_t)kx+(size_t)ky*nx+(size_t)kz*nx*ny;
        const real KX=(real)kidx(kx,nx), KY=(real)kidx(ky,ny), KZ=(real)kidx(kz,nz);
        const Cpx U1=uh[0][idx],U2=uh[1][idx],U3=uh[2][idx];
        const real dr=KX*U1.re+KY*U2.re+KZ*U3.re, di=KX*U1.im+KY*U2.im+KZ*U3.im;
        const real um=std::sqrt(U1.re*U1.re+U1.im*U1.im+U2.re*U2.re+U2.im*U2.im+U3.re*U3.re+U3.im*U3.im);
        const real km=std::sqrt(KX*KX+KY*KY+KZ*KZ);
        kdot=std::max(kdot,std::sqrt(dr*dr+di*di));
        knorm=std::max(knorm,km*um);
    }

    // 逆变换 + 实部提取 + 往返自检（分量 0；前向变换应得 n·û）
    std::vector<real> uu(n),vv(n),ww(n);
    real imMax=0.0,rtDev=0.0,rtNorm=0.0;
    for(int c=0;c<3;c++)
    {
        std::vector<Cpx> f=uh[c];
        fft3dRun(f,nx,ny,nz,true);
        for(size_t i=0;i<n;i++)
        {
            imMax=std::max(imMax,std::abs(f[i].im));
            f[i].im=0.0;                                                    // 实场（虚部为舍入噪声）
        }
        if(c==0)
        {
            std::vector<Cpx> g=f;
            fft3dRun(g,nx,ny,nz,false);
            const real nn=(real)n;
            for(size_t i=0;i<n;i++)
            {
                const real dr=g[i].re-nn*uh[0][i].re, di=g[i].im-nn*uh[0][i].im;
                rtDev=std::max(rtDev,std::sqrt(dr*dr+di*di));
                rtNorm=std::max(rtNorm,nn*std::sqrt(uh[0][i].re*uh[0][i].re+uh[0][i].im*uh[0][i].im));
            }
        }
        for(size_t i=0;i<n;i++)
        {
            if(c==0) uu[i]=f[i].re;
            else if(c==1) vv[i]=f[i].re;
            else ww[i]=f[i].re;
        }
    }

    // 全局归一到目标 u_rms
    real sum2=0.0;
    for(size_t i=0;i<n;i++) sum2+=uu[i]*uu[i]+vv[i]*vv[i]+ww[i]*ww[i];
    const real rms0=std::sqrt(sum2/(3.0*(real)n));
    const real sc=(rms0>0.0)? urms/rms0 : 0.0;
    u.resize(n); v.resize(n); w.resize(n);
    real maxu=0.0,sum2n=0.0;
    for(size_t i=0;i<n;i++)
    {
        u[i]=uu[i]*sc; v[i]=vv[i]*sc; w[i]=ww[i]*sc;
        sum2n+=u[i]*u[i]+v[i]*v[i]+w[i]*w[i];
        maxu=std::max(maxu,std::sqrt(u[i]*u[i]+v[i]*v[i]+w[i]*w[i]));
    }
    const real rms1=std::sqrt(sum2n/(3.0*(real)n));

    std::printf("[HIT-IC] grid=%dx%dx%d seed=%llu k0=%g u_rms target=%g\n",nx,ny,nz,seed,(double)k0,(double)urms);
    std::printf("[HIT-IC] hermitian max|Im(iFFT)|=%.3e\n",(double)imMax);
    std::printf("[HIT-IC] roundtrip max rel dev=%.3e\n",(double)(rtNorm>0.0? rtDev/rtNorm : 0.0));
    std::printf("[HIT-IC] solenoidality max|k.u|=%.3e (scale %.3e)\n",(double)kdot,(double)knorm);
    std::printf("[HIT-IC] u_rms raw=%.16e normalized=%.16e\n",(double)rms0,(double)rms1);
    std::printf("[HIT-IC] max|u|=%.6e\n",(double)maxu);
    std::printf("[HIT-IC] shell spectrum E_num = sc^2*shell-sum(0.5|u_hat|^2)*4*pi*k^2/N:\n");
    std::printf("[HIT-IC]   k     N      E_num          E_target\n");
    for(int ks=1;ks<=10;ks++)
    {
        real shellSum=0.0; int cnt=0;
        for(int kz=0;kz<nz;kz++)
        {
            const int KZ=kidx(kz,nz);
            for(int ky=0;ky<ny;ky++)
            {
                const int KY=kidx(ky,ny);
                for(int kx=0;kx<nx;kx++)
                {
                    const int KX=kidx(kx,nx);
                    const real kk=std::sqrt((real)(KX*KX+KY*KY+KZ*KZ));
                    if(kk<(real)ks-0.5 || kk>=(real)ks+0.5) continue;
                    const size_t idx=(size_t)kx+(size_t)ky*nx+(size_t)kz*nx*ny;
                    const Cpx U1=uh[0][idx],U2=uh[1][idx],U3=uh[2][idx];
                    shellSum+=0.5*(U1.re*U1.re+U1.im*U1.im+U2.re*U2.re+U2.im*U2.im+U3.re*U3.re+U3.im*U3.im);
                    cnt++;
                }
            }
        }
        if(cnt==0) continue;
        const real kk=(real)ks;
        const real Enum=sc*sc*shellSum*(4.0*PI*kk*kk/(real)cnt);
        std::printf("[HIT-IC]  %2d  %5d  %.6e  %.6e\n",ks,cnt,(double)Enum,(double)Ek(kk));
    }
}

// 多波剖面（Peng 2021 式(31)；x∈[−1,1]；供 1D 熵波载体算例使用）
inline real multiwaveProfile(real x)
{
    const real delta=0.005, alpha=10.0, a=0.5, z=-0.7;
    const real beta=std::log(2.0)/(36.0*delta*delta);
    auto G=[&](real xx){ const real d=xx-z; return std::exp(-beta*d*d); };
    auto F=[&](real xx){ const real t=1.0-alpha*alpha*(xx-a)*(xx-a); return t>0.0? std::sqrt(t):0.0; };
    if(x>=-0.8 && x<-0.6) return (G(x-delta)+G(x+delta)+4.0*G(x))/6.0;
    if(x>=-0.4 && x<-0.2) return 1.0;
    if(x>= 0.0 && x< 0.2) return 1.0-std::abs(10.0*(x-0.1));
    if(x>= 0.4 && x< 0.6) return (F(x-delta)+F(x+delta)+4.0*F(x))/6.0;
    return 0.0;
}
}

/**
 * @brief 解初始化
 * @param grid 网格块指针
 * @param sol 解数据指针
 * 
 * 根据不同的方程类型和算例编号初始化解数据。
 */
void Initializer::solInit(Block* grid,Data* sol)
{
    std::vector<real> tempsol;
    switch (info->eqType)
    {
        /*case begin*/
    case LINEARCONV1D:
        if (grid->dim!=1)
        { std::cout<<"initialize: dim error \n";return;}
        switch (info->nCase)
        {
        case 0:
            tempsol.reserve(grid->icMax[0]);
            for(int i=0;i<grid->icMax[0];i++)
            {
                real x=(*grid)(i,0);
                tempsol.push_back(1+sin(x));
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;
        case 1:
            //ADR
            tempsol.reserve(grid->icMax[0]);
            for(int i=0;i<grid->icMax[0];i++)
            {
                real x=(*grid)(i,0);
                tempsol.push_back(1+sin(x));
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;
        
        default:
            break;
        }
        break;
    /*case end*/
    
    /*case begin*/
    case BURGERS1D:
        if (grid->dim!=1)
        { std::cout<<"initialize: dim error \n";return;}
        switch (info->nCase)
        {
        /*case 0 begin*/
        case 0:
            tempsol.reserve(grid->icMax[0]);
            for(int i=0;i<grid->icMax[0];i++)
            {
                real x=(*grid)(i,0);
                tempsol.push_back(-sin(M_PI*x));
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;
        /*case 0 end*/
        
        default:
            break;
        }
        break;
    /*case end*/

    case EULER:
        /*case begin*/
        if (grid->dim==1)
        switch (info->nCase)
        {
        /*case 0 begin*/
        case 0:
            tempsol.reserve(grid->icMax[0]);
            for(int i=0;i<grid->icMax[0];i++)
            {
                real x=(*grid)(i,0);
                real gamma=1.4;
                if (x<0)
                {
                    tempsol.push_back(1);
                    tempsol.push_back(0);
                    tempsol.push_back(1.0/(gamma-1)*1);
                }
                else
                {
                    tempsol.push_back(0.125);
                    tempsol.push_back(0);
                    tempsol.push_back(1.0/(gamma-1)*0.1);
                }
                
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;

            case 1:
            //shu-osher
            tempsol.reserve(grid->icMax[0]);
            for(int i=0;i<grid->icMax[0];i++)
            {
                real x=(*grid)(i,0);
                real gamma=GAMMA;
                if (x<=0.5)
                {
                    tempsol.push_back(3.857143);
                    tempsol.push_back(3.857143*2.629369);
                    tempsol.push_back(1.0/(gamma-1.0)*(10.0+1.0/3.0)
                                     +3.857143*2.629369*2.629369/2);
                }
                else
                {
                    tempsol.push_back(1.0+0.2*sin(5.0*x));
                    tempsol.push_back(0.0);
                    tempsol.push_back(1.0/(gamma-1.0)*1.0);
                }
                
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;

            case 2:
            //Lax
            tempsol.reserve(grid->icMax[0]);
            for(int i=0;i<grid->icMax[0];i++)
            {
                real x=(*grid)(i,0);
                real gamma=GAMMA;
                real r,u,p;
                if (x<=0)
                {
                    r=0.445;u=0.698;p=3.528;
                }
                else
                {
                    r=0.5;u=0;p=0.571;
                }
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(1.0/(gamma-1)*p+r*u*u/2);
                
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;

            case 3:
            //sedov
            tempsol.reserve(grid->icMax[0]);
            for(int i=0;i<grid->icMax[0];i++)
            {
                real x=(*grid)(i,0);
                
                real gamma=GAMMA;
                real r,u,E;
                if (std::abs(x)<1e-10)
                {
                    real dx=(*grid)(i,0)-(*grid)(i-1,0);
                    r=1;u=0;E=3200000.0/dx;
                }
                else
                {
                    r=1;u=0;E=1e-12;
                }
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(r*E);
                
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;

            case 4:
            //Woodward-Colella
            tempsol.reserve(grid->icMax[0]);
            for(int i=0;i<grid->icMax[0];i++)
            {
                real x=(*grid)(i,0);
                real gamma=GAMMA;
                real r,u,p;
                if (x<0.1)
                {
                    r=1;u=-2;p=0.4;
                }
                else if (x>=0.9)
                {
                    r=1;u=0;p=100;
                }
                else
                {
                    r=1;u=0;p=0.01;
                }
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(1.0/(gamma-1)*p+r*u*u/2);
                
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;

        case 5:
            //双稀疏波
            tempsol.reserve(grid->icMax[0]);
            for(int i=0;i<grid->icMax[0];i++)
            {
                real x=(*grid)(i,0);
                real gamma=GAMMA;
                real r,u,p;
                if (x<0)
                {
                    r=1;u=-2;p=0.4;
                }
                else
                {
                    r=1;u=2;p=0.4;
                }
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(1.0/(gamma-1)*p+r*u*u/2);
                
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;
        case 6: // 多波熵波：ρ 剖面精确平流（Peng 2021 式(31)；u=1、p 均匀，周期域）
        {
            tempsol.reserve(grid->icMax[0]);
            for(int i=0;i<grid->icMax[0];i++)
            {
                real x=(*grid)(i,0);
                real gamma=GAMMA;
                real r=1.0+multiwaveProfile(x);
                real u=1.0, p=1.0;
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(1.0/(gamma-1.0)*p + r*u*u/2.0);
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;
        }
        case 7: // 多波熵波·缩放变体（扰动 ×1e-3；标度不变性检验用）
        {
            tempsol.reserve(grid->icMax[0]);
            for(int i=0;i<grid->icMax[0];i++)
            {
                real x=(*grid)(i,0);
                real gamma=GAMMA;
                real r=1.0+1.0e-3*multiwaveProfile(x);
                real u=1.0, p=1.0;
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(1.0/(gamma-1.0)*p + r*u*u/2.0);
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;
        }
        case 8: // 双爆轰（Peng 2021 式(35)；域 [0,1]）
        {
            tempsol.reserve(grid->icMax[0]);
            for(int i=0;i<grid->icMax[0];i++)
            {
                real x=(*grid)(i,0);
                real gamma=GAMMA;
                real r=1.0, u=0.0, p;
                if(x<0.1) p=1000.0;
                else if(x<0.9) p=0.01;
                else p=100.0;
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(1.0/(gamma-1.0)*p + r*u*u/2.0);
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;
        }
        /*case 0 end*/
        default:
            break;
        }
        else if (grid->dim==2)
        switch (info->nCase)
        {
        case 0:
            /**
             * @brief 二维黎曼问题 - Configuration 3 (Schulz-Rinne et al., 2004)
             * 
             * 特征: 四象限不同状态，产生双激波和接触间断结构
             * 计算域: x∈[-0.5, 0.5], y∈[-0.5, 0.5]
             * 初始条件分界线: x=0.3, y=0.3
             * 
             * 四象限状态 (ρ, u, v, p):
             * Q1 (右上): (1.5, 0, 0, 1.5)
             * Q2 (左上): (0.5323, 1.206, 0, 0.3)
             * Q3 (左下): (0.138, 1.206, 1.206, 0.029)
             * Q4 (右下): (0.5323, 0, 1.206, 0.3)
             */
            //2D Riemann Problem - Configuration 3;
            tempsol.reserve(grid->icMax[0]*grid->icMax[1]);
            for(int i=0;i<grid->icMax[1];i++)
            for(int j=0;j<grid->icMax[0];j++)
            {
                real x=(*grid)(i*grid->icMax[0]+j,0);
                real y=(*grid)(i*grid->icMax[0]+j,1);
                real gamma=GAMMA;
                if (x>0.3)
                {
                    if (y>0.3)
                    {
                        tempsol.push_back(1.5);
                        tempsol.push_back(0);
                        tempsol.push_back(0);
                        tempsol.push_back(1.0/(gamma-1)*1.5);
                    }
                    else
                    {
                        tempsol.push_back(0.5323);
                        tempsol.push_back(0);
                        tempsol.push_back(0.5323*1.206);
                        tempsol.push_back(1.0/(gamma-1)*0.3+1.206*1.206/2*0.5323);
                    }
                }
                else
                {
                    if (y>0.3)
                    {
                        tempsol.push_back(0.5323);
                        tempsol.push_back(0.5323*1.206);
                        tempsol.push_back(0);
                        tempsol.push_back(1.0/(gamma-1)*0.3+1.206*1.206/2*0.5323);
                    }
                    else
                    {
                        tempsol.push_back(0.138);
                        tempsol.push_back(0.138*1.206);
                        tempsol.push_back(0.138*1.206);
                        tempsol.push_back(1.0/(gamma-1)*0.029+1.206*1.206*0.138);
                    }
                }
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
                else std::cout<<"initialize: length error \n";
            break;
        case 1:
            /**
             * @brief 二维黎曼问题 - Configuration 4 (Vortex) (Schulz-Rinne et al., 2004)
             * 
             * 特征: 产生涡旋结构的黎曼问题配置
             * 计算域: x∈[-0.5, 0.5], y∈[-0.5, 0.5]
             * 初始条件分界线: x=0, y=0 (坐标轴)
             * 
             * 四象限状态 (ρ, u, v, p):
             * Q1 (右上): (1.0, 0.75, -0.5, 1.0)
             * Q2 (左上): (2.0, 0.75, 0.5, 1.0)
             * Q3 (左下): (1.0, -0.75, 0.5, 1.0)
             * Q4 (右下): (3.0, -0.75, -0.5, 1.0)
             * 
             * 注意: 所有象限压力相同(p=1.0)，速度场形成涡旋结构
             */
            //2D Riemann Problem - Configuration 4 (Vortex);
            //x [-0.5,0.5]
            //y [-0.5,0.5]
            tempsol.reserve(grid->icMax[0]*grid->icMax[1]);
            for(int i=0;i<grid->icMax[1];i++)
            for(int j=0;j<grid->icMax[0];j++)
            {
                real x=(*grid)(i*grid->icMax[0]+j,0);
                real y=(*grid)(i*grid->icMax[0]+j,1);
                real gamma=GAMMA;
                real r,u,v,p;
                if (x>0.0)
                {
                    if (y>0.0)
                    {r=1.0;u=0.75;v=-0.5;p=1.0;}
                    else
                    {r=3.0;u=-0.75;v=-0.5;p=1.0;}
                }
                else
                {
                    if (y>0)
                    {r=2.0;u=0.75;v=0.5;p=1.0;}
                    else
                    {r=1.0;u=-0.75;v=0.5;p=1.0;}
                }
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(r*v);
                tempsol.push_back(1.0/(gamma-1)*p+r*(u*u+v*v)/2);
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
                else std::cout<<"initialize: length error \n";
            break;

        case 2:
            //Implosion problem;
            //x [-0.3,0.3]
            //y [-0.3,0.3]
            tempsol.reserve(grid->icMax[0]*grid->icMax[1]);
            for(int i=0;i<grid->icMax[1];i++)
            for(int j=0;j<grid->icMax[0];j++)
            {
                real eps=1e-10;
                real x=(*grid)(i*grid->icMax[0]+j,0);
                real y=(*grid)(i*grid->icMax[0]+j,1);
                real gamma=GAMMA;
                real r,u,v,p;
                if (std::abs(x)+std::abs(y)<0.15-eps)
                {
                    r=0.125;
                    u=0.0;
                    v=0.0;
                    p=0.14;
                }
                else
                {
                    r=1.0;
                    u=0.0;
                    v=0.0;
                    p=1.0;
                }
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(r*v);
                tempsol.push_back(1.0/(gamma-1)*p+r*(u*u+v*v)/2);
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
                else std::cout<<"initialize: length error \n";
            break;

            
        case 3:
            //R-T instability;
            tempsol.reserve(grid->icMax[0]*grid->icMax[1]);
            for(int i=0;i<grid->icMax[1];i++)
            for(int j=0;j<grid->icMax[0];j++)
            {
                real x=(*grid)(i*grid->icMax[0]+j,0);
                real y=(*grid)(i*grid->icMax[0]+j,1);
                real gamma=GAMMA;
                real r,u,v,p;
                if(y<=0.5){
                    r=2.0;u=0;p=2.0*y+1.0;
                    real c=std::sqrt(GAMMA*p/r);
                    v=-0.025*c*cos(8.0*M_PI*(x<0.125?x:(0.25-x)));
                }
                else{
                    r=1.0;u=0;p=y+3.0/2.0;
                    real c=std::sqrt(GAMMA*p/r);
                    v=-0.025*c*cos(8.0*M_PI*(x<0.125?x:(0.25-x)));
                }
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(r*v);
                tempsol.push_back(1.0/(gamma-1)*p+r*(u*u+v*v)/2);
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
                else std::cout<<"initialize: length error \n";
            break;

        case 4:
            //Double Mach reflection;
            tempsol.reserve(grid->icMax[0]*grid->icMax[1]);
            for(int i=0;i<grid->icMax[1];i++)
            for(int j=0;j<grid->icMax[0];j++)
            {
                real x=(*grid)(i*grid->icMax[0]+j,0);
                real y=(*grid)(i*grid->icMax[0]+j,1);
                real gamma=GAMMA;
                real r,u,v,p;
                if(y>=std::sqrt(3)*(x-1.0/6.0)){
                    r=8.0;u=8.25*cos(M_PI/6);v=-8.25*sin(M_PI/6);p=116.5;
                }
                else{
                    r=1.4;u=0;v=0;p=1.0;
                }
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(r*v);
                tempsol.push_back(1.0/(gamma-1)*p+r*(u*u+v*v)/2);
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
                else std::cout<<"initialize: length error \n";
            break;

        case 5:
            /**
             * @brief 二维黎曼问题 - Configuration 12 (Schulz-Rinne et al., 2004)
             * 
             * 特征: 混合波系结构，包含激波、稀疏波和接触间断
             * 计算域: x∈[-0.5, 0.5], y∈[-0.5, 0.5]
             * 初始条件分界线: x=0, y=0 (坐标轴)
             * 
             * 四象限状态 (ρ, u, v, p):
             * Q1 (右上): (0.5313, 0, 0, 0.4)
             * Q2 (左上): (1.0, 0.7276, 0, 1.0)
             * Q3 (左下): (0.8, 0, 0, 1.0)
             * Q4 (右下): (1.0, 0, 0.7276, 1.0)
             */
            //2D Riemann Problem - Configuration 12;
            //x [-0.5,0.5]
            //y [-0.5,0.5]
            tempsol.reserve(grid->icMax[0]*grid->icMax[1]);
            for(int i=0;i<grid->icMax[1];i++)
            for(int j=0;j<grid->icMax[0];j++)
            {
                real x=(*grid)(i*grid->icMax[0]+j,0);
                real y=(*grid)(i*grid->icMax[0]+j,1);
                real gamma=GAMMA;
                real r,u,v,p;
                if (x>0.0)
                {
                    if (y>0.0)
                    {r=0.5313;u=0.0;v=0.0;p=0.4;}
                    else
                    {r=1.0;u=0.0;v=0.7276;p=1.0;}
                }
                else
                {
                    if (y>0)
                    {r=1.0;u=0.7276;v=0.0;p=1.0;}
                    else
                    {r=0.8;u=0.0;v=0.0;p=1.0;}
                }
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(r*v);
                tempsol.push_back(1.0/(gamma-1)*p+r*(u*u+v*v)/2);
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
                else std::cout<<"initialize: length error \n";
            break;
        
        case 6:
            /**
             * @brief 二维黎曼问题 - Configuration 3 变体 (Schulz-Rinne et al., 2004)
             * 
             * 特征: 与 Configuration 3 类似但参数略有不同的配置
             * 计算域: x∈[-0.5, 0.5], y∈[-0.5, 0.5]
             * 初始条件分界线: x=0, y=0 (坐标轴)
             * 
             * 四象限状态 (ρ, u, v, p):
             * Q1 (右上): (0.5323, 0, 0, 0.4)
             * Q2 (左上): (1.0, 0.7276, 0, 1.0)
             * Q3 (左下): (1.0, 0, 0.7276, 1.0)
             * Q4 (右下): (0.8, 0, 0, 1.0)
             * 
             * 注意: 与 case 5 (Config 12) 的区别在于 Q1 和 Q3 的状态互换
             */
            //2D Riemann Problem - Configuration 3 variant;
            //x [-0.5,0.5]
            //y [-0.5,0.5]
            tempsol.reserve(grid->icMax[0]*grid->icMax[1]);
            for(int i=0;i<grid->icMax[1];i++)
            for(int j=0;j<grid->icMax[0];j++)
            {
                real x=(*grid)(i*grid->icMax[0]+j,0);
                real y=(*grid)(i*grid->icMax[0]+j,1);
                real gamma=GAMMA;
                real r,u,v,p;
                if (x>0.0)
                {
                    if (y>0.0)
                    {r=0.5323;u=0.0;v=0.0;p=0.4;}
                    else
                    {r=0.8;u=0.0;v=0.0;p=1.0;}
                }
                else
                {
                    if (y>0)
                    {r=1.0;u=0.7276;v=0.0;p=1.0;}
                    else
                    {r=1.0;u=0.0;v=0.7276;p=1.0;}
                }
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(r*v);
                tempsol.push_back(1.0/(gamma-1)*p+r*(u*u+v*v)/2);
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
                else std::cout<<"initialize: length error \n";
            break;
        
        case 7:
            /**
             * @brief 二维黎曼问题 - Configuration 6 (Schulz-Rinne et al., 2004)
             * 
             * 特征: 反向涡旋结构，与 Configuration 4 速度场方向相反
             * 计算域: x∈[-0.5, 0.5], y∈[-0.5, 0.5]
             * 初始条件分界线: x=0, y=0 (坐标轴)
             * 
             * 四象限状态 (ρ, u, v, p):
             * Q1 (右上): (1.0, -0.75, -0.5, 1.0)
             * Q2 (左上): (2.0, -0.75, 0.5, 1.0)
             * Q3 (左下): (1.0, 0.75, 0.5, 1.0)
             * Q4 (右下): (3.0, 0.75, -0.5, 1.0)
             * 
             * 注意: 与 case 1 (Config 4) 的区别在于水平速度 u 的符号相反
             *       Config 4: Q1,Q2 的 u>0; Q3,Q4 的 u<0
             *       Config 6: Q1,Q2 的 u<0; Q3,Q4 的 u>0
             */
            //2D Riemann Problem - Configuration 6 (Reverse Vortex);
            //x [-0.5,0.5]
            //y [-0.5,0.5]
            tempsol.reserve(grid->icMax[0]*grid->icMax[1]);
            for(int i=0;i<grid->icMax[1];i++)
            for(int j=0;j<grid->icMax[0];j++)
            {
                real x=(*grid)(i*grid->icMax[0]+j,0);
                real y=(*grid)(i*grid->icMax[0]+j,1);
                real gamma=GAMMA;
                real r,u,v,p;
                if (x>0.0)
                {
                    if (y>0.0)
                    {r=1.0;u=-0.75;v=-0.5;p=1.0;}
                    else
                    {r=3.0;u= 0.75;v=-0.5;p=1.0;}
                }
                else
                {
                    if (y>0)
                    {r=2.0;u=-0.75;v= 0.5;p=1.0;}
                    else
                    {r=1.0;u= 0.75;v= 0.5;p=1.0;}
                }
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(r*v);
                tempsol.push_back(1.0/(gamma-1)*p+r*(u*u+v*v)/2);
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
                else std::cout<<"initialize: length error \n";
            break;

        case 8: // 2D TGV（周期；任务 8.6 新增——code-to-code 对照的二维载体）
            tempsol.reserve(grid->icMax[0]*grid->icMax[1]);
            for(int j=0;j<grid->icMax[1];j++)
            for(int i=0;i<grid->icMax[0];i++)
            {
                int idx=i+j*grid->icMax[0];
                real x=(*grid)(idx,0), y=(*grid)(idx,1);
                real gamma=GAMMA;
                real r=1.0;
                real u=std::sin(x)*std::cos(y);
                real v=-std::cos(x)*std::sin(y);
                real p=1.0+(1.0/8.0)*(std::cos(2.0*x)+std::cos(2.0*y));
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(r*v);
                tempsol.push_back(1.0/(gamma-1.0)*p + r*(u*u+v*v)/2.0);
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;
        case 9: // 激波/剪切层相互作用（Peng 2019 §4.2.3 式(72)-(74)；域 [0,200]x[-20,20]；下支 ρ=0.3626，Tang 2024 转录作 0.3636，从原文献）
            tempsol.reserve(grid->icMax[0]*grid->icMax[1]);
            for(int j=0;j<grid->icMax[1];j++)
            for(int i=0;i<grid->icMax[0];i++)
            {
                int idx=i+j*grid->icMax[0];
                real y=(*grid)(idx,1);
                real gamma=GAMMA;
                real r=(y>=0.0)? 1.6374 : 0.3626;
                real u=2.5+0.5*std::tanh(2.0*y);
                real v=0.0;
                real p=0.3327;
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(r*v);
                tempsol.push_back(1.0/(gamma-1.0)*p + r*(u*u+v*v)/2.0);
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;
        default:
            break;
        }
        else if (grid->dim==3)
        switch (info->nCase)
        {
        case 0: // Taylor-Green 涡（无粘档，Fu et al. 2018 式 (40)；任务 6 新增）
            tempsol.reserve(grid->icMax[0]*grid->icMax[1]*grid->icMax[2]);
            for(int k=0;k<grid->icMax[2];k++)
            for(int j=0;j<grid->icMax[1];j++)
            for(int i=0;i<grid->icMax[0];i++)
            {
                int idx=i+j*grid->icMax[0]+k*grid->icMax[0]*grid->icMax[1];
                real x=(*grid)(idx,0), y=(*grid)(idx,1), z=(*grid)(idx,2);
                real gamma=GAMMA;
                real r=1.0;
                real u=std::sin(x)*std::cos(y)*std::cos(z);
                real v=-std::cos(x)*std::sin(y)*std::cos(z);
                real w=0.0;
                real p=100.0+(1.0/16.0)*((std::cos(2*x)+std::cos(2*y))*(2.0+std::cos(2*z))-2.0);
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(r*v);
                tempsol.push_back(r*w);
                tempsol.push_back(1.0/(gamma-1.0)*p + r*(u*u+v*v+w*w)/2.0);
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;
        case 1: // 均匀场（三维核验用；任务 6 新增）
            tempsol.reserve(grid->icMax[0]*grid->icMax[1]*grid->icMax[2]);
            for(int k=0;k<grid->icMax[2];k++)
            for(int j=0;j<grid->icMax[1];j++)
            for(int i=0;i<grid->icMax[0];i++)
            {
                real gamma=GAMMA;
                real r=1.0, u=1.0, v=2.0, w=3.0, p=100.0;
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(r*v);
                tempsol.push_back(r*w);
                tempsol.push_back(1.0/(gamma-1.0)*p + r*(u*u+v*v+w*w)/2.0);
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;
        case 2: // 三维等熵涡对流（解析解基准；涡轴与平流沿体对角，三方向链路全激活）
            tempsol.reserve(grid->icMax[0]*grid->icMax[1]*grid->icMax[2]);
            for(int k=0;k<grid->icMax[2];k++)
            for(int j=0;j<grid->icMax[1];j++)
            for(int i=0;i<grid->icMax[0];i++)
            {
                int idx=i+j*grid->icMax[0]+k*grid->icMax[0]*grid->icMax[1];
                real x=(*grid)(idx,0), y=(*grid)(idx,1), z=(*grid)(idx,2);
                real gamma=GAMMA;
                real pi=3.14159265358979323846;
                real epsV=5.0;                                                  // 涡强度参数（2D 标准测试常用值）
                real n3=1.0/std::sqrt(3.0);                                     // 体对角单位向量分量
                real Ux=n3, Uy=n3, Uz=n3;                                       // 平流速度（|U|=1，与涡轴同向）
                real dx=x-pi, dy=y-pi, dz=z-pi;                                 // 相对涡心（域中心）坐标
                real dn=(dx+dy+dz)*n3;                                          // d·n：沿涡轴分量
                real r2=dx*dx+dy*dy+dz*dz-dn*dn;                                // 到涡轴距离平方
                real g=std::exp((1.0-r2)*0.5);                                  // 高斯核 e^{(1-r^2)/2}
                real amp=epsV/(2.0*pi)*g;
                real u=Ux+amp*n3*(dz-dy);                                       // 周向速度扰动 = (eps/2pi) g (n × d)
                real v=Uy+amp*n3*(dx-dz);
                real w=Uz+amp*n3*(dy-dx);
                // 等熵涡平衡：T 修正 + 等熵关系给出 rho、p（解析平移解 q(x,t)=q0(x-Ut)）
                real T=1.0-(gamma-1.0)*epsV*epsV/(8.0*gamma*pi*pi)*std::exp(1.0-r2);
                real r=std::pow(T,1.0/(gamma-1.0));
                real p=std::pow(T,gamma/(gamma-1.0));
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(r*v);
                tempsol.push_back(r*w);
                tempsol.push_back(1.0/(gamma-1.0)*p + r*(u*u+v*v+w*w)/2.0);
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;
        case 3: // 三维熵波（密度正弦沿斜向匀速平流；欧拉方程严格精确解，且波矢取整数分量保证周期相容）
            tempsol.reserve(grid->icMax[0]*grid->icMax[1]*grid->icMax[2]);
            for(int k=0;k<grid->icMax[2];k++)
            for(int j=0;j<grid->icMax[1];j++)
            for(int i=0;i<grid->icMax[0];i++)
            {
                int idx=i+j*grid->icMax[0]+k*grid->icMax[0]*grid->icMax[1];
                real x=(*grid)(idx,0), y=(*grid)(idx,1), z=(*grid)(idx,2);
                real gamma=GAMMA;
                real u=1.0/std::sqrt(3.0);                                      // 平流速度 U=(1,1,1)/sqrt(3)
                real ph=x+y+z;                                                  // 波矢 n=(1,1,1)，整数分量 → 2pi 域上严格周期
                real r=1.0+0.2*std::sin(ph);                                    // 密度正弦扰动
                real v=u, w=u, p=1.0;
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(r*v);
                tempsol.push_back(r*w);
                tempsol.push_back(1.0/(gamma-1.0)*p + r*(u*u+v*v+w*w)/2.0);
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;
        case 4: // 2D TGV 的 z 均匀嵌入（任务 8.6：code-to-code 对照的三维侧；初值与 2D nCase=8 逐式相同）
            tempsol.reserve(grid->icMax[0]*grid->icMax[1]*grid->icMax[2]);
            for(int k=0;k<grid->icMax[2];k++)
            for(int j=0;j<grid->icMax[1];j++)
            for(int i=0;i<grid->icMax[0];i++)
            {
                int idx=i+j*grid->icMax[0]+k*grid->icMax[0]*grid->icMax[1];
                real x=(*grid)(idx,0), y=(*grid)(idx,1);
                real gamma=GAMMA;
                real r=1.0;
                real u=std::sin(x)*std::cos(y);                                 // 与 2D nCase=8 相同（z 均匀）
                real v=-std::cos(x)*std::sin(y);
                real w=0.0;
                real p=1.0+(1.0/8.0)*(std::cos(2.0*x)+std::cos(2.0*y));
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(r*v);
                tempsol.push_back(r*w);
                tempsol.push_back(1.0/(gamma-1.0)*p + r*(u*u+v*v+w*w)/2.0);
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;
        case 5: // 可压缩均匀各向同性湍流衰减 HIT（Johnsen et al. 2010 附录 A；方案见 case_for_comment8/HIT）
        {
            const int nxc=grid->icMax[0], nyc=grid->icMax[1], nzc=grid->icMax[2];
            const size_t nc=(size_t)nxc*nyc*nzc;
            std::vector<real> velU,velV,velW;
            hitInitialVelocity(velU,velV,velW,nxc,nyc,nzc,4.0,1.0,20260920ull);   // k0=4、u_rms=1、固定种子
            if(velU.size()!=nc)
            {
                std::cout<<"HIT initialize: velocity field generation failed\n";
                break;
            }
            // ρ=1；T=p/ρ=3/(0.6²γ)=5.952380952…（实现 Ma_t,0=√3·u_rms/⟨c⟩=0.6，核定见定义档 §5）
            const real T0=3.0/(0.6*0.6*GAMMA);
            tempsol.reserve(nc);
            for(size_t idx=0;idx<nc;idx++)
            {
                const real r=1.0;
                const real u=velU[idx], v=velV[idx], w=velW[idx];
                const real p=r*T0;
                tempsol.push_back(r);
                tempsol.push_back(r*u);
                tempsol.push_back(r*v);
                tempsol.push_back(r*w);
                tempsol.push_back(1.0/(GAMMA-1.0)*p + r*(u*u+v*v+w*w)/2.0);
            }
            if (tempsol.size()==sol->size()) sol->setValue(tempsol);
            else std::cout<<"initialize: length error \n";
            break;
        }
        default:
            break;
        }

    /*case end*/
    default:
        break;
    }
}


void Initializer::initUniformBlock(Block* block)
{
    if(!block)
    {
        std::cout<<"Initializer error: empty shared_ptr Block\n";
        return;
    }
    auto iMax=info->iMax;
    int dim=info->dim;
    auto icMax=info->icMax();
    int nVer=1,nCel=1;
    for (int i = 0; i < 3; i++)
    {
        nVer*=iMax[i];
        nCel*=icMax[i];
    }

    block->dim=dim;
    block->nVer=nVer;
    block->nCel=nCel;
    block->icMax=icMax;
    block->iMax=iMax;
    block->coorVer.init(nVer,dim);
    block->coorCel.init(nCel,dim);
    block->intervalCel.init(nCel,dim);
    block->inited=true;

    //consth related
    info->interval=0;
    std::array<real,3> inters;
    double cmin=info->calZone[0];
    double cmax=info->calZone[1];
    // double interval=round((cmax-cmin)/(iMax[0]-1)*1e8)/1e8;
    double interval=(cmax-cmin)/(iMax[0]-1);
    info->interval=interval;
    if (std::abs(*std::max_element(inters.begin(),inters.end())-*std::min_element(inters.begin(),inters.end()))>1e-10)
    std::cout<<"Initalize error: interval incorrect\n";
    info->constH=true;

    //for vertex
    for (int idim = 0; idim < dim; idim++)
    {
        double cmin=info->calZone[idim*2];
        double cmax=info->calZone[idim*2+1];
        double interval=info->interval;

        int l,m,n;
        int* onedIndex=((idim == 0 ) ? &l : ( idim == 1 ? &m : &n));
        
        for (l = 0; l < iMax[0]; l++)
        for (m = 0; m < iMax[1]; m++)
        for (n = 0; n < iMax[2]; n++)
        {
            int globalIndex=l+m*iMax[0]+n*iMax[0]*iMax[1];
            block->coorVer(globalIndex,idim)=(cmin+(*onedIndex)*interval);
            //block->coorVer(globalIndex,idim)=cmin+(*onedIndex)*interval;
        }
    }
    
    //for cellcenter
    // for (int idim = 0; idim < dim; idim++)
    // {
    //     int l,m,n,iLen=(dim==1?2:dim==2? 4:8);
    //     std::vector<int> index;
    //     index.resize(iLen);
    //     for (l = 0; l < icMax[0]; l++)
    //     for (m = 0; m < icMax[1]; m++)
    //     for (n = 0; n < icMax[2]; n++)
    //     {
    //         int iVerGlobal=l+m*iMax[0]+n*iMax[0]*iMax[1];
    //         int iCelGlobal=l+m*icMax[0]+n*icMax[0]*icMax[1];
    //         double temp=0;
    //         index[0]=iVerGlobal;
    //         index[1]=iVerGlobal+1;
    //         if (dim>=2)
    //         {
    //             index[2]=iVerGlobal+iMax[0];
    //             index[3]=iVerGlobal+iMax[0]+1;
    //         }
    //         if (dim>=3)
    //         {
    //             index[4]=iVerGlobal+iMax[0]*iMax[1];
    //             index[5]=iVerGlobal+1+iMax[0]*iMax[1];
    //             index[6]=iVerGlobal+iMax[0]+iMax[0]*iMax[1];
    //             index[7]=iVerGlobal+iMax[0]+1+iMax[0]*iMax[1];
    //         }
    //         real tempInterval;
    //         for(auto iver:index) temp+=block->coorVer(iver,idim);
    //         block->coorCel(iCelGlobal,idim)=temp/iLen;
    //         //only lower line approx_1imate
    //         if(idim==0) block->intervalCel(iCelGlobal,idim)=block->coorVer(index[1],idim)-block->coorVer(index[0],idim);
    //         if(idim==1) block->intervalCel(iCelGlobal,idim)=block->coorVer(index[2],idim)-block->coorVer(index[0],idim);
    //         if(idim==2) block->intervalCel(iCelGlobal,idim)=block->coorVer(index[4],idim)-block->coorVer(index[0],idim);
    //     }
    // }

    for (int idim = 0; idim < dim; idim++)
    {
        double cmin=info->calZone[idim*2];      // 按维原点（缺失时各维沿用 x 原点，原点不同的域单元中心坐标错位）
        double cmax=info->calZone[idim*2+1];
        int l,m,n,iLen=(dim==1?2:dim==2? 4:8);
        std::vector<int> index;
        index.resize(iLen);
        for (l = 0; l < icMax[0]; l++)
        for (m = 0; m < icMax[1]; m++)
        for (n = 0; n < icMax[2]; n++)
        {
            int iCelGlobal=l+m*icMax[0]+n*icMax[0]*icMax[1];
            int* onedIndex=((idim == 0 ) ? &l : ( idim == 1 ? &m : &n));
            block->coorCel(iCelGlobal,idim)=(cmin+((*onedIndex)+0.5)*interval);
            //only lower line approx_1imate
            if(idim==0) block->intervalCel(iCelGlobal,idim)=interval;
            if(idim==1) block->intervalCel(iCelGlobal,idim)=interval;
            if(idim==2) block->intervalCel(iCelGlobal,idim)=interval;
        }
    }
}


Initializer::Initializer()
{

}

Initializer::Initializer(Info* info_)
{
    info=info_;
}


void Initializer::initEqution(Equation* eq,Block* block)
{
    if(!eq||!block)
    {
        std::cout<<"Initializer error: empty shared_ptr Equation or Block\n";
        return;
    }
    eq->n=block->icMax[0]*block->icMax[1]*block->icMax[2];
    eq->nPrim=info->nPrim();
    eq->nCons=info->nCons();
    eq->type=info->eqType;
    eq->dim=block->dim;

    
    eq->rhs=new Data(eq->n,eq->nCons);
    eq->cons=new Data(eq->n,eq->nCons);
    eq->prim=new Data(eq->n,eq->nPrim);
    eq->inited=true;
    eq->cons->setvarName(info->getVarNameListCons());
    eq->prim->setvarName(info->getVarNameListPrim());
    eq->rhs->setvarName(info->getVarNameListRhs());

    solInit(block,eq->cons);
}


void Initializer::initBnds(Bnds* bnds,Equation* eqn,std::array<int,3> iMax,Block* block)
{
    bnds->iMax=iMax;
    bnds->dim=info->dim;

    int nGhost=info->nGhostCell();
    int nCons=info->nCons();
    int nPrim=info->nPrim();

    std::array<int,3> nBnds  {bnds->iMax[1]*bnds->iMax[2],
                            bnds->iMax[0]*bnds->iMax[2],
                            bnds->iMax[0]*bnds->iMax[1]};
    int nBnd=0;
    for (int i = 0; i < bnds->dim; i++)
    {
        nBnd+=nBnds[i]*2;
    }
    
    

    bnds->oneDBnds.resize(nBnd);
    std::array<int,2> offsets;
    switch (info->eqType)
    {
    case LINEARCONV1D:
    case BURGERS1D:
        {
            for (int i = 0; i < 2; i++)
            {
                bnds->oneDBnds.at(i)=std::make_shared<OneDBnd>(nGhost,nPrim,PERIODIC1D);
            }
            offsets=calOffsetInverse(1,0,0,bnds->iMax);
            bnds->oneDBnds.at(0)->setUpdate(eqn->prim,offsets[0],offsets[1]);

            offsets=calOffset(1,0,0,bnds->iMax);
            bnds->oneDBnds.at(1)->setUpdate(eqn->prim,offsets[0],offsets[1]);
        }
        break;
    case EULER:
        if(eqn->dim==1)
        {
            BndType Xtype=SUPERSONICOUTLET;
            if(info->nCase==4) {Xtype=SYMMETRY1D;}
            if(info->nCase==6||info->nCase==7) {Xtype=PERIODIC1D;}   // 多波熵波：1D 周期
            bnds->oneDBnds.at(0)=std::make_shared<OneDBnd>(nGhost,nPrim,Xtype);
            // 周期边界的拷贝方向与外推/对称相反：左边界取数组尾端、右边界取数组首端
            if(Xtype==PERIODIC1D) offsets=calOffsetInverse(1,0,0,bnds->iMax);
            else                  offsets=calOffset(1,0,0,bnds->iMax);
            bnds->oneDBnds.at(0)->setUpdate(eqn->prim,offsets[0],offsets[1]);


            bnds->oneDBnds.at(1)=std::make_shared<OneDBnd>(nGhost,nPrim,Xtype);
            if(Xtype==PERIODIC1D) offsets=calOffset(1,0,0,bnds->iMax);
            else                  offsets=calOffsetInverse(1,0,0,bnds->iMax);
            bnds->oneDBnds.at(1)->setUpdate(eqn->prim,offsets[0],offsets[1]);
        }
        else if (eqn->dim==2)
        {
            BndType Xtype=SUPERSONICOUTLET,Ytype=SUPERSONICOUTLET;
            if(info->nCase==2) {Xtype=SYMMETRYX;Ytype=SYMMETRYY;}
            if(info->nCase==3) {Xtype=SYMMETRYX;Ytype=DIRICLET;}
            if(info->nCase==4)
            {
                initDoubleMachBnds(bnds,eqn,iMax,block);
                break;
            }
            if(info->nCase==9)
            {
                initShockShearBnds(bnds,eqn,iMax,block);
                break;
            }
            if(info->nCase==8)   // 2D TGV：四面全周期（任务 8.6 新增；与 3D 六面周期同型、二维化）
            {
                for (int i = 0; i < iMax[1]; i++)
                {
                    int base = i*2;
                    bnds->oneDBnds.at(base)=std::make_shared<OneDBnd>(nGhost,nPrim,PERIODIC1D);
                    offsets=calOffsetInverse(1,i,0,bnds->iMax);
                    bnds->oneDBnds.at(base)->setUpdate(eqn->prim,offsets[0],offsets[1]);

                    bnds->oneDBnds.at(base+1)=std::make_shared<OneDBnd>(nGhost,nPrim,PERIODIC1D);
                    offsets=calOffset(1,i,0,bnds->iMax);
                    bnds->oneDBnds.at(base+1)->setUpdate(eqn->prim,offsets[0],offsets[1]);
                }
                for (int i = 0; i < iMax[0]; i++)
                {
                    int base = iMax[1]*2 + i*2;
                    bnds->oneDBnds.at(base)=std::make_shared<OneDBnd>(nGhost,nPrim,PERIODIC1D);
                    offsets=calOffsetInverse(2,i,0,bnds->iMax);
                    bnds->oneDBnds.at(base)->setUpdate(eqn->prim,offsets[0],offsets[1]);

                    bnds->oneDBnds.at(base+1)=std::make_shared<OneDBnd>(nGhost,nPrim,PERIODIC1D);
                    offsets=calOffset(2,i,0,bnds->iMax);
                    bnds->oneDBnds.at(base+1)->setUpdate(eqn->prim,offsets[0],offsets[1]);
                }
                break;
            }

            for (int i = 0; i < iMax[1]; i++)
            {
                bnds->oneDBnds.at(2*i)=std::make_shared<OneDBnd>(nGhost,nPrim,Xtype);
                offsets=calOffset(1,i,0,bnds->iMax);
                bnds->oneDBnds.at(2*i)->setUpdate(eqn->prim,offsets[0],offsets[1]);

                bnds->oneDBnds.at(2*i+1)=std::make_shared<OneDBnd>(nGhost,nPrim,Xtype);
                offsets=calOffsetInverse(1,i,0,bnds->iMax);
                bnds->oneDBnds.at(2*i+1)->setUpdate(eqn->prim,offsets[0],offsets[1]);
            }

            for (int i = 0; i < iMax[0]; i++)
            {
                bnds->oneDBnds.at(2*i+iMax[1]*2)=std::make_shared<OneDBnd>(nGhost,nPrim,Ytype);
                offsets=calOffset(2,i,0,bnds->iMax);
                bnds->oneDBnds.at(2*i+iMax[1]*2)->setUpdate(eqn->prim,offsets[0],offsets[1]);

                if(info->nCase==3) //For R-T instability
                {
                    std::array<real,4> dirVar={2.0,0,0,1.0};
                    std::vector<real> dirVars(nGhost*nPrim);
                    for(int j=0;j<nGhost;j++) {
                        dirVars.at(j*nPrim+0)=dirVar[0];dirVars.at(j*nPrim+1)=dirVar[1];
                        dirVars.at(j*nPrim+2)=dirVar[2];dirVars.at(j*nPrim+3)=dirVar[3];
                    }
                    bnds->oneDBnds.at(2*i+iMax[1]*2)->setValue(dirVars);
                }

                bnds->oneDBnds.at(2*i+1+iMax[1]*2)=std::make_shared<OneDBnd>(nGhost,nPrim,Ytype);
                offsets=calOffsetInverse(2,i,0,bnds->iMax);
                bnds->oneDBnds.at(2*i+1+iMax[1]*2)->setUpdate(eqn->prim,offsets[0],offsets[1]);

                if(info->nCase==3) //For R-T instability
                {
                    std::array<real,4> dirVar={1.0,0,0,2.5};
                    std::vector<real> dirVars(nGhost*nPrim);
                    for(int j=0;j<nGhost;j++) {
                        dirVars.at(j*nPrim+0)=dirVar[0];dirVars.at(j*nPrim+1)=dirVar[1];
                        dirVars.at(j*nPrim+2)=dirVar[2];dirVars.at(j*nPrim+3)=dirVar[3];
                    }
                    bnds->oneDBnds.at(2*i+iMax[1]*2+1)->setValue(dirVars);
                }
            }

        }
        else if (eqn->dim==3)
        {
            // 三维：六面全周期（任务 6 新增；装配方式与 2D 同型：
            // 每条线两端各一个 OneDBnd，偶数下标=左（取域尾反向值），奇数下标=右）
            std::array<int,2> offsets;
            // X 向面
            for (int i = 0; i < iMax[1]; i++)
            for (int j = 0; j < iMax[2]; j++)
            {
                int base = (i + j*iMax[1])*2;                                // 与 getOneDBnd 一致（任务 8.6 修复 stride）
                bnds->oneDBnds.at(base)=std::make_shared<OneDBnd>(nGhost,nPrim,PERIODIC1D);
                offsets=calOffsetInverse(1,i,j,bnds->iMax);
                bnds->oneDBnds.at(base)->setUpdate(eqn->prim,offsets[0],offsets[1]);

                bnds->oneDBnds.at(base+1)=std::make_shared<OneDBnd>(nGhost,nPrim,PERIODIC1D);
                offsets=calOffset(1,i,j,bnds->iMax);
                bnds->oneDBnds.at(base+1)->setUpdate(eqn->prim,offsets[0],offsets[1]);
            }
            // Y 向面
            for (int i = 0; i < iMax[0]; i++)
            for (int j = 0; j < iMax[2]; j++)
            {
                int base = iMax[1]*iMax[2]*2 + (i + j*iMax[0])*2;            // 与 getOneDBnd 一致（任务 8.6 修复 stride）
                bnds->oneDBnds.at(base)=std::make_shared<OneDBnd>(nGhost,nPrim,PERIODIC1D);
                offsets=calOffsetInverse(2,i,j,bnds->iMax);
                bnds->oneDBnds.at(base)->setUpdate(eqn->prim,offsets[0],offsets[1]);

                bnds->oneDBnds.at(base+1)=std::make_shared<OneDBnd>(nGhost,nPrim,PERIODIC1D);
                offsets=calOffset(2,i,j,bnds->iMax);
                bnds->oneDBnds.at(base+1)->setUpdate(eqn->prim,offsets[0],offsets[1]);
            }
            // Z 向面
            for (int i = 0; i < iMax[0]; i++)
            for (int j = 0; j < iMax[1]; j++)
            {
                int base = (iMax[1]*iMax[2] + iMax[0]*iMax[2])*2 + (i + j*iMax[0])*2;   // 与 getOneDBnd 一致（任务 8.6 修复 stride）
                bnds->oneDBnds.at(base)=std::make_shared<OneDBnd>(nGhost,nPrim,PERIODIC1D);
                offsets=calOffsetInverse(3,i,j,bnds->iMax);
                bnds->oneDBnds.at(base)->setUpdate(eqn->prim,offsets[0],offsets[1]);

                bnds->oneDBnds.at(base+1)=std::make_shared<OneDBnd>(nGhost,nPrim,PERIODIC1D);
                offsets=calOffset(3,i,j,bnds->iMax);
                bnds->oneDBnds.at(base+1)->setUpdate(eqn->prim,offsets[0],offsets[1]);
            }
        }
        

        break;
    
    default:
        break;
    }
    
}

void Initializer::initDoubleMachBnds(Bnds* bnds,Equation* eqn,std::array<int,3> iMax,Block* block)
{
        std::array<int,2> offsets;
        int nGhost=info->nGhostCell();
        int nCons=info->nCons();
        int nPrim=info->nPrim();
        for (int i = 0; i < iMax[1]; i++)
        {
            offsets=calOffset(1,i,0,bnds->iMax);
            real x=block->coorCel(offsets[0],0);
            bnds->oneDBnds.at(2*i)=std::make_shared<OneDBnd>(nGhost,nPrim,DIRICLET);
            bnds->oneDBnds.at(2*i)->setUpdate(eqn->prim,offsets[0],offsets[1]);

            //set values for dirichlet boundary
            std::array<real,4> dirVar={8.0,8.25*cos(M_PI/6),-8.25*sin(M_PI/6),116.5};
            std::vector<real> dirVars(nGhost*nPrim);
            for(int j=0;j<nGhost;j++) {
                dirVars.at(j*nPrim+0)=dirVar[0];dirVars.at(j*nPrim+1)=dirVar[1];
                dirVars.at(j*nPrim+2)=dirVar[2];dirVars.at(j*nPrim+3)=dirVar[3];
            }
            bnds->oneDBnds.at(2*i)->setValue(dirVars);

            bnds->oneDBnds.at(2*i+1)=std::make_shared<OneDBnd>(nGhost,nPrim,SUPERSONICOUTLET);
            offsets=calOffsetInverse(1,i,0,bnds->iMax);
            bnds->oneDBnds.at(2*i+1)->setUpdate(eqn->prim,offsets[0],offsets[1]);
        }

        for (int i = 0; i < iMax[0]; i++)
        {
            offsets=calOffset(2,i,0,bnds->iMax);
            real x=block->coorCel(offsets[0],0);
            BndType Ytype=SYMMETRYY;
            if(x<1.0/6.0) Ytype=DIRICLET;
            bnds->oneDBnds.at(2*i+iMax[1]*2)=std::make_shared<OneDBnd>(nGhost,nPrim,Ytype);
            bnds->oneDBnds.at(2*i+iMax[1]*2)->setUpdate(eqn->prim,offsets[0],offsets[1]);

            if(Ytype==DIRICLET)
            {
                std::array<real,4> dirVar={8.0,8.25*cos(M_PI/6),-8.25*sin(M_PI/6),116.5};
                std::vector<real> dirVars(nGhost*nPrim);
                for(int j=0;j<nGhost;j++) {
                    dirVars.at(j*nPrim+0)=dirVar[0];dirVars.at(j*nPrim+1)=dirVar[1];
                    dirVars.at(j*nPrim+2)=dirVar[2];dirVars.at(j*nPrim+3)=dirVar[3];
                }
                bnds->oneDBnds.at(2*i+iMax[1]*2)->setValue(dirVars);
            }

            bnds->oneDBnds.at(2*i+1+iMax[1]*2)=std::make_shared<OneDBnd>(nGhost,nPrim,DoubleMachUp);
            offsets=calOffsetInverse(2,i,0,bnds->iMax);
            bnds->oneDBnds.at(2*i+1+iMax[1]*2)->setUpdate(eqn->prim,offsets[0],offsets[1]);

            std::array<real,3> coor={block->coorCel(offsets[0],0),
                                    block->coorCel(offsets[0],1),
                                    block->coorCel(offsets[0],2)};
            std::array<real,3> dh={block->coorCel(offsets[0],0)-block->coorCel(offsets[0]+offsets[1],0),
                                   block->coorCel(offsets[0],1)-block->coorCel(offsets[0]+offsets[1],1),
                                   block->coorCel(offsets[0],2)-block->coorCel(offsets[0]+offsets[1],2)};
            bnds->oneDBnds.at(2*i+1+iMax[1]*2)->setInfo(info);
            bnds->oneDBnds.at(2*i+1+iMax[1]*2)->setCoor(coor,dh);
        }
}

void Initializer::initShockShearBnds(Bnds* bnds,Equation* eqn,std::array<int,3> iMax,Block* block)
{
        std::array<int,2> offsets;
        int nGhost=info->nGhostCell();
        int nPrim=info->nPrim();
        // X 向：左=带时变扰动的入流，右=超声速出口
        for (int i = 0; i < iMax[1]; i++)
        {
            offsets=calOffset(1,i,0,bnds->iMax);
            bnds->oneDBnds.at(2*i)=std::make_shared<OneDBnd>(nGhost,nPrim,ShockShearIn);
            bnds->oneDBnds.at(2*i)->setUpdate(eqn->prim,offsets[0],offsets[1]);

            std::array<real,3> coor={block->coorCel(offsets[0],0),
                                    block->coorCel(offsets[0],1),
                                    block->coorCel(offsets[0],2)};
            std::array<real,3> dh={block->coorCel(offsets[0],0)-block->coorCel(offsets[0]+offsets[1],0),
                                   block->coorCel(offsets[0],1)-block->coorCel(offsets[0]+offsets[1],1),
                                   block->coorCel(offsets[0],2)-block->coorCel(offsets[0]+offsets[1],2)};
            bnds->oneDBnds.at(2*i)->setInfo(info);
            bnds->oneDBnds.at(2*i)->setCoor(coor,dh);

            offsets=calOffsetInverse(1,i,0,bnds->iMax);
            bnds->oneDBnds.at(2*i+1)=std::make_shared<OneDBnd>(nGhost,nPrim,SUPERSONICOUTLET);
            bnds->oneDBnds.at(2*i+1)->setUpdate(eqn->prim,offsets[0],offsets[1]);
        }

        // Y 向：下=滑移固壁（对称镜像），上=后激波态（固定值）
        for (int i = 0; i < iMax[0]; i++)
        {
            offsets=calOffset(2,i,0,bnds->iMax);
            bnds->oneDBnds.at(2*i+iMax[1]*2)=std::make_shared<OneDBnd>(nGhost,nPrim,SYMMETRYY);
            bnds->oneDBnds.at(2*i+iMax[1]*2)->setUpdate(eqn->prim,offsets[0],offsets[1]);

            offsets=calOffsetInverse(2,i,0,bnds->iMax);
            bnds->oneDBnds.at(2*i+1+iMax[1]*2)=std::make_shared<OneDBnd>(nGhost,nPrim,DIRICLET);
            bnds->oneDBnds.at(2*i+1+iMax[1]*2)->setUpdate(eqn->prim,offsets[0],offsets[1]);

            std::array<real,4> dirVar={2.1101,2.9709,-0.1367,0.4754};
            std::vector<real> dirVars(nGhost*nPrim);
            for(int j=0;j<nGhost;j++) {
                dirVars.at(j*nPrim+0)=dirVar[0];dirVars.at(j*nPrim+1)=dirVar[1];
                dirVars.at(j*nPrim+2)=dirVar[2];dirVars.at(j*nPrim+3)=dirVar[3];
            }
            bnds->oneDBnds.at(2*i+1+iMax[1]*2)->setValue(dirVars);
        }
}

void Initializer::initSpDistributor(SpDistributor* spDis,Equation* eqn
                                    ,Block* block,Bnds* bnds)
{
    
    spDis->nCons=info->nCons();
    spDis->nPrim=info->nPrim();
    spDis->iMax=block->icMax;
    spDis->prim=eqn->prim;
    spDis->cons=eqn->cons;
    spDis->dim=info->dim;
    spDis->bnds=bnds;
    spDis->rhs=eqn->rhs;
    spDis->info=info;
    

}