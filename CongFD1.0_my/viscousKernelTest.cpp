/**
 * @file viscousKernelTest.cpp
 * @brief 粘性模块内核测试（任务 4）
 *
 * 组 1：物性数值（μ=1/Re、κ=μc_p/Pr）
 * 组 2：CD6 对五次多项式精确（机器精度）
 * 组 3：CD6 梯度场 sin 收敛 6 阶（整场含周期回绕，16/32/64）
 * 组 4：粘性 rhs 解析锚点（单向 u=sin x：动量与能量分量逐点对照）
 * 组 4b：幂律粘性解析锚点（均匀温度场，因子=1 与 2^0.75；HIT 开发新增）
 * 组 5：均匀场 → 粘性贡献为零
 * （组 6 至 8：平面波族衰减阶数，S8 追加）
 */
// 注：STL 头须先于项目头（macro.hpp 的 #define real double 会破坏 <complex> 的成员名 real）
#include <complex>
#include <cstdio>
#include <cmath>
#include <vector>
#include "ViscousTerm.hpp"
#include "blockSolver.hpp"
#include "diagnostics.hpp"

static int nFail=0;
static void check(bool ok,const char* name)
{
    std::printf("%s : %s\n", ok? "PASS":"FAIL", name);
    if(!ok) nFail++;
}

static const real PI=3.14159265358979323846;

// 构造 N^3 场并按给定函数填 prim（格心坐标；域 [0,2π]^3）
template<class F>
static void fillField(Data& prim,int N,real h,F f)
{
    for(int k=0;k<N;k++)
    for(int j=0;j<N;j++)
    for(int i=0;i<N;i++)
    {
        int idx=i+j*N+k*N*N;
        real x=(i+0.5)*h,y=(j+0.5)*h,z=(k+0.5)*h;
        real r,u,v,w,p;
        f(x,y,z,r,u,v,w,p);
        prim(idx,0)=r; prim(idx,1)=u; prim(idx,2)=v; prim(idx,3)=w; prim(idx,4)=p;
    }
}

static void setupInfo(Info& info,int N,real L=2*PI)
{
    info.eqType=EULER;
    info.dim=3;
    info.iMax={N+1,N+1,N+1};
    info.calZone={0,L,0,L,0,L};
    info.viscous=true;
    info.Re=1600;
    info.Pr=0.72;
}

// 平面波复投影：A = |Σ f e^{−i k·x}| / (n/2)；k_d = 2π m_d/L_d（格心；minusOne 用于 ρ′）
static real modeAmp(Data* f,int ivar,const std::array<int,3>& Nvec,
                    const std::array<int,3>& m,bool minusOne)
{
    int ni=Nvec[0]*Nvec[1]*Nvec[2];
    std::complex<real> s(0.0,0.0);
    for(int k=0;k<Nvec[2];k++)
    for(int j=0;j<Nvec[1];j++)
    for(int i=0;i<Nvec[0];i++)
    {
        int idx=i+j*Nvec[0]+k*Nvec[0]*Nvec[1];
        real ph=2*PI*(m[0]*(i+0.5)/Nvec[0]+m[1]*(j+0.5)/Nvec[1]+m[2]*(k+0.5)/Nvec[2]);
        real val=(*f)(idx,ivar)-(minusOne?1.0:0.0);
        s+=val*std::complex<real>(std::cos(ph),-std::sin(ph));
    }
    return std::abs(s)/(ni/2.0);
}

int main()
{
    // ---- 组 1：物性数值 ----
    {
        int N=16; real h=2*PI/N;
        Info info; setupInfo(info,N);
        Data prim(N*N*N,5), rhs(N*N*N,5);
        ViscousTerm vt(&prim,&rhs,&info);
        real e=0;
        e=std::max(e,std::abs(vt.getMu()-1.0/1600.0));
        e=std::max(e,std::abs(vt.getKappa()-(1.0/1600.0)*(GAMMA/(GAMMA-1.0))/0.72));
        e=std::max(e,std::abs(vt.getCp()-3.5));
        check(e<1e-15,"组1 物性: mu=6.25e-4, kappa=3.038e-3, cp=3.5");
    }

    // ---- 组 2：CD6 五次多项式精确 ----
    {
        real h=0.1,x0=0.7;
        auto f=[&](real m){ real x=x0+m*h; return x*x*x*x*x; };
        real d=visc::cd6(f(-3),f(-2),f(-1),f(1),f(2),f(3),h);
        check(std::abs(d-5*x0*x0*x0*x0)<1e-10,"组2 CD6 五次多项式导数（机器精度级）");
    }

    // ---- 组 3：CD6 梯度场 sin 收敛 6 阶 ----
    {
        int Ns[3]={16,32,64};
        real err[3];
        for(int t=0;t<3;t++)
        {
            int N=Ns[t]; real h=2*PI/N;
            Data prim(N*N*N,5);
            fillField(prim,N,h,[&](real x,real y,real z,real&r,real&u,real&v,real&w,real&p)
                     { r=1; u=std::sin(x)*std::cos(y)*std::cos(z); v=0; w=0; p=100; });
            std::vector<real> gu((size_t)3*N*N*N);
            std::array<int,3> icM{N,N,N};
            for(int d=0;d<3;d++) visc::gradField(&prim,1,gu,d,icM,h);
            real e=0;
            for(int k=0;k<N;k++)
            for(int j=0;j<N;j++)
            for(int i=0;i<N;i++)
            {
                int idx=i+j*N+k*N*N;
                real x=(i+0.5)*h,y=(j+0.5)*h,z=(k+0.5)*h;
                real ex=std::cos(x)*std::cos(y)*std::cos(z);
                real ey=-std::sin(x)*std::sin(y)*std::cos(z);
                real ez=-std::sin(x)*std::cos(y)*std::sin(z);
                e=std::max(e,std::abs(gu[0*(size_t)N*N*N+idx]-ex));
                e=std::max(e,std::abs(gu[1*(size_t)N*N*N+idx]-ey));
                e=std::max(e,std::abs(gu[2*(size_t)N*N*N+idx]-ez));
            }
            err[t]=e;
        }
        std::printf("   组3 梯度误差 N=16/32/64: %.3e / %.3e / %.3e\n",err[0],err[1],err[2]);
        check(err[0]/err[1]>40.0 && err[1]/err[2]>40.0,"组3 CD6 梯度 sin 场 6 阶收敛（比值>40）");
    }

    // ---- 组 4：粘性 rhs 解析锚点（单向 u=sin x） ----
    {
        int N=64; real h=2*PI/N;
        Info info; setupInfo(info,N);
        Data prim(N*N*N,5), rhs(N*N*N,5);
        rhs.setZeros();
        fillField(prim,N,h,[&](real x,real y,real z,real&r,real&u,real&v,real&w,real&p)
                 { r=1; u=std::sin(x); v=0; w=0; p=100; });
        ViscousTerm vt(&prim,&rhs,&info);
        vt.calViscous();
        real mu=1.0/1600.0, coef=mu*4.0/3.0;
        real em=0,ee=0,eo=0;
        for(int k=0;k<N;k++)
        for(int j=0;j<N;j++)
        for(int i=0;i<N;i++)
        {
            int idx=i+j*N+k*N*N;
            real x=(i+0.5)*h;
            real want_m=coef*std::sin(x);                  // rhs_u = −∂_x τ_xx = μ(4/3) sin x
            real want_e=-coef*std::cos(2.0*x);             // rhs_E = −∂_x(u·τ_xx) = −μ(4/3) cos 2x
            em=std::max(em,std::abs(rhs(idx,1)-want_m));
            ee=std::max(ee,std::abs(rhs(idx,4)-want_e));
            eo=std::max(eo,std::abs(rhs(idx,0)));
            eo=std::max(eo,std::abs(rhs(idx,2)));
            eo=std::max(eo,std::abs(rhs(idx,3)));
        }
        std::printf("   组4 误差 动量=%.3e 能量=%.3e 其余=%.3e\n",em,ee,eo);
        check(em<1e-9 && ee<1e-9 && eo<1e-14,"组4 粘性 rhs 解析锚点（单向 u=sin x）");
    }

    // ---- 组 4b：幂律粘性解析锚点（均匀温度场；HIT 开发新增） ----
    // T=100=Tref → 因子 1；T=200 → 因子 2^0.75。物性均匀，右端项应精确等于常数因子放大的组 4 解析值。
    {
        int N=64; real h=2*PI/N;
        const real coefBase=(1.0/1600.0)*4.0/3.0;
        const real ps[2]={100.0,200.0};
        const real facs[2]={1.0,std::pow(2.0,0.75)};
        for(int c=0;c<2;c++)
        {
            Info info; setupInfo(info,N);
            info.powerLawExp=0.75; info.Tref=100.0;
            Data prim(N*N*N,5), rhs(N*N*N,5);
            rhs.setZeros();
            real p0=ps[c];
            fillField(prim,N,h,[&](real x,real y,real z,real&r,real&u,real&v,real&w,real&p)
                     { r=1; u=std::sin(x); v=0; w=0; p=p0; });
            ViscousTerm vt(&prim,&rhs,&info);
            vt.calViscous();
            real coef=coefBase*facs[c];
            real e=0;
            for(int k=0;k<N;k++)
            for(int j=0;j<N;j++)
            for(int i=0;i<N;i++)
            {
                int idx=i+j*N+k*N*N;
                real x=(i+0.5)*h;
                real want_m=coef*std::sin(x);              // rhs_u = μ·fac·(4/3)·sin x
                real want_e=-coef*std::cos(2.0*x);         // rhs_E = −μ·fac·(4/3)·cos 2x
                e=std::max(e,std::abs(rhs(idx,1)-want_m));
                e=std::max(e,std::abs(rhs(idx,4)-want_e));
            }
            char nm[128];
            std::snprintf(nm,sizeof(nm),"组4b 幂律粘性解析锚点（T=%g，因子=%.4f）",(double)p0,(double)facs[c]);
            check(e<1e-9,nm);
        }
    }

    // ---- 组 5：均匀场 → 粘性贡献为零 ----
    {
        int N=16; real h=2*PI/N;
        Info info; setupInfo(info,N);
        Data prim(N*N*N,5), rhs(N*N*N,5);
        rhs.setZeros();
        fillField(prim,N,h,[&](real x,real y,real z,real&r,real&u,real&v,real&w,real&p)
                 { r=1; u=1; v=2; w=3; p=100; });
        ViscousTerm vt(&prim,&rhs,&info);
        vt.calViscous();
        real e=0;
        for(int i=0;i<N*N*N;i++)
        for(int c=0;c<5;c++) e=std::max(e,std::abs(rhs(i,c)));
        check(e<1e-12,"组5 均匀场粘性贡献为零");
    }

    // ---- 组 6：剪切波（模块级斜波；ν|k|²）：各向同性 + 非等边网格 ----
    {
        real A=1e-3;
        real betaTh=(1.0/1600.0)*3.0;                 // ν|k|²，k=(1,1,1)（L=2π、m=(1,1,1)）
        auto run=[&](const std::array<int,3>& Nvec,real& beta)
        {
            std::array<real,3> hh={2*PI/Nvec[0],2*PI/Nvec[1],2*PI/Nvec[2]};
            int ni=Nvec[0]*Nvec[1]*Nvec[2];
            Info info; info.eqType=EULER; info.dim=3;
            info.iMax={Nvec[0]+1,Nvec[1]+1,Nvec[2]+1};
            info.calZone={0,2*PI,0,2*PI,0,2*PI};
            info.viscous=true; info.Re=1600; info.Pr=0.72;
            Data prim(ni,5), rhs(ni,5); rhs.setZeros();
            for(int k=0;k<Nvec[2];k++)
            for(int j=0;j<Nvec[1];j++)
            for(int i=0;i<Nvec[0];i++)
            {
                int idx=i+j*Nvec[0]+k*Nvec[0]*Nvec[1];
                real x=(i+0.5)*hh[0],y=(j+0.5)*hh[1],z=(k+0.5)*hh[2];
                real wv=std::sin(x+y+z);              // 斜波 k=(1,1,1)
                prim(idx,0)=1.0;
                prim(idx,1)= A/std::sqrt(2.0)*wv;     // 极化 (1,−1,0)/√2：∇·u=0
                prim(idx,2)=-A/std::sqrt(2.0)*wv;
                prim(idx,3)=0.0;
                prim(idx,4)=100.0;
            }
            ViscousTerm vt(&prim,&rhs,&info);
            vt.calViscous();
            // 斜波为粘性算子本征模：衰减率 = ⟨rhs,u⟩/⟨u,u⟩（免时间步偏差）
            real num=0,den=0;
            for(int idx=0;idx<ni;idx++)
            {
                real uu=prim(idx,1), vv=prim(idx,2);
                num+=rhs(idx,1)*uu+rhs(idx,2)*vv;
                den+=uu*uu+vv*vv;
            }
            beta=num/den;
        };
        int Ns[3]={16,32,64}; real err[3];
        for(int t=0;t<3;t++)
        {
            real beta; std::array<int,3> Nv={Ns[t],Ns[t],Ns[t]};
            run(Nv,beta);
            err[t]=std::abs(beta-betaTh)/betaTh;
        }
        std::printf("   组6a 剪切衰减率相对误差: %.3e / %.3e / %.3e\n",err[0],err[1],err[2]);
        check(err[0]/err[1]>40.0 && err[1]/err[2]>40.0,"组6a 剪切波 ν|k|² 6 阶（各向同性斜波）");
        std::array<std::array<int,3>,3> tri={{{16,12,8},{32,24,16},{64,48,32}}};
        for(int t=0;t<3;t++)
        {
            real beta;
            run(tri[t],beta);
            err[t]=std::abs(beta-betaTh)/betaTh;
        }
        std::printf("   组6b 非等边网格相对误差: %.3e / %.3e / %.3e\n",err[0],err[1],err[2]);
        check(err[0]/err[1]>40.0 && err[1]/err[2]>40.0,"组6b 剪切波 ν|k|² 6 阶（斜波＋非等边网格）");
    }

    // ---- 组 7：熵波（求解器级；长时开−关差分剥离格式耗散）：κ/(ρc_p)·k² ----
    {
        real eps=1e-3;
        real betaTh=(1.0/1600.0)/0.72*4.0;            // κ/(ρc_p)·k²，k=2
        int Ns[3]={16,24,32}; real err[3], ord[2];
        auto slopeFit=[&](const std::vector<real>& As,const std::vector<real>& Ts,real t1,real t2)
        {
            real Sx=0,Sy=0,Sxx=0,Sxy=0,nn=0;
            for(size_t i=0;i<As.size();i++)
            {
                if(Ts[i]<t1||Ts[i]>t2) continue;
                real xx=Ts[i],yy=std::log(As[i]);
                Sx+=xx;Sy+=yy;Sxx+=xx*xx;Sxy+=xx*yy;nn+=1.0;
            }
            return (nn*Sxy-Sx*Sy)/(nn*Sxx-Sx*Sx);
        };
        auto run=[&](int N,bool visc,std::vector<real>& As,std::vector<real>& Ts)
        {
            real c=std::sqrt(GAMMA*100.0);
            real h=2*PI/N;
            real dt=0.4*h/c;                          // 稳定余量内固定步长
            real Ttot=50.0;
            Info* info=new Info;
            info->eqType=EULER; info->dim=3; info->nCase=3;
            info->iMax={N+1,9,9}; info->calZone={0,2*PI,0,2*PI,0,2*PI};
            info->interMethod=(InterMethod)31;
            info->diffMethod=MND6;                // 与生产一致（zoneMain.cpp:312）
            info->viscous=visc; info->Re=1600; info->Pr=0.72;
            BlockSolver bs(info);
            Data* cons=bs.getConsPtr();
            for(int k=0;k<8;k++)
            for(int j=0;j<8;j++)
            for(int i=0;i<N;i++)
            {
                int idx=i+j*N+k*N*8;
                real x=(i+0.5)*(2*PI/N);
                real r=1.0+eps*std::sin(2.0*x);       // 熵模初场：p 均匀、u=0（k=2）
                (*cons)(idx,0)=r;
                (*cons)(idx,1)=0.0; (*cons)(idx,2)=0.0; (*cons)(idx,3)=0.0;
                (*cons)(idx,4)=100.0/(GAMMA-1.0);
            }
            std::array<int,3> Nv={N,8,8}, m={2,0,0};
            long nstep=(long)(Ttot/dt);
            for(long s=0;s<nstep;s++)
            {
                bs.solve(dt);
                As.push_back(modeAmp(cons,0,Nv,m,true));
                Ts.push_back(info->t);
            }
        };
        for(int t=0;t<3;t++)
        {
            std::vector<real> Ao,To,Aon,Ton;
            run(Ns[t],false,Ao,To);
            run(Ns[t],true ,Aon,Ton);
            real bv=slopeFit(Ao,To,25.0,50.0)-slopeFit(Aon,Ton,25.0,50.0);
            err[t]=std::abs(bv-betaTh)/betaTh;
            std::printf("   组7 N=%d: β_visc=%.6e (理论 %.6e) 相对误差 %.3e\n",Ns[t],bv,betaTh,err[t]);
        }
        ord[0]=std::log(err[0]/err[1])/std::log(2.0);
        ord[1]=std::log(err[1]/err[2])/std::log(1.5);
        std::printf("   组7 观测阶: %.2f / %.2f\n",ord[0],ord[1]);
        check(ord[0]>3.0 && ord[1]>3.0,"组7 熵波 κ/(ρc_p)k² 收敛（≥3.0；受格式有效阶与测量窗口残差限制，见 04 实施记录 §5）");
    }

    // ---- 组 8：声波（求解器级；长时开−关差分）：Stokes 衰减 ----
    {
        real eps=1e-3;
        real c=std::sqrt(GAMMA*100.0);                // c²=γp0/ρ0
        real betaTh=(1.0/1600.0)/2.0*(4.0/3.0+(GAMMA-1.0)/0.72)*4.0;  // (νk²/2)(4/3+(γ−1)/Pr)，k=2
        int Ns[3]={16,24,32}; real err[3], ord[2];
        auto slopeFit=[&](const std::vector<real>& As,const std::vector<real>& Ts,real t1,real t2)
        {
            real Sx=0,Sy=0,Sxx=0,Sxy=0,nn=0;
            for(size_t i=0;i<As.size();i++)
            {
                if(Ts[i]<t1||Ts[i]>t2) continue;
                real xx=Ts[i],yy=std::log(As[i]);
                Sx+=xx;Sy+=yy;Sxx+=xx*xx;Sxy+=xx*yy;nn+=1.0;
            }
            return (nn*Sxy-Sx*Sy)/(nn*Sxx-Sx*Sx);
        };
        auto run=[&](int N,bool visc,std::vector<real>& As,std::vector<real>& Ts)
        {
            real h=2*PI/N;
            real dt=0.4*h/c;
            real Ttot=50.0;
            real bk=betaTh/2.0;                       // 声模特征矢 u 分量的黏性相位修正（β/k，k=2）
            Info* info=new Info;
            info->eqType=EULER; info->dim=3; info->nCase=3;
            info->iMax={N+1,9,9}; info->calZone={0,2*PI,0,2*PI,0,2*PI};
            info->interMethod=(InterMethod)31;
            info->diffMethod=MND6;                // 与生产一致（zoneMain.cpp:312）
            info->viscous=visc; info->Re=1600; info->Pr=0.72;
            BlockSolver bs(info);
            Data* cons=bs.getConsPtr();
            for(int k=0;k<8;k++)
            for(int j=0;j<8;j++)
            for(int i=0;i<N;i++)
            {
                int idx=i+j*N+k*N*8;
                real x=(i+0.5)*(2*PI/N);
                real dr=eps*std::sin(2.0*x);          // 右行声模（k=2）：ρ'=ε sin2x、p'=c²ρ'、u=ε(c sin2x − (β/k)cos2x)
                real r=1.0+dr;
                real p=100.0+c*c*dr;
                real u=eps*(c*std::sin(2.0*x)-bk*std::cos(2.0*x));
                (*cons)(idx,0)=r;
                (*cons)(idx,1)=r*u;
                (*cons)(idx,2)=0.0; (*cons)(idx,3)=0.0;
                (*cons)(idx,4)=p/(GAMMA-1.0)+0.5*r*u*u;
            }
            std::array<int,3> Nv={N,8,8}, m={2,0,0};
            long nstep=(long)(Ttot/dt);
            for(long s=0;s<nstep;s++)
            {
                bs.solve(dt);
                As.push_back(modeAmp(cons,0,Nv,m,true));
                Ts.push_back(info->t);
            }
        };
        for(int t=0;t<3;t++)
        {
            std::vector<real> Ao,To,Aon,Ton;
            run(Ns[t],false,Ao,To);
            run(Ns[t],true ,Aon,Ton);
            real bv=slopeFit(Ao,To,25.0,50.0)-slopeFit(Aon,Ton,25.0,50.0);
            err[t]=std::abs(bv-betaTh)/betaTh;
            std::printf("   组8 N=%d: β_visc=%.6e (理论 %.6e) 相对误差 %.3e\n",Ns[t],bv,betaTh,err[t]);
        }
        ord[0]=std::log(err[0]/err[1])/std::log(2.0);
        ord[1]=std::log(err[1]/err[2])/std::log(1.5);
        std::printf("   组8 观测阶: %.2f / %.2f\n",ord[0],ord[1]);
        check(err[1]<5.0e-3,"组8 声波 Stokes 衰减定量核对（N=24 相对误差 <0.5%；网格收敛由组6与 MMS 承担，见 04 实施记录 §5）");
    }

    // ---- 组 9：恒等式 ε=−dK/dt ≈ 2νΩ（模块级；不可压恒等，低马赫偏差 O(Ma²)）----
    {
        int N=32; real h=2*PI/N;
        Info info; setupInfo(info,N);
        Data prim(N*N*N,5), rhs(N*N*N,5);
        rhs.setZeros();
        fillField(prim,N,h,[&](real x,real y,real z,real&r,real&u,real&v,real&w,real&p)
                 { r=1; u=std::sin(x)*std::cos(y)*std::cos(z);
                   v=-std::cos(x)*std::sin(y)*std::cos(z); w=0;
                   p=100.0+((std::cos(2*x)+std::cos(2*y))*(2.0+std::cos(2*z))-2.0)/16.0; });
        real K0,Om0,U2d,Tvard,Thvard;
        calcDiagnostics(&prim,{N,N,N},h,K0,Om0,U2d,Tvard,Thvard);
        ViscousTerm vt(&prim,&rhs,&info);
        vt.calViscous();
        real sum=0;
        for(int idx=0;idx<N*N*N;idx++)
            sum+=(*&prim)(idx,1)*rhs(idx,1)+(*&prim)(idx,2)*rhs(idx,2)+(*&prim)(idx,3)*rhs(idx,3);
        real epsHat=sum*h*h*h;                     // ε=−dK/dt=Σ u·rhs·h³（∂_t(ρu)=−rhs；含格体积）
        real twoNuOm=2.0*(1.0/1600.0)*Om0;         // ρ=1：ν=μ
        real rel=std::abs(epsHat-twoNuOm)/twoNuOm;
        std::printf("   组9 ε=%.6e vs 2νΩ=%.6e  相对差=%.3e（O(Ma²)≈7e-3）\n",epsHat,twoNuOm,rel);
        check(rel<2.0e-2,"组9 恒等式 ε=2νΩ（模块级，低马赫 O(Ma²) 内）");
    }

    std::printf("== %s (failures=%d) ==\n", nFail==0? "ALL PASS":"FAIL", nFail);
    return nFail==0? 0:1;
}
