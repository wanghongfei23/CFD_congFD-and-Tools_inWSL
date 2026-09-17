/**
 * @file dim3KernelTest.cpp
 * @brief 三维物理层内核测试（任务 3：通量组）
 *
 * 测试组 1–5（任务 7 起全部就位）：
 *   1. fEuler3D 物理通量
 *   2. roeFlux3DSym 左右同态 → 应等于物理通量
 *   3. roeFlux3DSym 的 2D 退化一致性（w=0, n=x）
 *   4. eigensystemEuler3D 往返与 L*R=I（任务 4 追加）
 *   5. K/Ω 诊断解析场锚点（任务 7 追加）
 */
#include "macro.hpp"
#include "fluxScheme.hpp"
#include "eigenSystem.hpp"
#include "diagnostics.hpp"
#include <cstdio>
#include <cmath>

static int nFail = 0;
static void check(bool ok, const char* name)
{
    std::printf("%s : %s\n", ok? "PASS":"FAIL", name);
    if(!ok) nFail++;
}

int main()
{
    {   // 1. 3D 物理通量（x 向）：逐分量对照解析式
        std::array<real,3> n{1,0,0};
        real r=1.2,u=0.3,v=-0.2,w=0.5,p=2.0;
        auto f=fEuler3D({r,u,v,w,p},n);
        real e=0;
        e=std::max(e,std::abs(f[0]-r*u));
        e=std::max(e,std::abs(f[1]-(r*u*u+p)));
        e=std::max(e,std::abs(f[2]-r*v*u));
        e=std::max(e,std::abs(f[3]-r*w*u));
        e=std::max(e,std::abs(f[4]-((u*u+v*v+w*w)/2*r+GAMMA/(GAMMA-1)*p)*u));
        check(e<1e-14, "fEuler3D physical flux");
    }
    {   // 2. Roe 3D：左右同态 → 耗散项为零，通量应等于物理通量
        std::array<real,3> n{1,0,0};
        real r=1.2,u=0.3,v=-0.2,w=0.5,p=2.0;
        auto F=roeFlux3DSym(r,r,u,u,v,v,w,w,p,p,n);
        auto f=fEuler3D({r,u,v,w,p},n);
        real e=0; for(int i=0;i<5;i++) e=std::max(e,std::abs(F[i]-f[i]));
        check(e<1e-12, "roeFlux3DSym identical states == physical flux");
    }
    {   // 3. 2D 退化一致性：w=0、n=x 时 3D 通量与 2D 通量逐值一致（ρw 通量应为 0）
        std::array<real,3> n{1,0,0};
        real rl=1.0,ul=0.75,vl=-0.5,pl=1.0;
        real rr=2.0,ur=-0.3,vr=0.2,pr=0.6;
        auto F3=roeFlux3DSym(rl,rr,ul,ur,vl,vr,0.0,0.0,pl,pr,n);
        auto F2=roeFlux2DSym(rl,rr,ul,ur,vl,vr,pl,pr,n);
        real e=0;
        e=std::max(e,std::abs(F3[0]-F2[0]));
        e=std::max(e,std::abs(F3[1]-F2[1]));
        e=std::max(e,std::abs(F3[2]-F2[2]));
        e=std::max(e,std::abs(F3[4]-F2[3]));
        e=std::max(e,std::abs(F3[3]));
        check(e<1e-12, "roeFlux3DSym w=0 degenerates to 2D");
    }
    {   // 4. 特征系统：往返 + L*R=I
        std::array<real,3> n{0,0,1};
        std::array<real,5> wl{1.0,0.3,-0.4,0.2,1.5};
        std::array<real,5> wr{1.4,-0.2,0.1,-0.6,0.8};
        eigensystemEuler3D eig(wl,wr,n);
        std::array<real,5> wm{1.2,0.05,-0.15,-0.2,1.1};
        auto c=eig.primToChar(wm);
        auto w2=eig.charToPrim(c);
        real e=0;
        for(int i=0;i<5;i++) e=std::max(e,std::abs(w2[i]-wm[i])/(std::abs(wm[i])+1e-30));
        check(e<1e-12, "eigensystem3D roundtrip prim->char->prim");
        real e2=0;
        for(int a=0;a<5;a++)
        for(int b=0;b<5;b++)
        {
            real s=0;
            for(int k=0;k<5;k++) s+=eig.leftEig[a*5+k]*eig.rightEig[k*5+b];
            e2=std::max(e2,std::abs(s-((a==b)?1.0:0.0)));
        }
        check(e2<1e-12, "eigensystem3D L*R == I");
    }
    {   // 5. K/Ω 诊断：解析 TGV 场（64 格、格心采样）→ K=pi^3、Ω=3*pi^3
        const int N=64;
        real h=2.0*M_PI/N;
        Data prim(N*N*N,5);
        for(int k=0;k<N;k++)
        for(int j=0;j<N;j++)
        for(int i=0;i<N;i++)
        {
            int id=i+j*N+k*N*N;
            real x=(i+0.5)*h,y=(j+0.5)*h,z=(k+0.5)*h;
            real u=std::sin(x)*std::cos(y)*std::cos(z);
            real v=-std::cos(x)*std::sin(y)*std::cos(z);
            real p=100.0+(std::cos(2*x)+std::cos(2*y))*(2.0+std::cos(2*z))/16.0-2.0/16.0;
            prim(id,0)=1.0; prim(id,1)=u; prim(id,2)=v; prim(id,3)=0.0; prim(id,4)=p;
        }
        real K,Omega;
        calcDiagnostics(&prim,{N,N,N},h,K,Omega);
        check(std::abs(K-M_PI*M_PI*M_PI)/(M_PI*M_PI*M_PI)<1e-10, "K(0)==pi^3");
        check(std::abs(Omega-3.0*M_PI*M_PI*M_PI)/(3.0*M_PI*M_PI*M_PI)<1e-4, "Omega(0)==3*pi^3");
    }
    std::printf("== %s (failures=%d) ==\n", nFail==0? "ALL PASS":"FAIL", nFail);
    return nFail==0? 0:1;
}
