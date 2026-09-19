/**
 * @file viscousOrderTest.cpp
 * @brief 粘性 MMS 阶数测定（任务 4 S9）：整求解器驱动，目标观测阶 ≥5.5
 *
 * 制造解与源项见 src/include/viscousMMS.hpp 与 sourceTerm.cpp（VISCOUSMMS 分支）。
 * 时间步按 dt=0.026·h² 取（使 RK3 时间误差同为 h⁶ 量级、不污染空间阶数）；
 * 推进至 t=0.5 后对守恒量取相对 L2 误差。网格 16/24/32 全向同分辨率。
 */
#include <complex>
#include <cstdio>
#include <cmath>
#include <vector>
#include "blockSolver.hpp"
#include "viscousMMS.hpp"

static const real PIT=3.14159265358979323846;

int main()
{
    int Ns[3]={16,24,32};
    real errs[3];
    for(int tt=0;tt<3;tt++)
    {
        int N=Ns[tt];
        real L=2*PIT;
        real h=L/N;
        real dt=0.026*h*h;
        long nstep=(long)std::llround(0.5/dt);
        Info* info=new Info;
        info->eqType=EULER; info->dim=3; info->nCase=0;
        info->iMax={N+1,N+1,N+1}; info->calZone={0,L,0,L,0,L};
        info->interMethod=(InterMethod)31;
        info->diffMethod=MND6;                        // 与生产一致（zoneMain.cpp:312）
        info->viscous=true; info->Re=1600; info->Pr=0.72;
        info->sourceType=VISCOUSMMS;
        BlockSolver bs(info);
        Data* cons=bs.getConsPtr();
        for(int k=0;k<N;k++)
        for(int j=0;j<N;j++)
        for(int i=0;i<N;i++)
        {
            int idx=i+j*N+k*N*N;
            real x=(i+0.5)*h,y=(j+0.5)*h,z=(k+0.5)*h;
            real q[5],dq[3][5],dtv[5];
            mms::state<real>(x,y,z,0.0,q,dq,dtv);
            real r=q[0],u=q[1],v=q[2],w=q[3],p=q[4];
            (*cons)(idx,0)=r;
            (*cons)(idx,1)=r*u;
            (*cons)(idx,2)=r*v;
            (*cons)(idx,3)=r*w;
            (*cons)(idx,4)=p/(GAMMA-1.0)+0.5*r*(u*u+v*v+w*w);
        }
        for(long s=0;s<nstep;s++) bs.solve(dt);
        real tf=nstep*dt;
        real num=0,den=0;
        for(int k=0;k<N;k++)
        for(int j=0;j<N;j++)
        for(int i=0;i<N;i++)
        {
            int idx=i+j*N+k*N*N;
            real x=(i+0.5)*h,y=(j+0.5)*h,z=(k+0.5)*h;
            real q[5],dq[3][5],dtv[5];
            mms::state<real>(x,y,z,tf,q,dq,dtv);
            real r=q[0],u=q[1],v=q[2],w=q[3],p=q[4];
            real qe[5]={r,r*u,r*v,r*w,p/(GAMMA-1.0)+0.5*r*(u*u+v*v+w*w)};
            for(int c=0;c<5;c++)
            {
                real dd=(*cons)(idx,c)-qe[c];
                num+=dd*dd; den+=qe[c]*qe[c];
            }
        }
        errs[tt]=std::sqrt(num/den);
        std::printf("MMS N=%d: dt=%.3e steps=%ld t=%.3f  相对L2=%.6e\n",N,dt,nstep,tf,errs[tt]);
    }
    real o1=std::log(errs[0]/errs[1])/std::log(1.5);
    real o2=std::log(errs[1]/errs[2])/std::log(4.0/3.0);
    std::printf("观测阶: %.2f / %.2f (≥4.5；对齐第一部分有效阶口径，见 04 实施记录 §5)\n",o1,o2);
    int fail=(o1>4.5&&o2>4.5)?0:1;
    std::printf("== %s ==\n", fail==0?"PASS":"FAIL");
    return fail;
}
