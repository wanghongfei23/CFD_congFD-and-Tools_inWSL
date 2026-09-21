/**
 * @file viscousTerm.cpp
 * @brief 粘性项模块实现（任务 4；对应 04 实施分解 S2 物性 / S3 梯度 / S4 通量与散度）
 *
 * 符号约定：面通量差分链向 rhs 累加 +∇·F_inv；主循环 cons = temp − dt·rhs；
 * 故粘性项按 rhs −= ∇·Ev 累加（等效 ∂_t q = −∇·F_inv + ∇·F_vis）。
 * 守恒性：CD6 为反对称差分，周期回绕下全域求和精确为零。
 */
#include "ViscousTerm.hpp"
#include <iostream>

ViscousTerm::ViscousTerm(Data* prim_,Data* rhs_,Info* info_)
{
    prim=prim_; rhs=rhs_; info=info_;

    icM=info->icMax();
    n=icM[0]*icM[1]*icM[2];
    nPrim=info->nPrim();
    nCons=info->nCons();

    // 逐向网格间距（各向同性时与 info->interval 同值；支持各向异性均匀网格）
    hx=(info->calZone[1]-info->calZone[0])/(info->iMax[0]-1);
    hy=(info->calZone[3]-info->calZone[2])/(info->iMax[1]-1);
    hz=(info->calZone[5]-info->calZone[4])/(info->iMax[2]-1);

    mu=1.0/info->Re;                               // 基准粘性：μ=1/Re（ρ0=U0=L0=1）
    cp=GAMMA/(GAMMA-1.0);                          // R=1 约定：c_p=γ/(γ−1)=3.5
    kappa=mu*cp/info->Pr;                          // κ=μc_p/Pr
    plExp=info->powerLawExp;                       // 幂律指数（≤0 → 常粘，路径与原完全一致）
    Tref=(info->Tref>0.0)? info->Tref : 1.0;       // 幂律参考温度（非正值回退 1.0）

    if(info->viscous)
    {
        if(info->dim==3 && nPrim==5 && nCons==5)
        {
            calMethod=&ViscousTerm::calViscous3D;
            gu.resize((size_t)3*n); gv.resize((size_t)3*n);
            gw.resize((size_t)3*n); gT.resize((size_t)3*n);
        }
        else
        {
            std::cout<<"ViscousTerm warn: only 3D Euler (nPrim=5) supported this round; viscous disabled\n";
            calMethod=&ViscousTerm::nothingHappened;
        }
    }
    else calMethod=&ViscousTerm::nothingHappened;
}

void ViscousTerm::nothingHappened(){}

void ViscousTerm::calViscous(){ (this->*calMethod)(); }

void ViscousTerm::calViscous3D()
{
    const real h[3]={hx,hy,hz};
    const int off[3]={1,icM[0],icM[0]*icM[1]};
    const int nd[3]={icM[0],icM[1],icM[2]};

    // ---- 1) 12 个梯度场：∂_d u、∂_d v、∂_d w、∂_d T（d=0,1,2） ----
    for(int d=0;d<3;d++)
    {
        visc::gradField(prim,1,gu,d,icM,h[d]);
        visc::gradField(prim,2,gv,d,icM,h[d]);
        visc::gradField(prim,3,gw,d,icM,h[d]);
        visc::gradT(prim,gT,d,icM,h[d]);
    }

    // ---- 2) 逐方向：装配 Ev_d，同一 CD6 求散度并累加 rhs（−∇·Ev） ----
    std::vector<real> Ev((size_t)5*n);
    for(int d=0;d<3;d++)
    {
        #pragma omp parallel for
        for(int idx=0;idx<n;idx++)
        {
            real u=(*prim)(idx,1), v=(*prim)(idx,2), w=(*prim)(idx,3);
            real ux=gu[0*(size_t)n+idx], uy=gu[1*(size_t)n+idx], uz=gu[2*(size_t)n+idx];
            real vx=gv[0*(size_t)n+idx], vy=gv[1*(size_t)n+idx], vz=gv[2*(size_t)n+idx];
            real wx=gw[0*(size_t)n+idx], wy=gw[1*(size_t)n+idx], wz=gw[2*(size_t)n+idx];
            real div=ux+vy+wz;
            real muP=mu, kapP=kappa;               // 物性：常粘缺省
            if(plExp>0.0)                          // 幂律：按当地温度 T=p/ρ 取物性
            {
                const real Tc=(*prim)(idx,4)/(*prim)(idx,0);
                const real fac=std::pow(Tc/Tref,plExp);
                muP=mu*fac;
                kapP=kappa*fac;
            }
            real s11=muP*(2.0*ux-(2.0/3.0)*div);   // τ_xx（含 2/3 Stokes）
            real s22=muP*(2.0*vy-(2.0/3.0)*div);   // τ_yy
            real s33=muP*(2.0*wz-(2.0/3.0)*div);   // τ_zz
            real s12=muP*(uy+vx);                  // τ_xy
            real s13=muP*(uz+wx);                  // τ_xz
            real s23=muP*(vz+wy);                  // τ_yz
            real td=gT[d*(size_t)n+idx];
            real s_d0,s_d1,s_d2;
            if(d==0)      { s_d0=s11; s_d1=s12; s_d2=s13; }
            else if(d==1) { s_d0=s12; s_d1=s22; s_d2=s23; }
            else          { s_d0=s13; s_d1=s23; s_d2=s33; }
            Ev[0*(size_t)n+idx]=0.0;               // 质量分量：粘性通量为 0
            Ev[1*(size_t)n+idx]=s_d0;
            Ev[2*(size_t)n+idx]=s_d1;
            Ev[3*(size_t)n+idx]=s_d2;
            Ev[4*(size_t)n+idx]=u*s_d0+v*s_d1+w*s_d2+kapP*td;    // u·τ + κ∂_dT
        }

        #pragma omp parallel for
        for(int k=0;k<icM[2];k++)
        for(int j=0;j<icM[1];j++)
        for(int i=0;i<icM[0];i++)
        {
            int cd[3]={i,j,k};
            int idx=i+j*icM[0]+k*icM[0]*icM[1];
            int c=cd[d];
            int m3=visc::wrap(c-3,nd[d]), m2=visc::wrap(c-2,nd[d]), m1=visc::wrap(c-1,nd[d]);
            int p1=visc::wrap(c+1,nd[d]), p2=visc::wrap(c+2,nd[d]), p3=visc::wrap(c+3,nd[d]);
            int im3=idx+(m3-c)*off[d], im2=idx+(m2-c)*off[d], im1=idx+(m1-c)*off[d];
            int ip1=idx+(p1-c)*off[d], ip2=idx+(p2-c)*off[d], ip3=idx+(p3-c)*off[d];
            for(int comp=1;comp<5;comp++)
            {
                real dEv=visc::cd6(Ev[comp*(size_t)n+im3],Ev[comp*(size_t)n+im2],Ev[comp*(size_t)n+im1],
                                   Ev[comp*(size_t)n+ip1],Ev[comp*(size_t)n+ip2],Ev[comp*(size_t)n+ip3],h[d]);
                (*rhs)(idx,comp)-=dEv;
            }
        }
    }
}
