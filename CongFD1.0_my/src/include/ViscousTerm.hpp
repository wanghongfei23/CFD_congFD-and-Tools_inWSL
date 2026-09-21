/**
 * @file ViscousTerm.hpp
 * @brief 粘性项模块（任务 4）：物性、CD6 梯度、粘性通量 Ev 与散度累加
 *
 * 路线①（单元中心 + 散度，主参照 OpenCFD-SC）：整场显式计算，不经线式路径；
 * 周期回绕以模索引自含实现（块级场仅存内部格；ghost 由逐线 OneDBnd 承担）。
 * 公式家：《02-粘性实现要点清单》§2.1；实施分解见 04 实施分解 S2 至 S4。
 */
#pragma once
#include "data.hpp"
#include "info.hpp"
#include <cmath>
#include <vector>

namespace visc
{
// CD6 一维一阶导（七点对称；系数 3/4、−3/20、1/60，见 02 §2.1）
inline real cd6(real fm3,real fm2,real fm1,real fp1,real fp2,real fp3,real h)
{
    return (0.75*(fp1-fm1)-0.15*(fp2-fm2)+(1.0/60.0)*(fp3-fm3))/h;
}

// 周期模回绕索引
inline int wrap(int i,int n){ int r=i%n; return r<0? r+n : r; }

// idim 方向一阶导场（CD6、周期模回绕）；dst[d*n+idx]
inline void gradField(Data* f,int ivar,std::vector<real>& dst,int d,
                      const std::array<int,3>& icM,real h)
{
    const int off[3]={1,icM[0],icM[0]*icM[1]};
    const int nd=icM[d];
    const int n=icM[0]*icM[1]*icM[2];
    for(int k=0;k<icM[2];k++)
    for(int j=0;j<icM[1];j++)
    for(int i=0;i<icM[0];i++)
    {
        int cd[3]={i,j,k};
        int idx=i+j*icM[0]+k*icM[0]*icM[1];
        int c=cd[d];
        int m3=wrap(c-3,nd), m2=wrap(c-2,nd), m1=wrap(c-1,nd);
        int p1=wrap(c+1,nd), p2=wrap(c+2,nd), p3=wrap(c+3,nd);
        dst[(size_t)d*n+idx]=cd6((*f)(idx+(m3-c)*off[d],ivar),(*f)(idx+(m2-c)*off[d],ivar),
                                 (*f)(idx+(m1-c)*off[d],ivar),(*f)(idx+(p1-c)*off[d],ivar),
                                 (*f)(idx+(p2-c)*off[d],ivar),(*f)(idx+(p3-c)*off[d],ivar),h);
    }
}

// 温度场 T=p/ρ 的 idim 方向梯度；dst[d*n+idx]
inline void gradT(Data* f,std::vector<real>& dst,int d,
                  const std::array<int,3>& icM,real h)
{
    const int off[3]={1,icM[0],icM[0]*icM[1]};
    const int nd=icM[d];
    const int n=icM[0]*icM[1]*icM[2];
    auto T=[f](int idx){ return (*f)(idx,4)/(*f)(idx,0); };
    for(int k=0;k<icM[2];k++)
    for(int j=0;j<icM[1];j++)
    for(int i=0;i<icM[0];i++)
    {
        int cd[3]={i,j,k};
        int idx=i+j*icM[0]+k*icM[0]*icM[1];
        int c=cd[d];
        int m3=wrap(c-3,nd), m2=wrap(c-2,nd), m1=wrap(c-1,nd);
        int p1=wrap(c+1,nd), p2=wrap(c+2,nd), p3=wrap(c+3,nd);
        dst[(size_t)d*n+idx]=cd6(T(idx+(m3-c)*off[d]),T(idx+(m2-c)*off[d]),T(idx+(m1-c)*off[d]),
                                 T(idx+(p1-c)*off[d]),T(idx+(p2-c)*off[d]),T(idx+(p3-c)*off[d]),h);
    }
}
}

class ViscousTerm
{
    public:
    ViscousTerm(Data* prim_,Data* rhs_,Info* info_);
    void calViscous();                 // 粘性关 → 空操作（Euler 路径不经粘性计算）

    real getMu(){ return mu; }         // 测试用访问器
    real getKappa(){ return kappa; }
    real getCp(){ return cp; }

    private:
    void calViscous3D();
    void nothingHappened();
    void (ViscousTerm::*calMethod)()=nullptr;

    Data *prim,*rhs;
    Info *info;
    int n,nPrim,nCons;
    std::array<int,3> icM;
    real hx,hy,hz;                     // 逐向网格间距（各向异性均匀网格；同向均匀）
    real mu,kappa,cp;                  // 基准物性与比热（R=1 约定；mu=1/Re）
    real plExp=0.0,Tref=1.0;           // 幂律粘性（μ=mu·(T/Tref)^plExp，κ 随动；plExp≤0 常粘）
    std::vector<real> gu,gv,gw,gT;     // 12 个梯度场（[d*n+idx]）
};
