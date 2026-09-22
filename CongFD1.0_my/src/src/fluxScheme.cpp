/**
 * @file fluxScheme.cpp
 * @brief 通量格式相关函数的实现文件
 */

#include "macro.hpp"
#include "fluxScheme.hpp"
#include <algorithm>


enum{
    L,
    R
};

/**
 * @brief 计算一维Roe通量
 * @param rl 左侧密度
 * @param rr 右侧密度
 * @param ul 左侧速度
 * @param ur 右侧速度
 * @param pl 左侧压力
 * @param pr 右侧压力
 * @return 通量向量
 */
std::vector<real> roeFlux1D(real rl,real rr,real ul,real ur,real pl,real pr)
{
    std::vector<real> res;
    res.resize(3);
    
    real gamma=GAMMA;
    arr2 H={pl/rl*GAMMA/(GAMMA-1)+(ul*ul)/2
        ,pr/rr*GAMMA/(GAMMA-1)+(ur*ur)/2};

    real rBar=std::sqrt(rl*rr);
    real uBar=(ul*std::sqrt(rl)+ur*std::sqrt(rr))/(std::sqrt(rl)+std::sqrt(rr));
    real HBar=(H[L]*std::sqrt(rl)+H[R]*std::sqrt(rr))/(std::sqrt(rl)+std::sqrt(rr));
    real cBar=std::sqrt((gamma-1)*(HBar-uBar*uBar/2));
    real cBar2=(gamma-1)*(HBar-uBar*uBar/2);

    real dr=rr-rl;
    real du=ur-ul;
    real dp=pr-pl;

    real K1[3]={1,uBar,uBar*uBar/2};
    real K2[3]={1,uBar-cBar,HBar-uBar*cBar};
    real K3[3]={1,uBar+cBar,HBar+uBar*cBar};

    real alpha[3]={
        dr-dp/cBar2,
        1.0/(2.0*cBar2)*(dp-rBar*cBar*du),
        1.0/(2.0*cBar2)*(dp+rBar*cBar*du)
    };
    real lambda[3]={
        std::abs(uBar),std::abs(uBar-cBar),std::abs(uBar+cBar)
    };
    /*
    real eps=0.1*(std::abs(uBar)+std::abs(uBar)+cBar);
    for(int i=0;i<3;i++)
    {
        if(lambda[i]<eps)
        {
            lambda[i]=(lambda[i]*lambda[i]+eps*eps)/(2.0*eps);
        }
    }*/

    real FL[3]={
        rl*ul,
        rl*ul*ul+pl,
        rl*H[L]*ul
    };
    real FR[3]={
        rr*ur,
        rr*ur*ur+pr,
        rr*H[R]*ur
    };

    for (int i = 0; i < 3; i++)
    {
        res[i]=(FL[i]+FR[i]
              -lambda[0]*alpha[0]*K1[i]
              -lambda[1]*alpha[1]*K2[i]
              -lambda[2]*alpha[2]*K3[i]
                )/2.0;
    }
    return res;


}

/**
 * @brief 计算一维Roe通量（第二种实现）
 * @param rl 左侧密度
 * @param rr 右侧密度
 * @param ul 左侧速度
 * @param ur 右侧速度
 * @param pl 左侧压力
 * @param pr 右侧压力
 * @return 通量向量
 */
std::vector<real> roeFlux1D2(real rl,real rr,real ul,real ur,real pl,real pr)
{
    //reference: https://blog.csdn.net/Tankrun1997/article/details/132743487
    std::vector<real> res;
    res.resize(3);
    
    real gamma=GAMMA;

    arr2 H={pl/rl*GAMMA/(GAMMA-1)+(ul*ul)/2
           ,pr/rr*GAMMA/(GAMMA-1)+(ur*ur)/2};

    real rBar=std::sqrt(rl*rr);
    real uBar=(ul*std::sqrt(rl)+ur*std::sqrt(rr))/(std::sqrt(rl)+std::sqrt(rr));
    real HBar=(H[L]*std::sqrt(rl)+H[R]*std::sqrt(rr))/(std::sqrt(rl)+std::sqrt(rr));
    real cBar=std::sqrt((gamma-1)*(HBar-uBar*uBar/2));
    real cBar2=(gamma-1)*(HBar-uBar*uBar/2);

    real dr=rr-rl;
    real du=ur-ul;
    real dp=pr-pl;

    real K1[3]={1,uBar,uBar*uBar/2};
    real K2[3]={1,uBar-cBar,HBar-uBar*cBar};
    real K3[3]={1,uBar+cBar,HBar+uBar*cBar};


    real lambda[3]={
        std::abs(uBar),std::abs(uBar+cBar),std::abs(uBar-cBar)
    };
    
    real eps=0.1*(std::abs(uBar)+cBar);
    for(int i=0;i<3;i++)
    {
        if(lambda[i]<eps)
        {
            lambda[i]=(lambda[i]*lambda[i]+eps*eps)/(2.0*eps);
        }
    }

    real alpha[5]={
        lambda[0]*(dr-dp/cBar2),
        lambda[1]/(2.0*cBar2)*(dp+rBar*cBar*du),
        lambda[2]/(2.0*cBar2)*(dp-rBar*cBar*du)
    };
    alpha[3]=alpha[0]+alpha[1]+alpha[2];
    alpha[4]=cBar*(alpha[1]-alpha[2]);

    real FL[3]={
        rl*ul,
        rl*ul*ul+pl,
        rl*H[L]*ul
    };
    real FR[3]={
        rr*ur,
        rr*ur*ur+pr,
        rr*H[R]*ur
    };

    for (int i = 0; i < 3; i++)
    {
        res[i]=(FL[i]+FR[i])/2.0;
    }
    res[0]-=0.5*alpha[3];
    res[1]-=0.5*(alpha[3]*uBar+alpha[4]);
    res[2]-=0.5*(HBar*alpha[3]+uBar*alpha[4]-cBar2*alpha[0]/(gamma-1));



    return res;


}