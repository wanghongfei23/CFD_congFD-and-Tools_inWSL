/**
 * @file fluxScheme.hpp
 * @brief 定义通量格式相关的函数和类
 * 
 * 该文件包含了各种通量计算方案的实现，如Roe通量、HLLC通量等，
 * 用于计算CFD中界面的数值通量。
 */

#pragma once

#include <array>
#include <vector>
#include <algorithm>

typedef std::array<real,2> arr2;
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
std::vector<real> roeFlux1D(real rl,real rr,real ul,real ur,real pl,real pr);

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
std::vector<real> roeFlux1D2(real rl,real rr,real ul,real ur,real pl,real pr);

//std::vector<real> HLLCFlux1D(arr2 r, arr2 u, arr2 p, arr2 H, arr2 RT);

/**
 * @brief 计算一维HLLC通量
 * @param rl 左侧密度
 * @param rr 右侧密度
 * @param ul 左侧速度
 * @param ur 右侧速度
 * @param pl 左侧压力
 * @param pr 右侧压力
 * @return 通量数组
 */
constexpr std::array<real,3> HLLCFlux1D(real rl,real rr,real ul,real ur,real pl,real pr)
{
    //reference:https://zhuanlan.zhihu.com/p/583555029
    enum{
    L,
    R
    };
    arr2 H={pl/rl*GAMMA/(GAMMA-1)+(ul*ul)/2
           ,pr/rr*GAMMA/(GAMMA-1)+(ur*ur)/2};
    arr2 RT={pl/rl
           ,pr/rr};
    std::array<real,3> res;
    real gamma=GAMMA;
    real cl=std::sqrt(gamma*RT[L]);
    real cr=std::sqrt(gamma*RT[R]);

    real uBar=(ul+ur)/2,cBar=(cl+cr)/2;
    real SL=std::min(ur-cr,ul-cl);
    real SR=std::max(ul+cl,ur+cr);

    real Sstar=((pr-pl)
                +(rl*ul*(SL-ul)-rr*ur*(SR-ur)))/
               (rl*(SL-ul)-rr*(SR-ur));
    if (SL>=0)
    {
        res[0]=rl*ul;
        res[1]=rl*ul*ul+pl;
        res[2]=rl*H[L]*ul;
    }
    else if (Sstar>=0)
    {
        real pStar=pl+rl*(SL-ul)*(Sstar-ul);
        real U[3]={rl,rl*ul,pl/(gamma-1)+rl*ul*ul/2};
        real F[3]={rl*ul,rl*ul*ul+pl,rl*H[L]*ul};
        real D[3]={0,1,Sstar};
        for (int i = 0; i < 3; i++)
        {
            res[i]=(Sstar*(SL*U[i]-F[i])+SL*pStar*D[i])/(SL-Sstar);
        }
        
    }
    else if (SR>=0)
    {
        real pStar=pr+rr*(SR-ur)*(Sstar-ur);
        real U[3]={rr,rr*ur,pr/(gamma-1)+rr*ur*ur/2};
        real F[3]={rr*ur,rr*ur*ur+pr,rr*H[R]*ur};
        real D[3]={0,1,Sstar};
        for (int i = 0; i < 3; i++)
        {
            res[i]=(Sstar*(SR*U[i]-F[i])+SR*pStar*D[i])/(SR-Sstar);
        }
    }
    else
    {
        res[0]=rr*ur;
        res[1]=rr*ur*ur+pr;
        res[2]=rr*H[R]*ur;
    }
    return res;
}

/**
 * @brief 计算二维HLLC通量
 * @param rl 左侧密度
 * @param rr 右侧密度
 * @param ul 左侧速度（x方向）
 * @param ur 右侧速度（x方向）
 * @param vl 左侧速度（y方向）
 * @param vr 右侧速度（y方向）
 * @param pl 左侧压力
 * @param pr 右侧压力
 * @param norm 单位法向量
 * @return 通量数组
 */
constexpr std::array<real,4> HLLCFlux2D(real rl,real rr,real ul,real ur,real vl,real vr,real pl,real pr,std::array<real,3> norm)
{
    enum{
    L,
    R
    };
    //reference:https://zhuanlan.zhihu.com/p/583555029
    std::array<real,4> res;
    // if(pl<0 || pr<0 || rr<0 || rl<0)
    // {
    //     std::cout<<"fluxScheme error: positivity break\n";
    // }

     arr2 H;
     H[L]=(ul*ul+vl*vl)/2+pl/rl*GAMMA/(GAMMA-1);
     H[R]=(ur*ur+vr*vr)/2+pr/rr*GAMMA/(GAMMA-1);

    arr2 Vn={(norm[1]>norm[0])?vl:ul,(norm[1]>norm[0])?vr:ur};
    real gamma=GAMMA;
    real cl=std::sqrt((gamma-1)*pl/rl*GAMMA/(GAMMA-1));
    real cr=std::sqrt((gamma-1)*pr/rr*GAMMA/(GAMMA-1));

    real uBar=(Vn[L]+Vn[R])/2,cBar=(cl+cr)/2;
    // real SL=std::min(uBar-cBar,Vn[L]-cl);
    // real SR=std::max(uBar+cBar,Vn[R]+cr);

    real SL=std::min(Vn[R]-cr,Vn[L]-cl);
    real SR=std::max(Vn[L]+cl,Vn[R]+cr);

    real Sstar=((pr-pl)+(rl*Vn[L]*(SL-Vn[L])
                -rr*Vn[R]*(SR-Vn[R])))/
               (rl*(SL-Vn[L])-rr*(SR-Vn[R]));
    if (SL>0)
        {
            res[0]=rl         *Vn[L];
            res[1]=ul*Vn[L]*rl+pl*norm[0];
            res[2]=vl*Vn[L]*rl+pl*norm[1];
            res[3]=rl*H[L]*Vn[L];
        }
    else if (Sstar>0)
    {
        real pStar=pl+rl*(SL-Vn[L])*(Sstar-Vn[L]);
        real U[4]={rl,
                rl*ul,
                rl*vl,
                rl*H[L]-pl};
        real F[4]={rl     *Vn[L]
                ,ul*Vn[L]*rl+pl*norm[0]
                ,vl*Vn[L]*rl+pl*norm[1]
                ,rl*H[L]*Vn[L]};
        real D[4]={0,norm[0],norm[1],Sstar};
        for (int i = 0; i < 4; i++)
        {
            res[i]=(SL*(Sstar*U[i]+pStar*D[i])-Sstar*F[i])/(SL-Sstar);
        }
        
    }
    else if (SR>0)
    {
        real pStar=pr+rr*(SR-Vn[R])*(Sstar-Vn[R]);
        real U[4]={rr
                ,rr*ur
                ,rr*vr
                ,rr*H[R]-pr};

        real F[4]={rr     *Vn[R]
                ,ur*Vn[R]*rr+pr*norm[0]
                ,vr*Vn[R]*rr+pr*norm[1]
                ,rr*H[R]*Vn[R]};

        real D[4]={0,norm[0],norm[1],Sstar};
        for (int i = 0; i < 4; i++)
        {
            res[i]=(SR*(Sstar*U[i]+pStar*D[i])-Sstar*F[i])/(SR-Sstar);
        }
        
        
    }
    else
    {
        res[0]=rr     *Vn[R];
        res[1]=ur*Vn[R]*rr+pr*norm[0];
        res[2]=vr*Vn[R]*rr+pr*norm[1];
        res[3]=rr*H[R]*Vn[R];
    }
    
    return res;
}

/**
 * @brief 计算二维HLLC通量（第二种实现）
 * @param rl 左侧密度
 * @param rr 右侧密度
 * @param ul 左侧速度（x方向）
 * @param ur 右侧速度（x方向）
 * @param vl 左侧速度（y方向）
 * @param vr 右侧速度（y方向）
 * @param pl 左侧压力
 * @param pr 右侧压力
 * @param norm 单位法向量
 * @return 通量数组
 */
constexpr std::array<real,4> HLLCFlux2D2(real rl,real rr,real ul,real ur,real vl,real vr,real pl,real pr,std::array<real,3> norm)
{
    enum{
    L,
    R
    };
    std::array<real,4> res;
    // if(pl<0 || pr<0 || rr<0 || rl<0)
    // {
    //     std::cout<<"fluxScheme error: positivity break\n";
    // }

     arr2 H;
     H[L]=(ul*ul+vl*vl)/2+pl/rl*GAMMA/(GAMMA-1);
     H[R]=(ur*ur+vr*vr)/2+pr/rr*GAMMA/(GAMMA-1);

    arr2 Vn={(norm[1]>norm[0])?vl:ul,(norm[1]>norm[0])?vr:ur};
    real gamma=GAMMA;
    real cl=std::sqrt(pl/rl*gamma);
    real cr=std::sqrt(pr/rr*gamma);
    
    real rBar=(rl+rr)/2,cBar=(cl+cr)/2;
    real ps=std::max(0.0,(pl+pr)/2-rBar*cBar*(Vn[R]-Vn[L])/2);
    real ql=(ps<=pl)? 1.0 : std::sqrt(1.0+(gamma+1.0)/(2.0*gamma)*(ps/pl-1.0));
    real qr=(ps<=pr)? 1.0 : std::sqrt(1.0+(gamma+1.0)/(2.0*gamma)*(ps/pr-1.0));
    real SL=Vn[L]-cl*ql,SR=Vn[R]+cr*qr;

    // real uBar=(Vn[L]+Vn[R])/2,cBar=(cl+cr)/2;
    // real SL=std::min(uBar-cBar,Vn[L]-cl);
    // real SR=std::max(uBar+cBar,Vn[R]+cr);

    // real SL=std::min(Vn[R]-cr,Vn[L]-cl);
    // real SR=std::max(Vn[L]+cl,Vn[R]+cr);

    real Sstar=((pr-pl)+
                (rl*Vn[L]*(SL-Vn[L])-rr*Vn[R]*(SR-Vn[R])))/
               (rl*(SL-Vn[L])-rr*(SR-Vn[R]));
    if (SL>0)
    {
        res[0]=rl         *Vn[L];
        res[1]=rl*ul*Vn[L]+pl*norm[0];
        res[2]=rl*vl*Vn[L]+pl*norm[1];
        res[3]=rl*H[L]*Vn[L];
    }
    else if (Sstar>0)
    {
        real U[4]={rl,
                rl*ul,
                rl*vl,
                rl*H[L]-pl};
        real F[4]={rl     *Vn[L]
                ,rl*ul*Vn[L]+pl*norm[0]
                ,rl*vl*Vn[L]+pl*norm[1]
                ,rl*H[L]*Vn[L]};
        real Usf=(SL-Vn[L])/(SL-Sstar);
        real Ustar[4]={Usf*rl,
                       Usf*(Sstar*norm[0]+ul*norm[1])*rl,
                       Usf*(Sstar*norm[1]+vl*norm[0])*rl,
                       Usf*(U[3]+(Sstar-Vn[L])*(rl*Sstar+pl/(SL-Vn[L])))
                      };
        for (int i = 0; i < 4; i++)
        {
            res[i]=F[i]+SL*(Ustar[i]-U[i]);
        }
    }
    else if (Sstar==0)
    {
        {
            real U[4]={rl,
                    rl*ul,
                    rl*vl,
                    rl*H[L]-pl};
            real F[4]={rl     *Vn[L]
                    ,rl*ul*Vn[L]+pl*norm[0]
                    ,rl*vl*Vn[L]+pl*norm[1]
                    ,rl*H[L]*Vn[L]};
            real Usf=(SL-Vn[L])/(SL-Sstar);
            real Ustar[4]={Usf*rl,
                        Usf*(Sstar*norm[0]+ul*norm[1])*rl,
                        Usf*(Sstar*norm[1]+vl*norm[0])*rl,
                        Usf*(U[3]+(Sstar-Vn[L])*(rl*Sstar+pl/(SL-Vn[L])))
                        };
            for (int i = 0; i < 4; i++)
            {
                res[i]=(F[i]+SL*(Ustar[i]-U[i]))/2;
            }
        }

        {
            real U[4]={rr
                    ,rr*ur
                    ,rr*vr
                    ,rr*H[R]-pr};

            real F[4]={rr     *Vn[R]
                    ,rr*ur*Vn[R]+pr*norm[0]
                    ,rr*vr*Vn[R]+pr*norm[1]
                    ,rr*H[R]*Vn[R]};
            real Usf=(SR-Vn[R])/(SR-Sstar);
            real Ustar[4]={Usf*rr,
                        Usf*(Sstar*norm[0]+ur*norm[1])*rr,
                        Usf*(Sstar*norm[1]+vr*norm[0])*rr,
                        Usf*(U[3]+(Sstar-Vn[R])*(rr*Sstar+pr/(SR-Vn[R])))
                        };
            for (int i = 0; i < 4; i++)
            {
                res[i]+=(F[i]+SR*(Ustar[i]-U[i]))/2;
            }
        }
        
    }
    else if (SR>0)
    {
        real U[4]={rr
                ,rr*ur
                ,rr*vr
                ,rr*H[R]-pr};

        real F[4]={rr     *Vn[R]
                ,rr*ur*Vn[R]+pr*norm[0]
                ,rr*vr*Vn[R]+pr*norm[1]
                ,rr*H[R]*Vn[R]};
        real Usf=(SR-Vn[R])/(SR-Sstar);
        real Ustar[4]={Usf*rr,
                       Usf*(Sstar*norm[0]+ur*norm[1])*rr,
                       Usf*(Sstar*norm[1]+vr*norm[0])*rr,
                       Usf*(U[3]+(Sstar-Vn[R])*(rr*Sstar+pr/(SR-Vn[R])))
                      };
        for (int i = 0; i < 4; i++)
        {
            res[i]=F[i]+SR*(Ustar[i]-U[i]);
        }
    }
    else
    {
        res[0]=rr     *Vn[R];
        res[1]=rr*ur*Vn[R]+pr*norm[0];
        res[2]=rr*vr*Vn[R]+pr*norm[1];
        res[3]=rr*H[R]*Vn[R];
    }
    
    return res;
}

/**
 * @brief 计算二维Roe通量
 * @param rl 左侧密度
 * @param rr 右侧密度
 * @param ul 左侧速度（x方向）
 * @param ur 右侧速度（x方向）
 * @param vl 左侧速度（y方向）
 * @param vr 右侧速度（y方向）
 * @param pl 左侧压力
 * @param pr 右侧压力
 * @param norm 单位法向量
 * @return 通量数组
 */
constexpr std::array<real,4> roeFlux2D(real rl,real rr,real ul,real ur,real vl,real vr,real pl,real pr,std::array<real,3> norm)
{
    enum{
    L,
    R
    };
    std::array<real,4> FcL,FcR;
    std::array<real,4> result;
    if(pl<0 || pr<0)
    {
        std::cout<<"fluxScheme error: Pressure positivity break\n";
    }
    if(rr<0 || rl<0)
    {
        std::cout<<"fluxScheme error: Density positivity break\n";
    }


    arr2 H;
    H[L]=(ul*ul+vl*vl)/2+pl/rl*GAMMA/(GAMMA-1);
    H[R]=(ur*ur+vr*vr)/2+pr/rr*GAMMA/(GAMMA-1);
    
    arr2 Vn={ul*norm[0]+vl*norm[1],ur*norm[0]+vr*norm[1]};
    FcL[0]=rl*Vn[L];
    FcL[1]=ul*Vn[L]*rl+pl*norm[0];
    FcL[2]=vl*Vn[L]*rl+pl*norm[1];
    FcL[3]=rl*H[L]*Vn[L];

    FcR[0]=rr*Vn[R];
    FcR[1]=ur*Vn[R]*rr+pr*norm[0];
    FcR[2]=vr*Vn[R]*rr+pr*norm[1];
    FcR[3]=rr*H[R]*Vn[R];

    double rhoAvg,uAvg,vAvg,HAvg,cAvg,VnAvg,q_2Avg,coef1,coef2;
    coef1=std::sqrt(rl);
    coef2=std::sqrt(rr);
    real divisor=1.0/(std::sqrt(rl)+std::sqrt(rr));
    rhoAvg=std::sqrt(rl*rr);
    uAvg=(coef1*ul+coef2*ur)*divisor;
    vAvg=(coef1*vl+coef2*vr)*divisor;
    HAvg=(coef1*H[L]+coef2*H[R])*divisor;
    q_2Avg=(uAvg*uAvg+vAvg*vAvg)/2;

    cAvg=std::sqrt((GAMMA-1)*(HAvg-q_2Avg));
    if((HAvg-q_2Avg)<0)
    {
        cAvg=std::sqrt(0.1);
    }
    VnAvg=uAvg*norm[0]+vAvg*norm[1];

    double lambda[3]={std::abs(VnAvg-cAvg),std::abs(VnAvg),std::abs(VnAvg+cAvg)};
    real eps=0.05*(std::abs(VnAvg)+cAvg);
    for(int i=0;i<3;i++)
    {
        if(lambda[i]<eps)
        {
            lambda[i]=(lambda[i]*lambda[i]+eps*eps)/(2.0*eps);
        }
    }

    double deltaP=pr-pl,deltaVn=Vn[R]-Vn[L],deltaU=ur-ul,deltaV=vr-vl,
           deltaRho=rr-rl,coef;
    double FDispassion[4];
    coef1=(deltaP-rhoAvg*cAvg*deltaVn)/(2*cAvg*cAvg);
    FDispassion[0]=lambda[0]*coef1*1;
    FDispassion[1]=lambda[0]*coef1*(uAvg-cAvg*norm[0]);
    FDispassion[2]=lambda[0]*coef1*(vAvg-cAvg*norm[1]);
    FDispassion[3]=lambda[0]*coef1*(HAvg-cAvg*VnAvg);

    coef2=deltaRho-deltaP/(cAvg*cAvg);
    FDispassion[0]+=lambda[1]*(coef2*1.0       +  rhoAvg*0.0);
    FDispassion[1]+=lambda[1]*(coef2*uAvg      +  rhoAvg*(deltaU-deltaVn*norm[0]));
    FDispassion[2]+=lambda[1]*(coef2*vAvg      +  rhoAvg*(deltaV-deltaVn*norm[1]));
    FDispassion[3]+=lambda[1]*(coef2*q_2Avg+  rhoAvg*(uAvg*deltaU+vAvg*deltaV-VnAvg*deltaVn));

    coef=(deltaP+rhoAvg*cAvg*deltaVn)/(2.0*cAvg*cAvg);
    FDispassion[0]+=lambda[2]*(coef*1);
    FDispassion[1]+=lambda[2]*(coef*(uAvg+cAvg*norm[0]));
    FDispassion[2]+=lambda[2]*(coef*(vAvg+cAvg*norm[1]));
    FDispassion[3]+=lambda[2]*(coef*(HAvg+cAvg*VnAvg));

    for(int i=0;i<4;i++)
    {
        result[i]=0.5*(FcL[i]+FcR[i]-FDispassion[i]);
    }
    return result;
}

/**
 * @brief 计算二维Roe通量（对称版本）
 * @param rl 左侧密度
 * @param rr 右侧密度
 * @param ul 左侧速度（x方向）
 * @param ur 右侧速度（x方向）
 * @param vl 左侧速度（y方向）
 * @param vr 右侧速度（y方向）
 * @param pl 左侧压力
 * @param pr 右侧压力
 * @param norm 单位法向量
 * @return 通量数组
 */
constexpr std::array<real,4> roeFlux2DSym(real rl,real rr,real ul,real ur,real vl,real vr,real pl,real pr,std::array<real,3> norm)
{
    enum{
    L,
    R
    };
    std::array<real,4> FcL,FcR;
    std::array<real,4> result;
    if(pl<0 || pr<0)
    {
        std::cout<<"fluxScheme error: Pressure positivity break\n";
    }
    if(rr<0 || rl<0)
    {
        std::cout<<"fluxScheme error: Density positivity break\n";
    }


    arr2 H;
    H[L]=(ul*ul+vl*vl)/2+pl/rl*GAMMA/(GAMMA-1);
    H[R]=(ur*ur+vr*vr)/2+pr/rr*GAMMA/(GAMMA-1);
    
    arr2 Vn={(norm[0]>norm[1])? ul:vl,(norm[0]>norm[1])? ur:vr};
    FcL[0]=rl*Vn[L];
    FcL[1]=rl*ul*Vn[L]+pl*norm[0];
    FcL[2]=rl*vl*Vn[L]+pl*norm[1];
    FcL[3]=rl*H[L]*Vn[L];

    FcR[0]=rr*Vn[R];
    FcR[1]=rr*ur*Vn[R]+pr*norm[0];
    FcR[2]=rr*vr*Vn[R]+pr*norm[1];
    FcR[3]=rr*H[R]*Vn[R];

    double rhoAvg,uAvg,vAvg,HAvg,cAvg,VnAvg,q_2Avg,coef1,coef2;
    coef1=std::sqrt(rl);
    coef2=std::sqrt(rr);
    real divisor=1.0/(std::sqrt(rl)+std::sqrt(rr));
    rhoAvg=std::sqrt(rl*rr);
    uAvg=(coef1*ul+coef2*ur)*divisor;
    vAvg=(coef1*vl+coef2*vr)*divisor;
    HAvg=(coef1*H[L]+coef2*H[R])*divisor;
    q_2Avg=(uAvg*uAvg+vAvg*vAvg)/2;

    cAvg=std::sqrt((GAMMA-1)*(HAvg-q_2Avg));
    // if((HAvg-q_2Avg)<0)
    // {
    //     cAvg=std::sqrt(0.1);
    // }
    VnAvg=(norm[0]>norm[1])? uAvg:vAvg;

    double lambda[3]={std::abs(VnAvg-cAvg),std::abs(VnAvg),std::abs(VnAvg+cAvg)};
    // real eps=0.1*(std::abs(VnAvg)+cAvg);
    // for(int i=0;i<3;i++)
    // {
    //     if(lambda[i]<eps)
    //     {
    //         lambda[i]=(lambda[i]*lambda[i]+eps*eps)/(2.0*eps);
    //     }
    // }
    real eps=0.05*(std::abs(VnAvg)+cAvg);
    for(int i=0;i<3;i++)
    {
        if(lambda[i]<eps)
        {
            lambda[i]=(lambda[i]*lambda[i]+eps*eps)/(2.0*eps);
        }
    }

    double deltaP=pr-pl,deltaVn=Vn[R]-Vn[L],deltaU=ur-ul,deltaV=vr-vl,
           deltaRho=rr-rl,coef;
    double FDispassion[4];
    coef1=(deltaP-rhoAvg*cAvg*deltaVn)/(2.0*cAvg*cAvg);
    FDispassion[0]=lambda[0]*(coef1*1);
    FDispassion[1]=lambda[0]*(coef1*(uAvg-cAvg*norm[0]));
    FDispassion[2]=lambda[0]*(coef1*(vAvg-cAvg*norm[1]));
    FDispassion[3]=lambda[0]*(coef1*(HAvg-cAvg*VnAvg));

    coef=(deltaP+rhoAvg*cAvg*deltaVn)/(2.0*cAvg*cAvg);
    FDispassion[0]+=lambda[2]*(coef*1);
    FDispassion[1]+=lambda[2]*(coef*(uAvg+cAvg*norm[0]));
    FDispassion[2]+=lambda[2]*(coef*(vAvg+cAvg*norm[1]));
    FDispassion[3]+=lambda[2]*(coef*(HAvg+cAvg*VnAvg));

    coef2=deltaRho-deltaP/(cAvg*cAvg);
    FDispassion[0]+=lambda[1]*(coef2*1.0       +  rhoAvg*0.0);
    FDispassion[1]+=lambda[1]*(coef2*uAvg      +  rhoAvg*((norm[1]>norm[0])? deltaU:0));
    FDispassion[2]+=lambda[1]*(coef2*vAvg      +  rhoAvg*((norm[1]>norm[0])? 0:deltaV));
    FDispassion[3]+=lambda[1]*(coef2*q_2Avg+  rhoAvg*((norm[1]>norm[0])? uAvg*deltaU:vAvg*deltaV));



    for(int i=0;i<4;i++)
    {
        result[i]=(FcL[i]+FcR[i]-FDispassion[i])/2;
    }
    return result;
}

/**
 * @brief 计算二维Euler方程通量
 * @param q 状态向量
 * @param norm 单位法向量
 * @return 通量向量
 */
constexpr std::vector<real> fEuler2D(const std::vector<real>& q ,std::array<real,3> norm)
{
    real r=q[0],u=q[1],v=q[2],p=q[3];
    real Vn=norm[0]>norm[1]? u:v;//(u*norm[0]+v*norm[1]);
    //real H=GAMMA/(GAMMA-1)*p/r+(u*u+v*v)/2;
    std::vector<real> res={
        r*Vn,
        // u*Vn*r+p*norm[0],
        // v*Vn*r+p*norm[1],
        r*u*Vn+(norm[0]>norm[1]? p:0),
        r*v*Vn+(norm[0]>norm[1]? 0:p),
        ((u*u+v*v)/2*r+GAMMA/(GAMMA-1)*p)*Vn
    };
    return res;
}

/**
 * @brief 计算一维Euler方程通量
 * @param q 状态向量
 * @param norm 单位法向量
 * @return 通量向量
 */
constexpr std::vector<real> fEuler1D(const std::vector<real>& q ,std::array<real,3> norm)
{
    real r=q[0],u=q[1],p=q[2];
    //real H=GAMMA/(GAMMA-1)*p/r+(u*u)/2;
    std::vector<real> res={
        r*u,
        r*u*u+p,
        (GAMMA/(GAMMA-1)*p+r*(u*u)/2)*u
    };
    return res;
}

/**
 * @brief 默认通量计算函数
 * @param q 状态向量
 * @param norm 单位法向量
 * @return 通量向量
 */
constexpr std::vector<real> fDefault(const std::vector<real>& q ,std::array<real,3> norm)
{
    return {q[0]};
}

/* ================= 三维（任务 3 新增；与上面 2D 版一一对应） ================= */

/**
 * @brief 计算三维Euler方程物理通量
 * @param q 原始变量状态向量 [rho, u, v, w, p]
 * @param norm 单位法向量（本码仅使用轴向法向 {1,0,0}/{0,1,0}/{0,0,1}）
 * @return 通量向量（5 分量）
 */
constexpr std::vector<real> fEuler3D(const std::vector<real>& q ,std::array<real,3> norm)
{
    real r=q[0],u=q[1],v=q[2],w=q[3],p=q[4];
    // 轴向判定：把 2D 版的 norm[0]>norm[1] 技巧扩展为三选一
    const bool isX = (norm[0]>norm[1]) && (norm[0]>norm[2]);
    const bool isY = (!isX) && (norm[1]>norm[2]);
    const bool isZ = !(isX||isY);
    real Vn = isX? u : (isY? v : w);
    std::vector<real> res={
        r*Vn,
        r*u*Vn+(isX? p:0.0),
        r*v*Vn+(isY? p:0.0),
        r*w*Vn+(isZ? p:0.0),
        ((u*u+v*v+w*w)/2*r+GAMMA/(GAMMA-1)*p)*Vn
    };
    return res;
}

/**
 * @brief 计算三维Euler方程Roe通量（对称版本）
 *
 * 与 roeFlux2DSym 结构完全一致，仅扩展为 5 分量：
 *   守恒量与 Roe 平均增加 z 向分量（w）；剪切波贡献按法向轴向
 *   分配到两个切向动量行（法向行为 0）。
 */
constexpr std::array<real,5> roeFlux3DSym(real rl,real rr,real ul,real ur,real vl,real vr,real wl,real wr,real pl,real pr,std::array<real,3> norm)
{
    enum{
    L,
    R
    };
    std::array<real,5> FcL,FcR;
    std::array<real,5> result;
    if(pl<0 || pr<0)
    {
        std::cout<<"fluxScheme error: Pressure positivity break\n";
    }
    if(rr<0 || rl<0)
    {
        std::cout<<"fluxScheme error: Density positivity break\n";
    }

    // 轴向判定与 fEuler3D 相同
    const bool isX = (norm[0]>norm[1]) && (norm[0]>norm[2]);
    const bool isY = (!isX) && (norm[1]>norm[2]);
    const bool isZ = !(isX||isY);

    // 左右物理通量（H 与 Vn 均含 w 分量）
    arr2 H;
    H[L]=(ul*ul+vl*vl+wl*wl)/2+pl/rl*GAMMA/(GAMMA-1);
    H[R]=(ur*ur+vr*vr+wr*wr)/2+pr/rr*GAMMA/(GAMMA-1);

    arr2 Vn={ isX? ul : (isY? vl : wl), isX? ur : (isY? vr : wr)};
    FcL[0]=rl*Vn[L];
    FcL[1]=rl*ul*Vn[L]+(isX? pl:0.0);
    FcL[2]=rl*vl*Vn[L]+(isY? pl:0.0);
    FcL[3]=rl*wl*Vn[L]+(isZ? pl:0.0);
    FcL[4]=rl*H[L]*Vn[L];

    FcR[0]=rr*Vn[R];
    FcR[1]=rr*ur*Vn[R]+(isX? pr:0.0);
    FcR[2]=rr*vr*Vn[R]+(isY? pr:0.0);
    FcR[3]=rr*wr*Vn[R]+(isZ? pr:0.0);
    FcR[4]=rr*H[R]*Vn[R];

    // Roe 平均（sqrt(rho) 加权）
    double rhoAvg,uAvg,vAvg,wAvg,HAvg,cAvg,VnAvg,q_2Avg,coef1,coef2;
    coef1=std::sqrt(rl);
    coef2=std::sqrt(rr);
    real divisor=1.0/(std::sqrt(rl)+std::sqrt(rr));
    rhoAvg=std::sqrt(rl*rr);
    uAvg=(coef1*ul+coef2*ur)*divisor;
    vAvg=(coef1*vl+coef2*vr)*divisor;
    wAvg=(coef1*wl+coef2*wr)*divisor;
    HAvg=(coef1*H[L]+coef2*H[R])*divisor;
    q_2Avg=(uAvg*uAvg+vAvg*vAvg+wAvg*wAvg)/2;

    cAvg=std::sqrt((GAMMA-1)*(HAvg-q_2Avg));
    VnAvg= isX? uAvg : (isY? vAvg : wAvg);

    // 特征值 |Vn±c|、|Vn|（含熵修正，与 2D 版一致）
    double lambda[3]={std::abs(VnAvg-cAvg),std::abs(VnAvg),std::abs(VnAvg+cAvg)};
    real eps=0.05*(std::abs(VnAvg)+cAvg);
    for(int i=0;i<3;i++)
    {
        if(lambda[i]<eps)
        {
            lambda[i]=(lambda[i]*lambda[i]+eps*eps)/(2.0*eps);
        }
    }

    // 数值耗散项（5 分量）
    double deltaP=pr-pl,deltaVn=Vn[R]-Vn[L],deltaU=ur-ul,deltaV=vr-vl,deltaW=wr-wl,
           deltaRho=rr-rl,coef;
    double FDispassion[5];
    coef1=(deltaP-rhoAvg*cAvg*deltaVn)/(2.0*cAvg*cAvg);
    FDispassion[0]=lambda[0]*(coef1*1);
    FDispassion[1]=lambda[0]*(coef1*(uAvg-cAvg*norm[0]));
    FDispassion[2]=lambda[0]*(coef1*(vAvg-cAvg*norm[1]));
    FDispassion[3]=lambda[0]*(coef1*(wAvg-cAvg*norm[2]));
    FDispassion[4]=lambda[0]*(coef1*(HAvg-cAvg*VnAvg));

    coef=(deltaP+rhoAvg*cAvg*deltaVn)/(2.0*cAvg*cAvg);
    FDispassion[0]+=lambda[2]*(coef*1);
    FDispassion[1]+=lambda[2]*(coef*(uAvg+cAvg*norm[0]));
    FDispassion[2]+=lambda[2]*(coef*(vAvg+cAvg*norm[1]));
    FDispassion[3]+=lambda[2]*(coef*(wAvg+cAvg*norm[2]));
    FDispassion[4]+=lambda[2]*(coef*(HAvg+cAvg*VnAvg));

    // 剪切波项：法向动量行为 0，两个切向动量行与能量行承接速度跳变
    coef2=deltaRho-deltaP/(cAvg*cAvg);
    FDispassion[0]+=lambda[1]*(coef2*1.0 + rhoAvg*0.0);
    FDispassion[1]+=lambda[1]*(coef2*uAvg + rhoAvg*(isX? 0.0 : deltaU));
    FDispassion[2]+=lambda[1]*(coef2*vAvg + rhoAvg*(isY? 0.0 : deltaV));
    FDispassion[3]+=lambda[1]*(coef2*wAvg + rhoAvg*(isZ? 0.0 : deltaW));
    FDispassion[4]+=lambda[1]*(coef2*q_2Avg
                              + rhoAvg*(isX? (vAvg*deltaV+wAvg*deltaW)
                                      : (isY? (uAvg*deltaU+wAvg*deltaW)
                                            : (uAvg*deltaU+vAvg*deltaV))));

    for(int i=0;i<5;i++)
    {
        result[i]=(FcL[i]+FcR[i]-FDispassion[i])/2;
    }
    return result;
}