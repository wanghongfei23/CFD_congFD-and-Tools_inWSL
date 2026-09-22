/**
 * @file symmetrytest.cpp
 * @brief 对称性测试程序实现文件
 * 
 * 该程序用于测试二维欧拉方程特征系统和通量计算的对称性。
 */

#include "eigenSystem.hpp"
#include "fluxScheme.hpp"

/**
 * @brief 主函数
 * @return 程序退出状态
 * 
 * 测试不同方向上的特征变换和通量计算的对称性。
 */
int main()
{
    //left star
    // std::array<real,4> priml={1.5,0.05,0.1,1.5};
    // std::array<real,4> primr={0.5323,0.1,1.206,0.3};

    //right star
    std::array<real,4> priml={0.18217903560635557,0.8040523419013228,1.206,0.07123495167920027};
    std::array<real,4> primr={0.52330217813773561,0.0097235893105526473,1.2060000000000002,0.294926010838852};

    std::array<real,4> prim={1.0,0.05,0.8,0.8};

    std::array<real,3> normx={1,0,0};
    std::array<real,3> normy={0,1,0};
    eigensystemEuler2D testEig=eigensystemEuler2D(priml,priml,normx);
    auto Char=testEig.primToChar(prim);
    auto primTrans=testEig.charToPrim(Char);

    auto flux=roeFlux2DSym(priml[0],primr[0],
                         priml[1],primr[1],
                         priml[2],primr[2],
                         priml[3],primr[3],normx);

    auto fluxInv=HLLCFlux2D2(primr[0],priml[0],
                         primr[1],priml[1],
                         primr[2],priml[2],
                         primr[3],priml[3],normx);
    std::vector<real> primrvec={primr[0],primr[1],primr[2],primr[3]};
    auto fluxs=fEuler2D(primrvec,normx);
    
    // auto fluxinv=HLLCFlux2D2(priml[0],primr[0],
    //                      priml[1],primr[1],
    //                      priml[2],primr[2],
    //                      priml[3],primr[3],normy);

    std::array<real,4> primlm={priml[0],priml[2],priml[1],priml[3]};
    std::array<real,4> primrm={primr[0],primr[2],primr[1],primr[3]};
    std::array<real,4> primm={prim[0],prim[2],prim[1],prim[3]};

    eigensystemEuler2D testEigm=eigensystemEuler2D(primlm,primlm,normy);

    auto Charm=testEigm.primToChar(primm);
    auto primTransm=testEigm.charToPrim(Charm);

    // auto fluxm=roeFlux2DSym(primrm[0],primlm[0],
    //                      primrm[1],primlm[1],
    //                      primrm[2],primlm[2],
    //                      primrm[3],primlm[3],normy);
    auto fluxm=roeFlux2DSym(primrm[0],primlm[0],
                         primrm[1],primlm[1],
                         -primrm[2],-primlm[2],
                         primrm[3],primlm[3],normy);

    auto fluxInvm=HLLCFlux2D2(primrm[0],primlm[0],
                         primrm[1],primlm[1],
                         primrm[2],primlm[2],
                         primrm[3],primlm[3],normy);
    std::vector<real> primrmvec={primrm[0],primrm[1],primrm[2],primrm[3]};
    auto fluxsm=fEuler2D(primrmvec,normy);


    std::cout<<std::format("0:{}  1:{}  2:{}  3:{}\n"
    ,Char[0]==-Charm[0],Char[1]==Charm[1],Char[2]==Charm[2],Char[3]==Charm[3]);

    std::cout<<std::format("0:{}  1:{}  2:{}  3:{}\n"
    ,primTrans[0]==primTransm[0],primTrans[2]==primTransm[1]
    ,primTrans[2]==primTransm[1],primTrans[3]==primTransm[3]);

    std::cout<<std::format("0:{}  1:{}  2:{}  3:{}\n"
    ,flux[0]==-fluxm[0],flux[2]==-fluxm[1]
    ,flux[2]==-fluxm[1],flux[3]==-fluxm[3]);

    std::cout<<std::format("0:{}  1:{}  2:{}  3:{}\n"
    ,fluxInv[0]==fluxInvm[0],fluxInv[2]==fluxInvm[1]
    ,fluxInv[2]==fluxInvm[1],fluxInv[3]==fluxInvm[3]);

    std::cout<<std::format("0:{}  1:{}  2:{}  3:{}\n"
    ,fluxs[0]==fluxsm[0],fluxs[2]==fluxsm[1]
    ,fluxs[2]==fluxsm[1],fluxs[3]==fluxsm[3]);

    return 0;
}