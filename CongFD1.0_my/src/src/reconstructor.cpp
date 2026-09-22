/**
 * @file reconstructor.cpp
 * @brief 重构器类的实现文件
 */

#include "reconstructor.hpp"
#include "eigenSystem.hpp"
#include "interScheme.hpp"

/**
 * @brief 构造函数
 * @param n_ 网格点数
 * @param vals_ 原始数据指针
 * @param valsR_ 重构后数据指针
 * @param bndL_ 左边界数据指针
 * @param bndR_ 右边界数据指针
 * @param info_ Info对象指针
 * @param i0_ 起始索引
 * @param offset_ 偏移量
 */
Reconstructor::Reconstructor(int n_,Data* vals_,std::shared_ptr<Data> valsR_,std::shared_ptr<OneDBnd> bndL_,std::shared_ptr<OneDBnd> bndR_,Info* info_,int i0_, int offset_)
{
    vals=vals_;
    bndL=bndL_;
    bndR=bndR_;
    info=info_;
    i0=i0_;
    offset=offset_;
    n=n_;
    valsR=valsR_;
    nval=vals->getNVar();
}

/**
 * @brief 访问指定位置的数据
 * @param i 索引
 * @param ivar 变量索引
 * @return 数据引用
 */
real& Reconstructor::at(int i,int ivar)
{
    if (i<0) 
    {
        return (*bndL)(-(i+1),ivar);
    }
    if (i>=n) 
    {
        return (*bndR)(i-n,ivar);
    }
    return (*vals)[(i0+offset*i)*nval+ivar];
}

/**
 * @brief 实现特征重构方法
 * 
 * 对一维欧拉方程进行特征重构，使用WENO5方法。
 */
void ReconEigen1DEuler::recon()
{
    int nvalr=valsR->getN();
    int dn=(nvalr-n)/2;
    assert((nvalr-n)%2==0);
    int nVar=vals->getNVar();

    for (int i=0-dn;i<n+dn;i++)
    {
        std::array<real,3> primL,primR;
        memcpy(&primL[0],&at(i-1,0),nVar*sizeof(real));
        memcpy(&primR[0],&at(i,0),nVar*sizeof(real));
        eigensystemEuler1D eig=eigensystemEuler1D(primL,primR);
        std::array<real,5> q1L,q2L,q3L,q1R,q2R,q3R;
        for(int j=i-3;j<i+3;j++)
        {
            enum{R,U,P};
            auto charTemp=eig.primToChar({at(j,R),at(j,U),at(j,P)});

            int iLocal=j-i+3;
            if(iLocal<5){
            q1L[iLocal]=charTemp[0];
            q2L[iLocal]=charTemp[1];
            q3L[iLocal]=charTemp[2];}

            iLocal=i+2-j;
            if(iLocal<5){
            q1R[iLocal]=charTemp[0];
            q2R[iLocal]=charTemp[1];
            q3R[iLocal]=charTemp[2];}
        }
        // auto start = std::chrono::steady_clock::now(); 
        auto Q1LL=weno5_JSchen(q1L);
        auto Q1RR=weno5_JSchen(q1R);
        auto Q2LL=weno5_JSchen(q2L);
        auto Q2RR=weno5_JSchen(q2R);
        auto Q3LL=weno5_JSchen(q3L);
        auto Q3RR=weno5_JSchen(q3R);
        // auto stop = std::chrono::steady_clock::now(); 
        // auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count(); 
        // timep+=duration;

        auto resTempL=eig.charToPrim({Q1LL,Q2LL,Q3LL});
        auto resTempR=eig.charToPrim({Q1RR,Q2RR,Q3RR});
        // return {resTempL[0],resTempL[1],resTempL[2],resTempR[0],resTempR[1],resTempR[2]};
        // (*valsR)((i+dn)*2,0)
        memcpy(&(*valsR)((i+dn)*2,0),&resTempL[0],nVar*sizeof(real));
    }
    
}