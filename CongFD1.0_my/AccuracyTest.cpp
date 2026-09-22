/**
 * @file AccuracyTest.cpp
 * @brief 精度测试程序实现文件
 * 
 * 该程序用于测试数值方法的精度，通过计算已知解析解的导数问题来评估数值方法的收敛阶。
 */

#include "SpaceDis.hpp"
#include <fstream>
#include <format>

/**
 * @brief 计算测试函数值
 * @param x 自变量
 * @param k 指数参数
 * @return 函数值
 */
constexpr real gkx(real x,real k)
{
    //return k*x;
    return pow(x,k)*exp(0.75*(x-1));
}

/**
 * @brief 计算测试函数导数值
 * @param x 自变量
 * @param k 指数参数
 * @return 导数值
 */
constexpr real dgkx(real x,real k)
{
    //return k;
    return exp(0.75*(-1 + x))*k*pow(x,-1 + k) + 0.75*exp(0.75*(-1 + x))*pow(x,k);
}

/**
 * @brief 计算误差
 * @param h 网格步长
 * @param k 指数参数
 * @return 包含L1、L2和Linf范数误差的数组
 */
std::array<real,3> getError(real h,real k)
{
    
    real xmax=1.0,xmin=-1.0;
    int n=(xmax-xmin)/h;
    std::shared_ptr<OneDBnd> bndL,bndR;
    bndL=std::make_shared<OneDBnd>(10,1,TYPENULL);
    bndR=std::make_shared<OneDBnd>(10,1,TYPENULL);
    std::vector<real> valueL,valueR;
    valueL.resize(10);
    valueR.resize(10);
    for (size_t i = 0; i < 10; i++)
    {
        real x;
        //Left
        x=xmin-h/2-i*h;
        valueL.at(i)=gkx(x,k);
        //Right
        x=xmax+h/2+i*h;
        valueR.at(i)=gkx(x,k);
    }
    bndL->setValue(valueL);
    bndR->setValue(valueR);
    

    Data* prim,*rhs;
    prim=new Data(n,1);


    
    for (int i = 0; i < n; i++)
    {
        real x=xmin+h/2+i*h;
        (*prim)(i,0)=gkx(x,k);
    }
    

    rhs=new Data(n,1);
    rhs->setZeros();
    Info* info=new Info;

    info->eqType=ACCURACYTEST;
    info->diffMethod=MND6;
    info->spMethod=WCNS5;
    info->interMethod=TCNS5;
    info->calZone={-1.0,1.0,0,0,0,0};
    info->iMax={n+1,2,2};
    info->dim=1;

    SpaceDis spDis(n,prim,rhs,bndL,bndR,info);

    spDis.setConstNorm({1,0,0});
    spDis.setOffset(0,1);
    spDis.setIDim(0);
    spDis.difference();

    real L1=0,L2=0,Linf=0;
    for (int i = 0; i < n; i++)
    {
        real x=xmin+h/2+i*h;
        real err=std::abs((*rhs)(i,0)-dgkx(x,k));
        L1+=err;
        L2+=err*err;
        if(err>Linf) Linf=err;
    }
    L1/=n;
    L2=std::sqrt(L2/n);
    std::array<real,3> res={L1,L2,Linf};
    return res;
}

/**
 * @brief 主函数
 * @return 程序退出状态
 */
int main()
{
    real k=1;
    int n=12;
    std::vector<real> hs={2.0/10.0,2.0/20.0,2.0/30.0,2.0/40.0,2.0/50.0
                             ,2.0/60.0,2.0/70.0,2.0/80.0,2.0/90.0,2.0/100.0,2.0/200.0,2.0/300.0,2.0/400.0,2.0/800.0};
    if(hs.size() <n) return 0;
    std::vector<real> L1(n),L2(n),Linf(n);
    // for (int i = 0; i < n; i++)
    // {
    //     hs.at(i)=pow(2,-(i+1));
    // }
    for (int i = 0; i < n; i++)
    {
        auto err=getError(hs.at(i),k);
        L1.at(i)=err[0];
        L2.at(i)=err[1];
        Linf.at(i)=err[2];
        std::cout<<std::format("interval={}    L1 error={}    L2 error={}    Linf error={} \n",hs.at(i),err[0],err[1],err[2]);
    }
    std::cout<<"\n\n";
    std::vector<real> L1orders(n-1),L2orders(n-1),Linforders(n-1);
    for (int i = 0; i < n-1; i++)
    {
        L1orders.at(i)=log(L1.at(i)/L1.at(i+1))/log(hs.at(i)/hs.at(i+1));
        L2orders.at(i)=log(L2.at(i)/L2.at(i+1))/log(hs.at(i)/hs.at(i+1));
        Linforders.at(i)=log(Linf.at(i)/Linf.at(i+1))/log(hs.at(i)/hs.at(i+1));
        std::cout<<std::format("interval={}    L1 order={}    L2 order={}    Linf order={} \n"
                  ,hs.at(i),L1orders.at(i),L2orders.at(i),Linforders.at(i));
    }


    std::fstream file("accuracyTest.dat",std::ios::out);
    file<<std::format("Variables=dx,L1,L2,Linf \n Zone I={},f=point\n",n);

    for (int i = 0; i < n; i++)
    {
        file<<std::format("{} {} {} {} \n",hs.at(i),L1.at(i),L2.at(i),Linf.at(i));
    }
    file.close();
    file.open("order.txt",std::ios::out);
    file<<std::format("h L1order L2order Linforder \n");
    for (int i = 0; i < n-1; i++)
    {
        file<<std::format("{} {} {} {} \n",hs.at(i),L1orders.at(i),L2orders.at(i),Linforders.at(i));
    }
    file.close();
    
    
    
    
    
    //

}