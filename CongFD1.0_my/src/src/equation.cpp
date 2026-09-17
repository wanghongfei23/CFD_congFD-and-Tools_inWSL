/**
 * @file equation.cpp
 * @brief 方程类的实现文件
 */

#include "equation.hpp"

/**
 * @brief 将守恒变量转换为原始变量
 * 
 * 根据方程类型调用相应的转换函数，支持线性对流方程、Burgers方程和欧拉方程。
 */
void Equation::consToPrim()
{
    switch (type)
    {
    case LINEARCONV1D:  //一维线性对流方程
    case BURGERS1D:  //一维burgers方程
        std::copy(cons->begin(),cons->end(),prim->begin());
        break;
    case EULER:  //欧拉方程
        {
            if (dim==1) consToPrimEuler1D();
            else if(dim==2)
            {
                consToPrimEuler2D();
            }
            else if(dim==3)
            {
                consToPrimEuler3D();
            }
        }
        break;
    
    default:
        break;
    }
}

/**
 * @brief 一维欧拉方程的守恒变量到原始变量转换
 * 
 * 将一维欧拉方程的守恒变量（rho, rho*u, rho*E）转换为原始变量（rho, u, p）。
 */
void Equation::consToPrimEuler1D()
{
    if(nCons!=3,nPrim!=3)
    {
        std::cout<<"Equation error: Euler 1d equation variable number error \n";
    }

    for (int i = 0; i < n; i++)
    {
        real r=(*cons)(i,0);
        real ru=(*cons)(i,1);
        real rE=(*cons)(i,2);
        real u=ru/r;
        real E=rE/r;
        real e=-u*u/2+E;
        real gamma=GAMMA;
        real RT=(gamma-1)*e;
        real p=r*RT;
        (*prim)(i,0)=r;
        (*prim)(i,1)=u;
        (*prim)(i,2)=p;
    }
}

/**
 * @brief 二维欧拉方程的守恒变量到原始变量转换
 * 
 * 将二维欧拉方程的守恒变量（rho, rho*u, rho*v, rho*E）转换为原始变量（rho, u, v, p）。
 */
void Equation::consToPrimEuler2D()
{
    if(nCons!=3,nPrim!=4)
    {
        std::cout<<"Equation error: Euler 1d equation variable number error \n";
    }

    for (int i = 0; i < n; i++)
    {
        real r=(*cons)(i,0);
        real ru=(*cons)(i,1);
        real rv=(*cons)(i,2);
        real rE=(*cons)(i,3);
        real u=ru/r;
        real v=rv/r;
        real E=rE/r;
        real q2=(u*u+v*v)/2;
        real e=-q2+E;
        real gamma=GAMMA;
        real RT=(gamma-1)*e;
        real p=r*RT;
        (*prim)(i,0)=r;
        (*prim)(i,1)=u;
        (*prim)(i,2)=v;
        (*prim)(i,3)=p;
    }
}


/**
 * @brief 三维欧拉方程的守恒变量到原始变量转换
 *
 * 将守恒变量（rho, rho*u, rho*v, rho*w, rho*E）转换为原始变量（rho, u, v, w, p）。
 * 结构与 consToPrimEuler2D 一致，仅增加 z 向动量/速度分量（三维欧拉 5 变量）。
 */
void Equation::consToPrimEuler3D()
{
    if(nCons!=5||nPrim!=5)
    {
        std::cout<<"Equation error: Euler 3d equation variable number error \n";
    }

    for (int i = 0; i < n; i++)
    {
        real r=(*cons)(i,0);    // 密度
        real ru=(*cons)(i,1);   // x 向动量
        real rv=(*cons)(i,2);   // y 向动量
        real rw=(*cons)(i,3);   // z 向动量（三维新增）
        real rE=(*cons)(i,4);   // 总能
        real u=ru/r;
        real v=rv/r;
        real w=rw/r;
        real E=rE/r;
        real q2=(u*u+v*v+w*w)/2;   // 动能（含 w²）
        real e=-q2+E;
        real gamma=GAMMA;
        real RT=(gamma-1)*e;
        real p=r*RT;
        (*prim)(i,0)=r;
        (*prim)(i,1)=u;
        (*prim)(i,2)=v;
        (*prim)(i,3)=w;   // z 向速度（三维新增）
        (*prim)(i,4)=p;
    }
}


/**
 * @brief 获取原始变量数据指针
 * @return 原始变量数据指针
 */
Data* Equation::getPrim()
{
    return prim;
}

/**
 * @brief 获取守恒变量数据指针
 * @return 守恒变量数据指针
 */
Data* Equation::getCons()
{
    return cons;
}

/**
 * @brief 获取右手端数据指针
 * @return 右手端数据指针
 */
Data* Equation::getRhs()
{
    return rhs;
}