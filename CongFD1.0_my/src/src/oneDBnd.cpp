/**
 * @file oneDBnd.cpp
 * @brief 一维边界条件类的实现文件
 */

#include "oneDBnd.hpp"

/**
 * @brief 构造函数
 * @param n_ 网格点数
 * @param nVar_ 变量数
 * @param bType_ 边界类型
 */
OneDBnd::OneDBnd(int n_,int nVar_,BndType bType_)
{
    n=n_;
    nVar=nVar_;
    data.resize(n*nVar,0.0);
    type=bType_;
}

/**
 * @brief 重载括号运算符，用于访问边界数据
 * @param i 索引
 * @param ivar 变量索引
 * @return 数据引用
 */
real& OneDBnd::operator()(int i,int ivar)
{
    return data.at(i*nVar+ivar);
}

/**
 * @brief 获取网格点数
 * @return 网格点数
 */
int OneDBnd::getN()
{
    return n;
}


/**
 * @brief 设置边界值
 * @param value 值向量
 */
void OneDBnd::setValue(std::vector<real> value)
{
    if(value.size()!=data.size())
    {
        std::cout<<"OneDBnd error: setValue() incorrect size\n";
    }
    std::copy(value.begin(),value.end(),data.begin());
}

/**
 * @brief 获取边界类型
 * @return 边界类型
 */
BndType OneDBnd::getType()
{
    return type;
}

/**
 * @brief 设置更新参数
 * @param prim_ 原始变量数据指针
 * @param i0_ 起始索引
 * @param offset_ 偏移量
 */
void OneDBnd::setUpdate(Data* prim_,int i0_,int offset_)
{
    prim=prim_;
    i0=i0_;
    offset=offset_;
}

/**
 * @brief 更新边界数据
 * 
 * 根据边界类型更新边界数据，支持周期边界、超声速出口边界、
 * 对称边界和双马赫反射问题的上边界等。
 */
void OneDBnd::update()
{
    switch (type)
    {
    case PERIODIC1D:
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < nVar; j++)
            {
                data[i*nVar+j]=(*prim)(i0+i*offset,j);
            }
            
        }
        break;
    case SUPERSONICOUTLET:
        for (int i = 0; i < n; i++) 
        {
            for (int j = 0; j < nVar; j++)
            {
                data[i*nVar+j]=(*prim)(i0,j);
            }
        }
        break;
    case SYMMETRYX:
        //only for 2D
        for (int i = 0; i < n; i++) 
        {
            // data[i*nVar+0]=(*prim)(i0,0);
            // data[i*nVar+1]=-(*prim)(i0,1);
            // data[i*nVar+2]=(*prim)(i0,2);
            // data[i*nVar+3]=(*prim)(i0,3);

            data[i*nVar+0]=(*prim)(i0+i*offset,0);
            data[i*nVar+1]=-(*prim)(i0+i*offset,1);
            data[i*nVar+2]=(*prim)(i0+i*offset,2);
            data[i*nVar+3]=(*prim)(i0+i*offset,3);
        }
        break;
    case SYMMETRY1D:
        for (int i = 0; i < n; i++) 
        {
            // data[i*nVar+0]=(*prim)(i0,0);
            // data[i*nVar+1]=-(*prim)(i0,1);
            // data[i*nVar+2]=(*prim)(i0,2);
            // data[i*nVar+3]=(*prim)(i0,3);

            data[i*nVar+0]=(*prim)(i0+i*offset,0);
            data[i*nVar+1]=-(*prim)(i0+i*offset,1);
            data[i*nVar+2]=(*prim)(i0+i*offset,2);
        }
        break;
    case SYMMETRYY:
        //only for 2D
        for (int i = 0; i < n; i++) 
        {
            data[i*nVar+0]=(*prim)(i0+i*offset,0);
            data[i*nVar+1]=(*prim)(i0+i*offset,1);
            data[i*nVar+2]=-(*prim)(i0+i*offset,2);
            data[i*nVar+3]=(*prim)(i0+i*offset,3);

            // data[i*nVar+0]=(*prim)(i0,0);
            // data[i*nVar+1]=(*prim)(i0,1);
            // data[i*nVar+2]=-(*prim)(i0,2);
            // data[i*nVar+3]=(*prim)(i0,3);
        }
        break;
    case DoubleMachUp:
        //only for 2D
        {
        for (int i = 0; i < n; i++) 
        {
            real x=coor[0];
            real y=-dh[1]/2+i*dh[1];
            //real y=dh[1]/2;
            real gt=1.0/6.0+std::sqrt(3.0)/3.0*(1.0+20*info->t);
            std::array<real,4> exactValues;
            if(y>std::sqrt(3.0)*(x-gt))
            {
                exactValues={8.0,8.25*cos(M_PI/6),-8.25*sin(M_PI/6),116.5};
            }
            else{
                exactValues={1.4,0,0,1.0};
            }
            
                data[i*nVar+0]=exactValues[0];
                data[i*nVar+1]=exactValues[1];
                data[i*nVar+2]=exactValues[2];
                data[i*nVar+3]=exactValues[3];
        }
        }
        break;
    case ShockShearIn:
        //only for 2D：激波/剪切层算例左入流（Peng 2019 §4.2.3 式(72)-(74)：基态 + 时变扰动 v'）
        {
        for (int i = 0; i < n; i++)
        {
            real y=coor[1];
            real rho=(y>=0.0)? 1.6374 : 0.3626;
            real u=2.5+0.5*std::tanh(2.0*y);
            real T=30.0/2.68;
            real vp=0.05*std::cos(2.0*M_PI*info->t/T)*std::exp(-y*y/10.0)
                   +0.05*std::cos(4.0*M_PI*info->t/T+M_PI/2.0)*std::exp(-y*y/10.0);
                data[i*nVar+0]=rho;
                data[i*nVar+1]=u;
                data[i*nVar+2]=vp;
                data[i*nVar+3]=0.3327;
        }
        }
        break;

    default:
        break;
    }
    
    
}

/**
 * @brief 设置Info对象指针
 * @param info_ Info对象指针
 */
void OneDBnd::setInfo(Info* info_){info=info_;}

/**
 * @brief 设置坐标信息
 * @param coor_ 坐标数组
 * @param dh_ 网格间距数组
 */
void OneDBnd::setCoor(std::array<real,3> coor_,std::array<real,3> dh_){coor=coor_;dh=dh_;}