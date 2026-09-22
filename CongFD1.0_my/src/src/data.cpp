/**
 * @file data.cpp
 * @brief Data类的实现文件
 */

#include "data.hpp"

/**
 * @brief 带参数的构造函数
 * @param n_ 数据点数
 * @param nVar_ 变量数
 */
Data::Data(int n_,int nVar_)
{
    n=n_;
    nVar=nVar_;
    data.resize(nVar*n,0.0);
}

/**
 * @brief 初始化解数据
 * @param n_ 数据点数
 * @param nvar_ 变量数
 */
void Data::solInit(int n_,int nvar_)
{
    n=n_;
    nVar=nvar_;
    data.resize(nVar*n,0.0);

    for(int i=0;i<n;i++)
    {
        real h=2.0/n;
        real xi=h/2.0+i*h-1.0;
        //data[i]=-sin(M_PI*xi);//for burgers equation

        //for sod tube 1D
        
        real gamma=GAMMA;
        if(xi<0)
        {
            (*this)(i,0)=1.0;
            (*this)(i,1)=0;
            (*this)(i,2)=1.0/(gamma-1)*1;
        }
        else
        {
            (*this)(i,0)=0.125;
            (*this)(i,1)=0;
            (*this)(i,2)=1.0/(gamma-1)*0.1;
        }
    }
}

/**
 * @brief 初始化数据
 * @param n_ 数据点数
 * @param nvar_ 变量数
 */
void Data::init(int n_,int nvar_)
{
    n=n_;
    nVar=nvar_;
    data.resize(nVar*n,0.0);
}

/**
 * @brief 重载括号运算符，用于访问特定位置的特定变量
 * @param i 数据点索引
 * @param ivar 变量索引
 * @return 对应位置的数据引用
 */
real& Data::operator() (int i,int ivar)
{
    return data[i*nVar+ivar];
}

/**
 * @brief 重载方括号运算符，用于访问特定位置的数据
 * @param i 数据索引
 * @return 对应位置的数据引用
 */
real& Data::operator[] (int i)
{
    //return (this->data)[i];
    return data.at(i);
}



/**
 * @brief 获取指定变量的最大元素值
 * @param ivar 变量索引
 * @return 最大元素值
 */
real Data::maxElement(int ivar)
{
    if(data.empty()) return 0;
    real res=data[0];
    for(int i=1;i<n;i++)
    {
        res=(res<data[i]) ? data[i] : res;
    }
    return res;
}



/**
 * @brief 设置数据值
 * @param value 数据值向量
 */
void Data::setValue(std::vector<real> value)
{
    assert(value.size()==data.size());
    std::copy(value.begin(),value.end(),data.begin());
}

// void Data::operator= (Data& dat)
// {
//     if(this->n==dat.n&&this->nVar==dat.nVar)
//     {
//         for(ind i = 0; i < n*nVar; i++)
//         {
//             (*this)[i]=dat[i];
//         }
//         
//     }
//     else
//     {
//         std::cout<<"incorrect size at class Data=Data\n";
//     }
// }

/**
 * @brief 重载加赋值运算符
 * @param arr 要加上的数据向量
 */
void Data::operator+= (std::vector<real> arr)
{
    assert(arr.size()==n*nVar);
    for (int i = 0; i < n*nVar; i++)
    {
        (*this)[i]+=arr[i];
    }
        
}

/**
 * @brief 将数据设置为零
 */
void Data::setZeros()
{
    std::fill(data.begin(),data.end(),0.0);
}



std::vector<real>::iterator Data::begin()
{
    return data.begin();
}
std::vector<real>::iterator Data::end()
{
    return data.end();
}

real Data::size()
{
    return data.size();
}


int Data::getNVar()
{
    return nVar;
}
int Data::getN()
{
    return n;
}

std::vector<real> Data::getIVar(int ivar)
{
    std::vector<real> res;
    assert(ivar<=nVar);
    
    
    res.reserve(n);
    for (int i = 0; i < n; i++)
    {
        res.push_back((*this)(i,ivar));
    }
    return res;
    
}

Data::Data(const Data& origin_)
{
    n=origin_.n;
    nVar=origin_.nVar;
    data.resize(n*nVar);
    std::copy(origin_.data.begin(),origin_.data.end(),data.begin());

}

void Data::setvarName(std::vector<std::string>&& varName_)
{
    assert(nVar==varName_.size());
    varName.resize(nVar);
    std::copy(varName_.begin(),varName_.end(),varName.begin());
}

real Data::getL2(int ivar)
{
    if(ivar>= nVar)
    {
        std::cout<<"Data error: ivar >= nVar";
        return 0;
    }
    real res=0;
    for (int i = 0; i < n; i++)
    {
        res+=pow((*this)(i,ivar),2);
    }
    res=std::sqrt(res/n);
    return res;
    
}
real Data::getLinf(int ivar)
{
    if(ivar>= nVar)
    {
        std::cout<<"Data error: ivar >= nVar";
        return 0;
    }
    real res=std::abs((*this)(0,ivar));
    for (int i = 1; i < n; i++)
    {
        if(res<std::abs((*this)(i,ivar))) res=std::abs((*this)(i,ivar));
    }
    return res;
    
}