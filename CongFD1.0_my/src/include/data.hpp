#pragma once
#include "macro.hpp"

class Data
{
    public:

    Data(){};
    Data(int,int);
    Data(const Data&);

    void solInit(int,int);
    void init(int,int);
    void setValue(std::vector<real>);
    void setvarName(std::vector<std::string>&&);
    //void setDim(int,std::vector<int>);
    real& operator() (int,int);
    real& operator[] (int);
    
    //void operator= (Data&);
    void operator+= (std::vector<real>);
    void setZeros();
    real size();
    int getNVar();
    int getN();


    std::vector<real> getIVar(int);
    std::vector<real>::iterator begin();
    std::vector<real>::iterator end();
    
    //for global LF flux in burgers equation
    real maxElement(int);

    real getL2(int ivar);
    real getLinf(int ivar);
    std::vector<std::string> varName;

    private:
    std::vector<real> data;
    int n=200;
    int nVar=1;
};