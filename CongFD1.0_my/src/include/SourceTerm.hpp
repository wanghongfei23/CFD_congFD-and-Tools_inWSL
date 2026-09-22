#pragma once 
#include "data.hpp"
#include "info.hpp"

class SourceTerm
{
    public:

    SourceTerm(Data*,Data*,Info*);
    void calSource();
    

    private:
    friend class Initializer;
    void calGravitySource();
    void calViscousMMS();
    void nothingHappened();
    void (SourceTerm::*calSourceMethod)()=nullptr;
    Data *rhs;
    Data *prim;
    Info *info;
    int n,nprim,nCons;
    int type=0;
};