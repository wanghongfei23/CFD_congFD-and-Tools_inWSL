/**
 * @file info.cpp
 * @brief Info类的实现文件
 */

#include "info.hpp"

/**
 * @brief 获取虚拟点单元数
 * @return 虚拟点单元数
 */
int Info::nGhostCell()
{
    
    if(WCNS5==spMethod && TRAD6 == diffMethod) return 5;
    else if(WCNS5==spMethod && HDS6 == diffMethod) return 3;
    else if(WCNS5==spMethod && MND6 == diffMethod) return 4;
    else if(MUSCL==spMethod && TRAD2 == diffMethod) return 2;
    else if(FIRSTORDER==spMethod && TRAD2 == diffMethod) return 1;

    else
    {
        std::cout<<"Info error: undifined spMethod and diffMethod combination\n";
    }
    return 0;
}

/**
 * @brief 获取通量点数
 * @return 通量点数
 */
int Info::nFluxPoint()
{
    switch (diffMethod)
    {
    case TRAD2:
    case HDS6:
        return 0;
        break;
    case MND6:
        return 1;
        break;
    case TRAD6:
        return 2;
        break;
    
    default:
        return 0;
        break;
    }
}

/**
 * @brief 获取原始变量数
 * @return 原始变量数
 */
int Info::nPrim()
{
    switch (eqType)
    {
    case ACCURACYTEST:
    case LINEARCONV1D:
    case BURGERS1D:
        return 1;
        break;
    case EULER:
        {
            if(dim==1) return 3;
            if(dim==2) return 4;
            if(dim==3) return 5;   // 三维欧拉：共 5 个变量（在二维基础上增加 z 向分量）
            else return 0;
        }
        break;
    
    default:
        std::cout<<"Info error: undifined eqType in nPrim()\n";
        return 0;
        
        break;
    }
}

/**
 * @brief 获取守恒变量数
 * @return 守恒变量数
 */
int Info::nCons()
{
    switch (eqType)
    {
    case ACCURACYTEST:
    case LINEARCONV1D:
    case BURGERS1D:
        return 1;
        break;
    case EULER:
        {
            if(dim==1) return 3;
            if(dim==2) return 4;
            if(dim==3) return 5;   // 三维欧拉：共 5 个变量（在二维基础上增加 z 向分量）
            else return 0;
        }
        break;
    
    default:
        std::cout<<"Info error: undifined eqType in nCons()\n";
        return 0;
        
        break;
    }
}


/**
 * @brief 获取问题维度
 * @return 问题维度
 */
int Info::getDim()
{
    int dim=0;
    for (int i = 0; i < 3; i++)
    {
        if (iMax[i]>2) dim++;
    }
    return dim;
}

/**
 * @brief 获取内部网格最大索引
 * @return 内部网格最大索引数组
 */
std::array<int,3> Info::icMax()
{
    std::array<int,3> icMax;
    for (int i = 0; i < 3; i++)
    {
        if (iMax[i]<2) icMax[i]=1;
        else
        icMax[i]=iMax[i]-1;
    }
    return icMax;
}


/**
 * @brief 获取默认边界类型
 * @return 默认边界类型
 */
BndType Info::defaultBndType()
{
    switch (eqType)
    {
    case LINEARCONV1D:
    case BURGERS1D:
        return PERIODIC1D;
        break;
    
    default:
        return TYPENULL;
        break;
    }
}

static std::map<EquationType,std::string> fluxStr= {{LINEARCONV1D,"LINEARCONV1D"}
                                        ,{BURGERS1D,"BURGERS1D"}
                                        ,{EULER,"EULER"}
                                        ,{ACCURACYTEST,"ACCURACYTEST"}};
static std::map<InterMethod,std::string> disStr={
    {FIRSTORDER,"FIRSTORDER"},
    {MUSCL,"MUSCL"},
    // 【王鸿飞】begin-1命名
    {WCNS5,"WENO-JS"},
    {WCNSZ5,"WENO-Z"},
    {TCNS5,"TENO"},
    {WCNS5CONGZ,"TENO-S"},
    {WHFTCNSA,"TENO-A"},
    {WHFTCNSASF002,"TENO-AS-myF002"},
    {WHFTCNSAH002,"TENO-AS-myH002"},
    {WHFTCNSASF102,"TENO-AS-myF102"},
    {WHFTCNSASF103,"TENO-AS-myF103"},
    {WHFTCNSASF102_reciprocal,"TENO-AS-myF102_reciprocal"},
    {WHFTCNSASF103_reciprocal,"TENO-AS-myF103_reciprocal"},
    {WHFTCNSAS_fx,"TENO-AS-fx"},
    {WHFTCNSAS_initial,"TENO-AS-initial"},
    {WHFTCNSAS_approx_1,"TENO-AS-approx_1"},
    {WHFTCNSAS_fx_real,"TENO-AS-fx-real"},
    {WHFTCNSAS_approx_2,"TENO-AS-approx_2"},
    {WHFTCNSASF202_2S,"TENO-AS-myF202_2S"},
    {WHFTCNSASF202,"TENO-AS-myF202"},
    {WHFTCNSASFf_5_10,"TENO-AS-Ff_5_10"},
    {WHFTCNSASFf_5_9,"TENO-AS-Ff_5_9"},
    {WHFTCNSASFf2_test,"TENO-AS-Ff2_test"},
    {WHFTCNSASFf3_test,"TENO-AS-Ff3_test"},
    {WHFTCNSASFf3_5_9_time_improve,"TENO-AS-Ff3_5_9_time_improve"},
    {temp015,"temp_name_015"},
    {WHFTCNSLADSFf2_5_10,"TENO-LADS-Ff2_5_10"},
    {WHFTCNSLAD,"TENO-LAD"},
    {WHFTCNSLADS_g1,"TENO-LADS-g1"},
    {temp019,"temp_name_019"},
    {temp020,"temp_name_020"},
    {temp021,"temp_name_021"},
    {temp022,"temp_name_022"},
    {temp023,"temp_name_023"},
    {temp024,"temp_name_024"},
    {temp025,"temp_name_025"},
    {temp026,"temp_name_026"},
    {temp027,"temp_name_027"},
    {temp028,"temp_name_028"},
    {temp029,"temp_name_029"}
    // 【王鸿飞】end-1命名
};

// 新命名系统

// 定义静态映射表，将算例映射到字符串
static std::map<int,std::string> exampleStr1D={
    {0,"Sod"},
    {1,"ShuOsher"},
    {2,"Lax"},
    {3,"sedov"},
    {4,"Woodward_Colella"},
    {5,"Double_sparse_wave"}
};

/**
 * @brief 二维算例标识到字符串的映射表
 * 
 * 参考文献:
 * - Schulz-Rinne, C. W., Collins, J. P., & Glaz, H. M. (2004). 
 *   "Numerical Solution of the Two-Dimensional Riemann Problems for Gas Dynamics"
 *   SIAM Journal on Scientific Computing, 25(6), 1895-1917.
 * 
 * 该文献系统研究了19种不同的二维黎曼问题配置(Configuration 1-19)。
 * 本代码实现了其中常用的几种配置。
 */
static std::map<int,std::string> exampleStr2D={
    {0,"2D_Riemann_Config3"},      // Configuration 3: 双激波+接触间断
    {1,"2D_Riemann_Config4"},      // Configuration 4: 涡旋结构(Vortex)
    {2,"implosion"},               // 内爆问题(Implosion Problem)
    {3,"RTI"},                     // 瑞利-泰勒不稳定性(Rayleigh-Taylor Instability)
    {4,"Double_Mach"},             // 双马赫反射(Double Mach Reflection)
    {5,"2D_Riemann_Config12"},     // Configuration 12: 混合波系
    {6,"KHI"},                     // 开尔文-亥姆霍兹不稳定性(Kelvin-Helmholtz Instability)
    {7,"2D_Riemann_Config6"}       // Configuration 6: 反向涡旋结构
};

/**
 * @brief 三维算例标识到字符串的映射表（任务 6 新增）
 */
static std::map<int,std::string> exampleStr3D={
    {0,"3D_TGV"},       // 三维 Taylor-Green 涡
    {1,"3D_uniform"}    // 均匀场（三维核验用）
};

/**
 * @brief 生成文件名
 * @return 生成的文件名
 */
std::string Info::filename()
{
    // 字符串拼接，包含空间离散方法、插值方法、算例信息和当前时间
    std::string caseName = "unknown";
    
    if (dim == 1) {
        auto it = exampleStr1D.find(nCase);
        if (it != exampleStr1D.end()) {
            caseName = it->second;
        }
    } else if (dim == 2) {
        auto it = exampleStr2D.find(nCase);
        if (it != exampleStr2D.end()) {
            caseName = it->second;
        }
    } else if (dim == 3) {
        auto it = exampleStr3D.find(nCase);
        if (it != exampleStr3D.end()) {
            caseName = it->second;
        }
    }
    
    // 添加网格信息到文件名
    std::string gridInfo = "";
    if (dim == 1) {
        gridInfo = std::to_string(iMax[0]);
    } else if (dim == 2) {
        gridInfo = std::to_string(iMax[0]) + "x" + std::to_string(iMax[1]);
    } else if (dim == 3) {
        gridInfo = std::to_string(iMax[0]) + "x" + std::to_string(iMax[1]) + "x" + std::to_string(iMax[2]);
    }
    
    return caseName + " - " + disStr[interMethod] + " - " + gridInfo + " - " + std::format("t={:.4f}.cgns",t);
}


// 原命名系统
// std::string Info::filename()
// {
//     return fluxStr[eqType]+disStr[spMethod]+std::format("t={:.4f}.cgns",t);
// }


std::vector<std::string> Info::getVarNameListCons()
{
    std::vector<std::string> res;
    res.reserve(nCons());
    if(dim==1)
    {
        if(eqType==EULER)
        {
            res.push_back("rho");
            res.push_back("rhoU");
            res.push_back("rhoE");
        }
        else res.push_back("u");
    }
    else if (dim==2)
    {
        if(eqType==EULER)
        {
            res.push_back("rho");
            res.push_back("rhoU");
            res.push_back("rhoV");
            res.push_back("rhoE");
        }
    }
    else if (dim==3)
    {
        if(eqType==EULER)
        {
            // 三维守恒变量：rho, rhoU, rhoV, rhoW, rhoE
            res.push_back("rho");
            res.push_back("rhoU");
            res.push_back("rhoV");
            res.push_back("rhoW");
            res.push_back("rhoE");
        }
    }
    return res;
}
std::vector<std::string> Info::getVarNameListPrim()
{
    std::vector<std::string> res;
    res.reserve(nPrim());
    if(dim==1)
    {
        if(eqType==EULER)
        {
            res.push_back("Density");
            res.push_back("XVelocity");
            res.push_back("Pressure");
            // res.push_back("Density");
            // res.push_back("u+c");
            // res.push_back("u-c");
        }
        else res.push_back("u");
    }
    else if (dim==2)
    {
        if(eqType==EULER)
        {
            res.push_back("Density");
            res.push_back("XVelocity");
            res.push_back("YVelocity");
            res.push_back("Pressure");
        }
    }
    else if (dim==3)
    {
        if(eqType==EULER)
        {
            // 三维原始变量：密度、三个方向速度、压力
            res.push_back("Density");
            res.push_back("XVelocity");
            res.push_back("YVelocity");
            res.push_back("ZVelocity");
            res.push_back("Pressure");
        }
    }
    return res;
}

std::vector<std::string> Info::getVarNameListRhs()
{
    std::vector<std::string> res;
    res.reserve(nCons());
    if(dim==1)
    {
        if(eqType==EULER)
        {
            res.push_back("RHS-rho");
            res.push_back("RHS-rhoU");
            res.push_back("RHS-rhoE");
        }
        else res.push_back("RHS-u");
    }
    else if (dim==2)
    {
        if(eqType==EULER)
        {
            res.push_back("RHS-rho");
            res.push_back("RHS-rhoU");
            res.push_back("RHS-rhoV");
            res.push_back("RHS-rhoE");
        }
    }
    else if (dim==3)
    {
        if(eqType==EULER)
        {
            // 三维右端项变量名（与守恒变量一一对应）
            res.push_back("RHS-rho");
            res.push_back("RHS-rhoU");
            res.push_back("RHS-rhoV");
            res.push_back("RHS-rhoW");
            res.push_back("RHS-rhoE");
        }
    }
    return res;
}


real Info::geth(int idim)
{
    if(constH) return interval;
    return (calZone[2*idim+1]-calZone[2*idim])/(iMax[idim]-1);
}

Info::Info()
{
    dim=getDim();
}