/**
 * @file SpaceDis.cpp
 * @brief 空间离散类的实现文件
 */

#include "SpaceDis.hpp"
#include "interScheme.hpp"

/**
 * @brief 默认构造函数
 */
SpaceDis::SpaceDis() { };

/**
 * @brief 构造函数
 * @param n_ 网格点数
 * @param data_ 数据指针
 * @param rhs_ 右手端数据指针
 * @param bndL_ 左边界指针
 * @param bndR_ 右边界指针
 * @param info_ Info对象指针
 */
SpaceDis::SpaceDis(int n_, Data* data_, Data* rhs_, std::shared_ptr<OneDBnd> bndL_, std::shared_ptr<OneDBnd> bndR_, Info* info_)
{
    data = data_;
    info = info_;
    nPrim = info->nPrim();
    n = n_;
    nVar = info->nCons();
    nHalf = n + 1;
    bndL = bndL_;
    bndR = bndR_;
    flux_d = std::make_shared<Data>(nHalf, nVar);

    fBndL = std::make_shared<OneDBnd>(info->nFluxPoint(), nVar, FLUXGHOST);
    fBndR = std::make_shared<OneDBnd>(info->nFluxPoint(), nVar, FLUXGHOST);
    rhs = rhs_;

    fluxType = info->eqType;
    diffMethod = info->diffMethod;
    interMethod = info->interMethod;
    switch (fluxType) {
    case ACCURACYTEST:
        calTypeFlux = &SpaceDis::calFluxAccuracyTest;
        break;
    case LINEARCONV1D:
        calTypeFlux = &SpaceDis::calFluxConv;
        break;
    case BURGERS1D:
        calTypeFlux = &SpaceDis::calFluxBurgers;
        break;
    case EULER:
        if (info->dim == 1) {
            calTypeFlux = &SpaceDis::calFluxEuler1D;
        } else if (info->dim == 2) {
            calTypeFlux = &SpaceDis::calFluxEuler2D;
        } else if (info->dim == 3) {   // 三维（任务 5 新增）
            calTypeFlux = &SpaceDis::calFluxEuler3D;
        }

        break;

    default:
        break;
    }

    switch (diffMethod) {
    case TRAD2:
        difMethod = &SpaceDis::dif2Order;
        break;
    case TRAD6:
        difMethod = &SpaceDis::difTraditional6;
        break;
    case HDS6:
        difMethod = &SpaceDis::difHCS;
        break;
    case MND6:
        difMethod = &SpaceDis::difMND6;
        break;

    default:
        break;
    }
    // FIRSTORDER,
    //  MUSCL,
    //  WCNS5,
    //  WCNSZ5,
    //  WCNS5Char,
    //  WCNSZ5Char,
    //  WCNS5CONG,
    //  TCNS5
    switch (interMethod) {
    case FIRSTORDER:
    case MUSCL:
        reconLMethod = &SpaceDis::reconL;
        reconRMethod = &SpaceDis::reconR;
        break;
    case WCNS5:
        inter5 = &weno5_JSchen;
        inter5Positive = &weno5_JSchen;
        reconLMethod = &SpaceDis::reconLprim;
        reconRMethod = &SpaceDis::reconRprim;
        break; 
    case WCNSZ5:
        inter5 = &weno5_Z;
        inter5Positive = &weno5_Z;
        reconLMethod = &SpaceDis::reconLprim;
        reconRMethod = &SpaceDis::reconRprim;
        break;
    case WCNS5Char:
        inter5 = &weno5_JSchen;
        inter5Positive = &weno5_JSchen;
        reconLMethod = (info->dim == 1) ? (&SpaceDis::reconLChar1D) : (&SpaceDis::reconLChar2D);
        reconRMethod = (info->dim == 1) ? (&SpaceDis::reconRChar1D) : (&SpaceDis::reconRChar2D);
        break;
    case WCNSZ5Char:
        inter5 = &weno5_Z;
        inter5Positive = &weno5_Z;
        reconLMethod = (info->dim == 1) ? (&SpaceDis::reconLChar1D) : (&SpaceDis::reconLChar2D);
        reconRMethod = (info->dim == 1) ? (&SpaceDis::reconRChar1D) : (&SpaceDis::reconRChar2D);
        break;
    case TCNS5:
        inter5 = &Teno5_Z;
        inter5Positive = &Teno5_Z;
        if (fluxType == EULER) {
            reconLMethod = (info->dim == 1) ? (&SpaceDis::reconLChar1D) : (&SpaceDis::reconLChar2D);
            reconRMethod = (info->dim == 1) ? (&SpaceDis::reconRChar1D) : (&SpaceDis::reconRChar2D);
        } else {
            reconLMethod = &SpaceDis::reconLprim;
            reconRMethod = &SpaceDis::reconRprim;
        }
        break;
    case TCNS5CT4:
        inter5 = &Teno5_ZCT4;
        inter5Positive = &Teno5_ZCT4;
        if (fluxType == EULER) {
            reconLMethod = (info->dim == 1) ? (&SpaceDis::reconLChar1D) : (&SpaceDis::reconLChar2D);
            reconRMethod = (info->dim == 1) ? (&SpaceDis::reconRChar1D) : (&SpaceDis::reconRChar2D);
        } else {
            reconLMethod = &SpaceDis::reconLprim;
            reconRMethod = &SpaceDis::reconRprim;
        }
        break;
    case TCNS5CT7:
        inter5 = &Teno5_ZCT7;
        inter5Positive = &Teno5_ZCT7;
        if (fluxType == EULER) {
            reconLMethod = (info->dim == 1) ? (&SpaceDis::reconLChar1D) : (&SpaceDis::reconLChar2D);
            reconRMethod = (info->dim == 1) ? (&SpaceDis::reconRChar1D) : (&SpaceDis::reconRChar2D);
        } else {
            reconLMethod = &SpaceDis::reconLprim;
            reconRMethod = &SpaceDis::reconRprim;
        }
        break;
    case WCNS5CONG:
        inter5 = &Teno5_Cong;
        inter5Positive = &Teno5_Cong;
        if (fluxType == EULER) {
            reconLMethod = (info->dim == 1) ? (&SpaceDis::reconLChar1D) : (&SpaceDis::reconLChar2D);
            reconRMethod = (info->dim == 1) ? (&SpaceDis::reconRChar1D) : (&SpaceDis::reconRChar2D);
        } else {
            reconLMethod = &SpaceDis::reconLprim;
            reconRMethod = &SpaceDis::reconRprim;
        }
        break;
    case TCNSCongA:
        inter5 = &Teno5_CongA;
        inter5Positive = &Teno5_CongA;
        if (fluxType == EULER) {
            reconLMethod = (info->dim == 1) ? (&SpaceDis::reconLChar1D) : (&SpaceDis::reconLChar2D);
            reconRMethod = (info->dim == 1) ? (&SpaceDis::reconRChar1D) : (&SpaceDis::reconRChar2D);
        } else {
            reconLMethod = &SpaceDis::reconLprim;
            reconRMethod = &SpaceDis::reconRprim;
        }
        break;
    case WCNS5CONGZ:
        inter5 = &Teno5_CongZ;
        inter5Positive = &Teno5_CongZ;
        if (fluxType == EULER) {
            reconLMethod = (info->dim == 1) ? (&SpaceDis::reconLChar1D) : (&SpaceDis::reconLChar2D);
            reconRMethod = (info->dim == 1) ? (&SpaceDis::reconRChar1D) : (&SpaceDis::reconRChar2D);
        } else {
            reconLMethod = &SpaceDis::reconLprim;
            reconRMethod = &SpaceDis::reconRprim;
        }
        break;
        // WCNS5CONGZCT4
    case WCNS5CONGZCT4:
        inter5 = &Teno5_CongZCT4;
        inter5Positive = &Teno5_CongZCT4; 
        if (fluxType == EULER) {
            reconLMethod = (info->dim == 1) ? (&SpaceDis::reconLChar1D) : (&SpaceDis::reconLChar2D);
            reconRMethod = (info->dim == 1) ? (&SpaceDis::reconRChar1D) : (&SpaceDis::reconRChar2D);
        } else {
            reconLMethod = &SpaceDis::reconLprim;
            reconRMethod = &SpaceDis::reconRprim;
        }
        break;
    case WCNS5CONGZCT7:
        inter5 = &Teno5_CongZCT7;
        inter5Positive = &Teno5_CongZCT7;
        if (fluxType == EULER) {
            reconLMethod = (info->dim == 1) ? (&SpaceDis::reconLChar1D) : (&SpaceDis::reconLChar2D);
            reconRMethod = (info->dim == 1) ? (&SpaceDis::reconRChar1D) : (&SpaceDis::reconRChar2D);
        } else {
            reconLMethod = &SpaceDis::reconLprim;
            reconRMethod = &SpaceDis::reconRprim;
        }
        break;
    case LINEAR5:
        inter5 = &Linear5;
        inter5Positive = &Linear5;
        if (fluxType == EULER) { 
            reconLMethod = (info->dim == 1) ? (&SpaceDis::reconLChar1D) : (&SpaceDis::reconLChar2D);
            reconRMethod = (info->dim == 1) ? (&SpaceDis::reconRChar1D) : (&SpaceDis::reconRChar2D);
        } else {
            reconLMethod = &SpaceDis::reconLprim;
            reconRMethod = &SpaceDis::reconRprim;
        }
        break;
    case NICEST5:
        inter5 = &Nicest5;
        inter5Positive = &Nicest5;
        if (fluxType == EULER) {
            reconLMethod = (info->dim == 1) ? (&SpaceDis::reconLChar1D) : (&SpaceDis::reconLChar2D);
            reconRMethod = (info->dim == 1) ? (&SpaceDis::reconRChar1D) : (&SpaceDis::reconRChar2D);
        } else {
            reconLMethod = &SpaceDis::reconLprim;
            reconRMethod = &SpaceDis::reconRprim;
        }
        break;

    //【王鸿飞】begin-1

    // 王鸿飞统计CT_A
    case WHFTCNSA:
        inter5 = &whf_TCNS_A;
        inter5Positive = &whf_TCNS_A;
        // inter5 = &whf_TCNS_A_count;
        // inter5Positive = &whf_TCNS_A_count;
        if (fluxType == EULER) {
            reconLMethod = (info->dim == 1) ? (&SpaceDis::reconLChar1D) : (&SpaceDis::reconLChar2D);
            reconRMethod = (info->dim == 1) ? (&SpaceDis::reconRChar1D) : (&SpaceDis::reconRChar2D);
        } else {
            reconLMethod = &SpaceDis::reconLprim;
            reconRMethod = &SpaceDis::reconRprim;
        }
        break;
    













    case WHFTCNSASFf3_5_9_time_improve:
        inter5 = &whf_TCNS_AS_Ff3_5_9_time_improve;
        inter5Positive = &whf_TCNS_AS_Ff3_5_9_time_improve;
        if (fluxType == EULER) {
            reconLMethod = (info->dim == 1) ? (&SpaceDis::reconLChar1D) : (&SpaceDis::reconLChar2D);
            reconRMethod = (info->dim == 1) ? (&SpaceDis::reconRChar1D) : (&SpaceDis::reconRChar2D);
        } else {
            reconLMethod = &SpaceDis::reconLprim;
            reconRMethod = &SpaceDis::reconRprim;
        }
        break;


    case WHFTCNSLADSFf2_5_10:
        inter5 = &whf_TCNS_LADS_Ff2_5_10;
        inter5Positive = &whf_TCNS_LADS_Ff2_5_10;
        if (fluxType == EULER) {
            reconLMethod = (info->dim == 1) ? (&SpaceDis::reconLChar1D) : (&SpaceDis::reconLChar2D);
            reconRMethod = (info->dim == 1) ? (&SpaceDis::reconRChar1D) : (&SpaceDis::reconRChar2D);
        } else {
            reconLMethod = &SpaceDis::reconLprim;
            reconRMethod = &SpaceDis::reconRprim;
        }
        break;

    case WHFTCNSLAD:
        inter5 = &whf_TCNS_LAD;
        inter5Positive = &whf_TCNS_LAD;
        if (fluxType == EULER) {
            reconLMethod = (info->dim == 1) ? (&SpaceDis::reconLChar1D) : (&SpaceDis::reconLChar2D);
            reconRMethod = (info->dim == 1) ? (&SpaceDis::reconRChar1D) : (&SpaceDis::reconRChar2D);
        } else {
            reconLMethod = &SpaceDis::reconLprim;
            reconRMethod = &SpaceDis::reconRprim;
        }
        break;

    case WHFTCNSLADS_g1:
        inter5 = &whf_TCNS_LADS_g1;
        inter5Positive = &whf_TCNS_LADS_g1;
        if (fluxType == EULER) {
            reconLMethod = (info->dim == 1) ? (&SpaceDis::reconLChar1D) : (&SpaceDis::reconLChar2D);
            reconRMethod = (info->dim == 1) ? (&SpaceDis::reconRChar1D) : (&SpaceDis::reconRChar2D);
        } else {
            reconLMethod = &SpaceDis::reconLprim;
            reconRMethod = &SpaceDis::reconRprim;
        }
        break;














    //【王鸿飞】end-1

    default:
        break;
    }
}

void SpaceDis::setOffset(int i0_, int offset_)
{
    i0 = i0_;
    offset = offset_;
}

void SpaceDis::difference()
{
    calFlux();             //计算通量
    (this->*difMethod)();  // 差分
}
void SpaceDis::calFlux()
{
    // 获取左边界和右边界ghost cell的数量，用于确定计算范围
    int fGhostL = fBndL->getN(), fGhostR = fBndR->getN();
    // 将通量数组初始化为零
    flux_d->setZeros();
    // 在包括ghost cells在内的所有单元格面上计算通量
    // 计算范围从最左边的ghost cell面(索引为0-fGhostL)到最右边的ghost cell面(索引为nHalf+fGhostR-1)
    for (int i = 0 - fGhostL; i < nHalf + fGhostR; i++) {
        // 使用函数指针调用相应方程类型的通量计算方法
        // 根据问题类型(calTypeFlux)调用对应的通量计算函数
        (this->*calTypeFlux)(i);
    }
}

void SpaceDis::setMethod(EquationType type_, DiffMethod method_)
{
    fluxType = type_;
    diffMethod = method_;
    switch (fluxType) {
    case LINEARCONV1D:
        calTypeFlux = &SpaceDis::calFluxConv;
        break;
    case BURGERS1D:
        calTypeFlux = &SpaceDis::calFluxBurgers;
        break;
    case EULER:
        calTypeFlux = &SpaceDis::calFluxEuler1D;
        break;

    default:
        break;
    }

    switch (diffMethod) {
    case TRAD2:
        difMethod = &SpaceDis::dif2Order;
        break;
    case TRAD6:
        difMethod = &SpaceDis::difTraditional6;
        break;
    case HDS6:
        difMethod = &SpaceDis::difHCS;
        break;
    case MND6:
        difMethod = &SpaceDis::difMND6;
        break;

    default:
        break;
    }
}

real& SpaceDis::at(int i, int ivar)
{
    // get the values in data
    if (i < 0) {
        return (*bndL)(-(i + 1), ivar);
    }
    if (i >= n) {
        return (*bndR)(i - n, ivar);
    }
    return (*data)[(i0 + offset * i) * nPrim + ivar];
}

real& SpaceDis::fluxAt(int i, int ivar)
{
    if (i < 0) {
        return (*fBndL)(-(i + 1), ivar);
    }
    if (i >= nHalf) {
        return (*fBndR)(i - nHalf, ivar);
    }
    return (*flux_d)[i * nVar + ivar];
}

void SpaceDis::setConstNorm(std::array<real, 3>&& norm_)
{
    norm = norm_;
}

void SpaceDis::setIDim(int idim_)
{
    // x:0,y:1,z:2 要直接用来作索引
    idim = idim_;
}
