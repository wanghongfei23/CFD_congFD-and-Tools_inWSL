/**
 * @file sp_recon.cpp
 * @brief 空间离散重构计算实现文件
 */

#include "SpaceDis.hpp"
#include "interScheme.hpp"

/**
 * @brief 右侧界面重构函数
 * @param i 网格点索引
 * @return 重构后的变量值向量
 */
std::vector<real> SpaceDis::reconR(int i)
{
    std::vector<real> res(nPrim);
    if (interMethod == MUSCL) {
        for (int ivar = 0; ivar < nPrim; ivar++) {
            res[ivar] = musclInterpolation(at(i + 1, ivar), at(i, ivar), at(i - 1, ivar));
        }
        return res;
    } else {
        for (int ivar = 0; ivar < nPrim; ivar++) {
            res[ivar] = at(i, ivar);
        }
        return res;
    }
}

/**
 * @brief 左侧界面重构函数
 * @param i 网格点索引
 * @return 重构后的变量值向量
 */
std::vector<real> SpaceDis::reconL(int i)
{
    std::vector<real> res(nPrim);
    if (interMethod == MUSCL) {
        for (int ivar = 0; ivar < nPrim; ivar++) {
            res[ivar] = musclInterpolation(at(i - 2, ivar), at(i - 1, ivar), at(i, ivar));
        }
        return res;
    } else {
        for (int ivar = 0; ivar < nPrim; ivar++) {
            res[ivar] = at(i - 1, ivar);
        }
        return res;
    }
}

/**
 * @brief 一维特征重构左侧函数
 * @param i 网格点索引
 * @return 重构后的变量值向量
 */
std::vector<real> SpaceDis::reconLChar1D(int i)
{
    assert(info->dim == 1);
    assert(info->eqType == EULER);
    enum { R,
        U,
        P };
    real rRef = at(i - 1, R);
    real pRef = at(i - 1, P);
    real cRef = std::sqrt(GAMMA * (pRef / rRef));
    std::array<real, 5> q1, q2, q3;
    for (int j = i - 3; j < i + 2; j++) {
        int iLocal = j - i + 3;
        q1[iLocal] = at(j, R) - at(j, P) / (cRef * cRef);
        // q2[iLocal]=(at(j,P)+at(j,U)*cRef*rRef)/2;
        // q3[iLocal]=(at(j,P)-at(j,U)*cRef*rRef)/2;
        q2[iLocal] = at(j, P);
        q3[iLocal] = at(j, U);
    }
    auto start = std::chrono::high_resolution_clock::now();
    real Q1 = (*this->inter5Positive)(q1);
    real Q2 = (*this->inter5)(q2);
    real Q3 = (*this->inter5)(q3);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
    timep += duration;
    // return {Q1+(Q2+Q3)/(cRef*cRef),(Q2-Q3)/cRef/rRef,Q2+Q3};
    return { Q1 + (Q2) / (cRef * cRef), Q3, Q2 };
}
/**
 * @brief 一维特征重构右侧函数
 * @param i 网格点索引
 * @return 重构后的变量值向量
 */
std::vector<real> SpaceDis::reconRChar1D(int i)
{
    enum { R,
        U,
        P };
    real rRef = at(i - 1, R);
    real pRef = at(i, P);
    real cRef = std::sqrt(GAMMA * (pRef / rRef));
    std::array<real, 5> q1, q2, q3;
    for (int j = i + 2; j > i - 3; j--) {
        int iLocal = i + 2 - j;
        q1[iLocal] = at(j, R) - at(j, P) / (cRef * cRef);
        // q2[iLocal]=(at(j,P)+at(j,U)*cRef*rRef)/2;
        // q3[iLocal]=(at(j,P)-at(j,U)*cRef*rRef)/2;
        q2[iLocal] = at(j, P);
        q3[iLocal] = at(j, U);
    }
    auto start = std::chrono::high_resolution_clock::now();
    real Q1 = (*this->inter5Positive)(q1);
    real Q2 = (*this->inter5)(q2);
    real Q3 = (*this->inter5)(q3);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
    timep += duration;
    // return {Q1+(Q2+Q3)/(cRef*cRef),(Q2-Q3)/cRef/rRef,Q2+Q3};
    return { Q1 + (Q2) / (cRef * cRef), Q3, Q2 };
}

/**
 * @brief 最小差值计算函数
 * @param u1L 左侧第一个值
 * @param u2L 左侧第二个值
 * @param u1R 右侧第一个值
 * @param u2R 右侧第二个值
 * @return 包含最小差值的数组
 */
static std::array<real, 2> minDif(real u1L, real u2L, real u1R, real u2R)
{
    std::array<real, 3> difs = { std::abs(u1L - u1R), std::abs(u1L - u2R), std::abs(u2L - u1R) };
    int index = std::min_element(difs.begin(), difs.end()) - difs.begin();
    if (index == 0)
        return { u1L, u1R };
    if (index == 1)
        return { u1L, u2R };
    if (index == 2)
        return { u2L, u1R };
    std::cout << "spRecon Error: minDif index out\n";
    return { 0, 0 };
}

/**
 * @brief 最小差值计算函数（三维数组版本）
 * @param uL 左侧值数组
 * @param uR 右侧值数组
 * @return 包含最小差值的数组
 */
static std::array<real, 2> minDif(std::array<real, 3> uL, std::array<real, 3> uR)
{
    std::array<real, 2> difs = { std::abs(uL[0] - uR[0]), std::abs(uL[1] - uR[1]) };
    int index = std::min_element(difs.begin(), difs.end()) - difs.begin();
    if (index == 0)
        return { uL[0], uR[0] };
    else if (index == 1)
        return { uL[1], uR[1] };
    std::cout << "spRecon Error: minDif index out\n";
    return { 0, 0 };
}

/**
 * @brief 最小差值计算函数（二维数组版本）
 * @param uL 左侧值数组
 * @param uR 右侧值数组
 * @return 包含最小差值的数组
 */
static std::array<real, 2> minDif(std::array<real, 2> uL, std::array<real, 2> uR)
{
    // std::array<real,4> difs={std::abs(uL[0]-uR[0]),std::abs(uL[1]-uR[1]),std::abs(uL[0]-uR[1]),std::abs(uL[1]-uR[0])};
    std::array<real, 2> difs = { std::abs(uL[0] - uR[0]), std::abs(uL[1] - uR[1]) };
    int index = std::min_element(difs.begin(), difs.end()) - difs.begin();
    if (index == 0)
        return { uL[0], uR[0] };
    else if (index == 1)
        return { uL[1], uR[1] };
    else if (index == 2)
        return { uL[0], uR[1] };
    else if (index == 3)
        return { uL[1], uR[0] };
    std::cout << "spRecon Error: minDif index out\n";
    return { 0, 0 };
}
// static std::array<real,2> minDif(std::array<real,3> uL,std::array<real,3> uR)
// {
//     std::array<real,4> difs={std::abs(uL[0]-uR[0]),std::abs(uL[1]-uR[1]),std::abs(uL[0]-uR[1]),std::abs(uL[1]-uR[0])};
//     int index=std::min_element(difs.begin(),difs.end())-difs.begin();
//          if(index==0) return{uL[0],uR[0]};
//     else if(index==1)
//     return{uL[1],uR[1]};
//     else if(index==2)
//     return{uL[0],uR[1]};
//     else if(index==3)
//     return{uL[1],uR[0]};
//     std::cout<<"spRecon Error: minDif index out\n";
//     return{0,0};
// }
/**
 * @brief 最小差值计算函数（另一种三维数组版本）
 * @param uL 左侧值数组
 * @param uR 右侧值数组
 * @return 包含最小差值的数组
 */
static std::array<real, 2> minDif2(std::array<real, 3> uL, std::array<real, 3> uR)
{
    std::array<real, 3> difs = { std::abs(uL[0] - uR[0]), std::abs(uL[0] - uR[1]), std::abs(uL[1] - uR[0]) };
    int index = std::min_element(difs.begin(), difs.end()) - difs.begin();
    if (index == 0)
        return { uL[0], uR[0] };
    else if (index == 1)
        return { uL[0], uR[1] };
    else if (index == 2)
        return { uL[1], uR[0] };
    std::cout << "spRecon Error: minDif index out\n";
    return { 0, 0 };
}
/**
 * @brief 一维界面中心重构函数
 * @param i 网格点索引
 * @return 重构后的变量值向量
 */
std::vector<real> SpaceDis::recon1DFaceCenter(int i)
{
    std::array<real, 3> primL, primR;
    memcpy(&primL[0], &at(i - 1, 0), nVar * sizeof(real));
    memcpy(&primR[0], &at(i, 0), nVar * sizeof(real));
    eigensystemEuler1D eig = eigensystemEuler1D(primL, primR);
    std::array<real, 5> q1L, q2L, q3L, q1R, q2R, q3R;
    for (int j = i - 3; j < i + 3; j++) {
        enum { R,
            U,
            P };
        auto charTemp = eig.primToChar({ at(j, R), at(j, U), at(j, P) });

        int iLocal = j - i + 3;
        if (iLocal < 5) {
            q1L[iLocal] = charTemp[0];
            q2L[iLocal] = charTemp[1];
            q3L[iLocal] = charTemp[2];
        }

        iLocal = i + 2 - j;
        if (iLocal < 5) {
            q1R[iLocal] = charTemp[0];
            q2R[iLocal] = charTemp[1];
            q3R[iLocal] = charTemp[2];
        }
    }
    // auto start = std::chrono::steady_clock::now();
    auto Q1LL = inter5(q1L);
    auto Q1RR = inter5(q1R);
    auto Q2LL = inter5(q2L);
    auto Q2RR = inter5(q2R);
    auto Q3LL = inter5(q3L);
    auto Q3RR = inter5(q3R);

    // auto Q1LLTHINC=THINC1(q1L[1],q1L[2],q1L[3]);
    // auto Q1RRTHINC=THINC1(q1R[1],q1R[2],q1R[3]);
    // auto Q2LLTHINC=THINC1(q2L[1],q2L[2],q2L[3]);
    // auto Q2RRTHINC=THINC1(q2R[1],q2R[2],q2R[3]);
    // auto Q3LLTHINC=THINC1(q3L[1],q3L[2],q3L[3]);
    // auto Q3RRTHINC=THINC1(q3R[1],q3R[2],q3R[3]);

    // std::array<real,2> Q1,Q2,Q3;
    // Q1=minDif((std::array<real,2>){Q1LL,Q1LLTHINC},(std::array<real,2>){Q1RR,Q1RRTHINC});
    // Q2=minDif((std::array<real,2>){Q2LL,Q2LLTHINC},(std::array<real,2>){Q2RR,Q2RRTHINC});
    // Q3=minDif((std::array<real,2>){Q3LL,Q3LLTHINC},(std::array<real,2>){Q3RR,Q3RRTHINC});
    // Q1LL=Q1[0];Q1RR=Q1[1];
    // Q2LL=Q2[0];Q2RR=Q2[1];
    // Q3LL=Q3[0];Q3RR=Q3[1];
    // auto stop = std::chrono::steady_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
    // timep+=duration;

    auto resTempL = eig.charToPrim({ Q1LL, Q2LL, Q3LL });
    auto resTempR = eig.charToPrim({ Q1RR, Q2RR, Q3RR });
    return { resTempL[0], resTempL[1], resTempL[2], resTempR[0], resTempR[1], resTempR[2] };
}

/**
 * @brief 二维界面中心重构函数
 * @param i 网格点索引
 * @return 重构后的变量值向量
 */
std::vector<real> SpaceDis::recon2DFaceCenter(int i)
{
    // 定义左右两侧的原始变量数组，用于特征变换系统计算
    std::array<real, 4> primL, primR;
    // 复制左侧相邻单元格的变量值到primL
    memcpy(&primL[0], &at(i - 1, 0), nVar * sizeof(real));
    // 复制当前单元格的变量值到primR
    memcpy(&primR[0], &at(i, 0), nVar * sizeof(real));
    
    // 基于左右两侧单元格的变量值和法向量norm构造欧拉方程的特征系统
    eigensystemEuler2D eig = eigensystemEuler2D(primL, primR, norm);
    
    // 定义用于左右两侧特征变量插值的数组，每个数组存储5个点的数据
    std::array<real, 5> q1L, q2L, q3L, q4L, q1R, q2R, q3R, q4R;
    
    // 遍历从i-3到i+2的单元格，提取特征变量用于插值计算
    for (int j = i - 3; j < i + 3; j++) {
        // 定义变量索引枚举，方便访问不同物理量
        enum { R,     // 密度
            U,        // x方向速度
            V,        // y方向速度
            P };      // 压力
        
        // 将物理变量转换为特征变量
        auto charTemp = eig.primToChar({ at(j, R), at(j, U), at(j, V), at(j, P) });

        // 计算左侧相关索引（相对于当前面i）
        int iLocal = j - i + 3;
        // 如果索引有效，则将特征变量存储到左侧数组中
        if (iLocal < 5) {
            q1L[iLocal] = charTemp[0];  // 第一个特征变量
            q2L[iLocal] = charTemp[1];  // 第二个特征变量
            q3L[iLocal] = charTemp[2];  // 第三个特征变量
            q4L[iLocal] = charTemp[3];  // 第四个特征变量
        }

        // 计算右侧相关索引（相对于当前面i，镜像处理）
        iLocal = i + 2 - j;
        // 如果索引有效，则将特征变量存储到右侧数组中
        if (iLocal < 5) {
            q1R[iLocal] = charTemp[0];  // 第一个特征变量
            q2R[iLocal] = charTemp[1];  // 第二个特征变量
            q3R[iLocal] = charTemp[2];  // 第三个特征变量
            q4R[iLocal] = charTemp[3];  // 第四个特征变量
        }
    }
    
    // 声明标志变量，虽然在此段代码中声明了但并未使用
    bool flagL1, flagR1, flagL2, flagR2, flagL3, flagR3, flagL4, flagR4;
    
    // 记录插值计算开始时间
    auto start = std::chrono::steady_clock::now();
    
    // 对左右两侧的四个特征变量分别进行5点插值计算
    auto Q1LL = inter5(q1L);  // 左侧第一个特征变量的插值结果
    auto Q1RR = inter5(q1R);  // 右侧第一个特征变量的插值结果
    auto Q2LL = inter5(q2L);  // 左侧第二个特征变量的插值结果
    auto Q2RR = inter5(q2R);  // 右侧第二个特征变量的插值结果
    auto Q3LL = inter5(q3L);  // 左侧第三个特征变量的插值结果
    auto Q3RR = inter5(q3R);  // 右侧第三个特征变量的插值结果
    auto Q4LL = inter5(q4L);  // 左侧第四个特征变量的插值结果
    auto Q4RR = inter5(q4R);  // 右侧第四个特征变量的插值结果
    
    // 记录插值计算结束时间
    auto stop = std::chrono::steady_clock::now();
    // 计算插值计算耗时并累加到timep变量中
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
    timep += duration;

    // 将左侧插值得到的特征变量转换回物理变量
    auto resTempL = eig.charToPrim({ Q1LL, Q2LL, Q3LL, Q4LL });
    // 将右侧插值得到的特征变量转换回物理变量
    auto resTempR = eig.charToPrim({ Q1RR, Q2RR, Q3RR, Q4RR });
    
    // 返回左右两侧的物理变量值，共8个值（每侧4个变量：密度、x速度、y速度、压力）
    return { resTempL[0], resTempL[1], resTempL[2], resTempL[3], resTempR[0], resTempR[1], resTempR[2], resTempR[3] };
}

/**
 * @brief 三维界面中心重构函数（任务 5 新增；镜像 recon2DFaceCenter，10 值返回）
 * @param i 网格点索引（面索引）
 * @return 左右两侧各 5 个原始变量：{rho,u,v,w,p}_L, {rho,u,v,w,p}_R
 *
 * 与 2D 版结构一致：以左右两格状态做 Roe 平均构造特征系统，
 * 在 6 点模板（j=i-3..i+2）上把 5 个特征场分别用 inter5 格式插值，
 * 再变换回原始量（格式本体由 inter5 指针决定，与维度无关）。
 */
std::vector<real> SpaceDis::recon3DFaceCenter(int i)
{
    std::array<real, 5> primL, primR;
    memcpy(&primL[0], &at(i - 1, 0), nVar * sizeof(real));
    memcpy(&primR[0], &at(i, 0), nVar * sizeof(real));

    eigensystemEuler3D eig = eigensystemEuler3D(primL, primR, norm);

    std::array<real, 5> q1L,q2L,q3L,q4L,q5L,q1R,q2R,q3R,q4R,q5R;
    for (int j = i - 3; j < i + 3; j++)
    {
        enum { R, U, V, W, P };
        auto charTemp = eig.primToChar({ at(j, R), at(j, U), at(j, V), at(j, W), at(j, P) });

        int iLocal = j - i + 3;
        if (iLocal < 5) {
            q1L[iLocal]=charTemp[0]; q2L[iLocal]=charTemp[1]; q3L[iLocal]=charTemp[2];
            q4L[iLocal]=charTemp[3]; q5L[iLocal]=charTemp[4];
        }
        iLocal = i + 2 - j;
        if (iLocal < 5) {
            q1R[iLocal]=charTemp[0]; q2R[iLocal]=charTemp[1]; q3R[iLocal]=charTemp[2];
            q4R[iLocal]=charTemp[3]; q5R[iLocal]=charTemp[4];
        }
    }

    auto Q1LL = inter5(q1L); auto Q1RR = inter5(q1R);
    auto Q2LL = inter5(q2L); auto Q2RR = inter5(q2R);
    auto Q3LL = inter5(q3L); auto Q3RR = inter5(q3R);
    auto Q4LL = inter5(q4L); auto Q4RR = inter5(q4R);
    auto Q5LL = inter5(q5L); auto Q5RR = inter5(q5R);

    auto resTempL = eig.charToPrim({ Q1LL, Q2LL, Q3LL, Q4LL, Q5LL });
    auto resTempR = eig.charToPrim({ Q1RR, Q2RR, Q3RR, Q4RR, Q5RR });

    return { resTempL[0],resTempL[1],resTempL[2],resTempL[3],resTempL[4],
             resTempR[0],resTempR[1],resTempR[2],resTempR[3],resTempR[4] };
}

/**
 * @brief 二维特征重构左侧函数
 * @param i 网格点索引
 * @return 重构后的变量值向量
 */
std::vector<real> SpaceDis::reconLChar2D(int i)
{
    assert(info->dim == 2);
    assert(info->eqType == EULER);
    enum { R,
        U,
        V,
        P };
    eigensystemEuler2D eig = eigensystemEuler2D({ at(i - 1, R), at(i - 1, U), at(i - 1, V), at(i - 1, P) }, norm);
    std::array<real, 5> q1, q2, q3, q4;
    for (int j = i - 3; j < i + 2; j++) {
        int iLocal = j - i + 3;
        auto charTemp = eig.primToChar({ at(j, R), at(j, U), at(j, V), at(j, P) });
        q1[iLocal] = charTemp[0];
        q2[iLocal] = charTemp[1];
        q3[iLocal] = charTemp[2];
        q4[iLocal] = charTemp[3];
    }
    auto start = std::chrono::high_resolution_clock::now();
    real Q1 = (*this->inter5)(q1);
    real Q2 = (*this->inter5Positive)(q2);
    real Q3 = (*this->inter5)(q3);
    real Q4 = (*this->inter5)(q4);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
    timep += duration;
    auto resTemp = eig.charToPrim({ Q1, Q2, Q3, Q4 });
    return { resTemp[0], resTemp[1], resTemp[2], resTemp[3] };
}

/**
 * @brief 二维特征重构右侧函数
 * @param i 网格点索引
 * @return 重构后的变量值向量
 */
std::vector<real> SpaceDis::reconRChar2D(int i)
{
    assert(info->dim == 2);

    enum { R,
        U,
        V,
        P };
    eigensystemEuler2D eig = eigensystemEuler2D({ at(i, R), at(i, U), at(i, V), at(i, P) }, norm);
    std::array<real, 5> q1, q2, q3, q4;
    for (int j = i + 2; j > i - 3; j--) {
        int iLocal = i + 2 - j;
        auto charTemp = eig.primToChar({ at(j, R), at(j, U), at(j, V), at(j, P) });
        q1[iLocal] = charTemp[0];
        q2[iLocal] = charTemp[1];
        q3[iLocal] = charTemp[2];
        q4[iLocal] = charTemp[3];
    }
    auto start = std::chrono::high_resolution_clock::now();
    real Q1 = (*this->inter5)(q1);
    real Q2 = (*this->inter5Positive)(q2);
    real Q3 = (*this->inter5)(q3);
    real Q4 = (*this->inter5)(q4);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
    timep += duration;
    auto resTemp = eig.charToPrim({ Q1, Q2, Q3, Q4 });
    return { resTemp[0], resTemp[1], resTemp[2], resTemp[3] };
}

// std::vector<real> SpaceDis::reconLChar2D(int i)
// {
//     assert(info->dim==2);
//     assert(info->eqType==EULER);
//     enum{R,U,V,P};
//     double rhoAvg,uAvg,vAvg,HAvg,cAvg,VnAvg,q_2Avg,coef1,coef2;
//     std::array<real,2> H;
//     real rl=at(i-1,R),ul=at(i-1,U),vl=at(i-1,V),pl=at(i-1,P);
//     real rr=at(i,R),ur=at(i,U),vr=at(i,V),pr=at(i,P);
//     enum{LL,RR};
//     H[LL]=(ul*ul+vl*vl)/2+pl/rl*GAMMA/(GAMMA-1);
//     H[RR]=(ur*ur+vr*vr)/2+pr/rr*GAMMA/(GAMMA-1);
//     coef1=std::sqrt(rl);
//     coef2=std::sqrt(rr);
//     real divisor=1.0/(std::sqrt(rl)+std::sqrt(rr));
//     real rRef=std::sqrt(rl*rr);
//     uAvg=(coef1*ul+coef2*ur)*divisor;
//     vAvg=(coef1*vl+coef2*vr)*divisor;
//     HAvg=(coef1*H[LL]+coef2*H[RR])*divisor;
//     q_2Avg=(uAvg*uAvg+vAvg*vAvg)*0.5;
//     real cRef=std::sqrt((GAMMA-1)*(HAvg-q_2Avg));
//     std::array<real,5> q1,q2,q3,q4;
//     for(int j=i-3;j<i+2;j++)
//     {
//         int iLocal=j-i+3;
//         real Vn=(norm[1]>norm[0])?at(j,V):at(j,U);
//         q1[iLocal]=(norm[1]>norm[0])?at(j,U):at(j,V);
//         q2[iLocal]=at(j,R)-at(j,P)/(cRef*cRef);
//         q3[iLocal]=Vn+at(j,P)/(cRef*rRef);
//         q4[iLocal]=-Vn+at(j,P)/(cRef*rRef);
//     }
//     auto start = std::chrono::high_resolution_clock::now();
//     real Q1=(*this->inter5)(q1);
//     real Q2=(*this->inter5Positive)(q2);
//     real Q3=(*this->inter5)(q3);
//     real Q4=(*this->inter5)(q4);
//     auto stop = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
//     timep+=duration;
//     std::vector<real> res={Q2+(Q3+Q4)*0.5*rRef/cRef
//             ,(norm[1]>norm[0]? Q1 : (Q3-Q4)*0.5 )
//             ,(norm[1]>norm[0]? (Q3-Q4)*0.5  : Q1)
//             ,(Q3+Q4)*(cRef*rRef)*0.5};
//     return res;
// }

// std::vector<real> SpaceDis::reconRChar2D(int i)
// {
//     assert(info->dim==2);
//     enum{R,U,V,P};

//     double rhoAvg,uAvg,vAvg,HAvg,cAvg,VnAvg,q_2Avg,coef1,coef2;
//     std::array<real,2> H;
//     real rl=at(i-1,R),ul=at(i-1,U),vl=at(i-1,V),pl=at(i-1,P);
//     real rr=at(i,R),ur=at(i,U),vr=at(i,V),pr=at(i,P);
//     enum{LL,RR};
//     H[LL]=(ul*ul+vl*vl)/2+pl/rl*GAMMA/(GAMMA-1);
//     H[RR]=(ur*ur+vr*vr)/2+pr/rr*GAMMA/(GAMMA-1);
//     coef1=std::sqrt(rl);
//     coef2=std::sqrt(rr);
//     real divisor=1.0/(std::sqrt(rl)+std::sqrt(rr));
//     real rRef=std::sqrt(rl*rr);
//     uAvg=(coef1*ul+coef2*ur)*divisor;
//     vAvg=(coef1*vl+coef2*vr)*divisor;
//     HAvg=(coef1*H[LL]+coef2*H[RR])*divisor;
//     q_2Avg=(uAvg*uAvg+vAvg*vAvg)*0.5;
//     real cRef=std::sqrt((GAMMA-1)*(HAvg-q_2Avg));

//     std::array<real,5> q1,q2,q3,q4;
//     for(int j=i+2;j>i-3;j--)
//     {
//         int iLocal=i+2-j;
//         real Vn=(norm[1]>norm[0])?at(j,V):at(j,U);
//         q1[iLocal]=(norm[1]>norm[0])?at(j,U):at(j,V);
//         q2[iLocal]=at(j,R)-at(j,P)/(cRef*cRef);
//         q3[iLocal]=Vn+at(j,P)/(cRef*rRef);
//         q4[iLocal]=-Vn+at(j,P)/(cRef*rRef);
//     }
//     auto start = std::chrono::high_resolution_clock::now();
//     real Q1=(*this->inter5)(q1);
//     real Q2=(*this->inter5Positive)(q2);
//     real Q3=(*this->inter5)(q3);
//     real Q4=(*this->inter5)(q4);
//     auto stop = std::chrono::high_resolution_clock::now();
//     auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
//     timep+=duration;
//     std::vector<real> res={Q2+(Q3+Q4)*0.5*rRef/cRef
//             ,(norm[1]>norm[0]? Q1 : (Q3-Q4)*0.5 )
//             ,(norm[1]>norm[0]? (Q3-Q4)*0.5  : Q1)
//             ,(Q3+Q4)*(cRef*rRef)*0.5};
//     return res;
// }

/**
 * @brief 右侧原始变量重构函数
 * @param i 网格点索引
 * @return 重构后的变量值向量
 */
std::vector<real> SpaceDis::reconRprim(int i)
{
    std::vector<real> res(nPrim);
    for (int ivar = 0; ivar < nPrim; ivar++) {
        std::array iprims = { at(i + 2, ivar), at(i + 1, ivar), at(i, ivar), at(i - 1, ivar), at(i - 2, ivar) };
        res[ivar] = (*this->inter5)(iprims);
    }
    return res;
}

/**
 * @brief 左侧原始变量重构函数
 * @param i 网格点索引
 * @return 重构后的变量值向量
 */
std::vector<real> SpaceDis::reconLprim(int i)
{
    std::vector<real> res(nPrim);
    for (int ivar = 0; ivar < nPrim; ivar++) {
        std::array iprims = { at(i - 3, ivar), at(i - 2, ivar), at(i - 1, ivar), at(i, ivar), at(i + 1, ivar) };
        res[ivar] = (*this->inter5)(iprims);
    }
    return res;
}
