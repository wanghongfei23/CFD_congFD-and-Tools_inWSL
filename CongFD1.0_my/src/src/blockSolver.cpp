/**
 * @file blockSolver.cpp
 * @brief BlockSolver类的实现文件
 */

#include "blockSolver.hpp"
#include "diagnostics.hpp"
#include <fstream>
#include <iomanip>
#include <filesystem>

/**
 * @brief 追加一行诊断输出（时间与 K/Ω/U2/Tvar/Thvar 五列）到当前目录 diagnostics.txt（仅三维调用）
 */
static void appendDiagnostics(Info* info, Data* prim, Block* block)
{
    real K,Omega,U2,Tvar,Thvar;
    calcDiagnostics(prim, block->getICMax(), info->interval, K, Omega, U2, Tvar, Thvar);
    bool fresh=!std::filesystem::exists("diagnostics.txt");
    std::ofstream diag("diagnostics.txt", std::ios::app);
    if(fresh) diag << "# t K Omega U2 Tvar Thvar\n";
    diag << std::scientific << std::setprecision(16) << info->t << " " << K << " " << Omega
         << " " << U2 << " " << Tvar << " " << Thvar << "\n";
}

/**
 * @brief BlockSolver类的默认构造函数
 * 
 * 初始化求解器所需的各个组件，包括：
 * - Info: 求解器参数信息
 * - Block: 网格数据块
 * - Initializer: 初始化器
 * - Equation: 方程求解器
 * - Bnds: 边界条件处理器
 * - SpDistributor: 空间离散分配器
 * - SourceTerm: 源项计算器
 * 
 * 并调用初始化器对这些组件进行初始化。
 */
BlockSolver::BlockSolver()
{
   info=new Info();
   block=new Block();
   initer=new Initializer(info);
   eqn=new Equation();
   bnds=new Bnds();
   spDis=new SpDistributor();

   initer->initUniformBlock(block);
   initer->initEqution(eqn,block);
   initer->initBnds(bnds,eqn,block->getICMax(),block);
   initer->initSpDistributor(spDis,eqn,block,bnds);
   sourceTerm=new SourceTerm(eqn->getPrim(),eqn->getRhs(),info);
   viscousTerm=new ViscousTerm(eqn->getPrim(),eqn->getRhs(),info);

   cons=eqn->getCons();
   rhs=eqn->getRhs();
}

/**
 * @brief BlockSolver类的带参构造函数
 * 
 * 使用传入的Info对象初始化求解器，其余初始化过程与默认构造函数相同。
 * 
 * @param info_ 指向Info对象的指针，包含求解器参数信息
 */
BlockSolver::BlockSolver(Info* info_)
{
    info=info_;
    block=new Block();
    initer=new Initializer(info);
    eqn=new Equation();
    bnds=new Bnds();
    spDis=new SpDistributor();
    

    initer->initUniformBlock(block);
    initer->initEqution(eqn,block);
    initer->initBnds(bnds,eqn,block->getICMax(),block);
    initer->initSpDistributor(spDis,eqn,block,bnds);
    sourceTerm=new SourceTerm(eqn->getPrim(),eqn->getRhs(),info);
    viscousTerm=new ViscousTerm(eqn->getPrim(),eqn->getRhs(),info);

    cons=eqn->getCons();
    rhs=eqn->getRhs();
}

/**
 * @brief 使用三阶SSP Runge-Kutta方法进行时间推进
 * 
 * 实现三阶段的三阶SSP Runge-Kutta时间积分方案：
 * 1. 第一阶段：计算临时解
 * 2. 第二阶段：基于第一阶段结果计算
 * 3. 第三阶段：得到最终解
 * 
 * @param dt 时间步长
 */
void BlockSolver::RK3_SSP(real dt)
{
    
    Data tempdata(*cons);
    int n=cons->size();
    

    // third order RK
    // stage 1
    // 将右端项(rhs)重置为零，为新一轮计算做准备
    rhs->setZeros();
    // 将守恒变量转换为原始变量（密度、速度、压力等），用于后续的通量计算
    eqn->consToPrim();
    // 更新所有边界条件，确保边界上的值反映最新的内部流场状态
    bnds->update();
    // 计算空间离散项，即对流项的数值导数，贡献到右端项中
    spDis->rhsSolve();
    // 计算粘性项（黏性应力＋热传导的散度）并累加进 rhs（任务 4；与源项并列）
    viscousTerm->calViscous();
    // 计算物理源项（如体积力、化学反应等）并将其贡献添加到右端项中
    sourceTerm->calSource();
    // cgnsIO.BlockCgnsOutput(block,info);
    // cgnsIO.solCgnsOutput(rhs,info);
    #pragma omp parallel for
    for(int i=0;i<n;i++)
    {
        (*cons)[i]=tempdata[i]-dt*(*rhs)[i];
    }
    info->t+=dt;

    //stage 2
    // 将右端项重置为零，开始新的时间步计算
    rhs->setZeros();
    // 将守恒变量转换为原始变量，用于计算通量
    eqn->consToPrim();
    // 更新边界条件，确保边界值是最新的
    bnds->update();
    // 计算空间离散项，得到对流项的数值导数贡献
    spDis->rhsSolve();
    // 计算粘性项并累加（任务 4）
    viscousTerm->calViscous();
    // 计算源项（如体积力、能量源等）并添加到右端项中
    sourceTerm->calSource();
    #pragma omp parallel for
    for(int i=0;i<n;i++)
    {
        (*cons)[i]=0.75*tempdata[i]-0.25*dt*(*rhs)[i]+0.25*(*cons)[i];
    }
    info->t-=dt/2;

    //stage 3
    rhs->setZeros();
    eqn->consToPrim();
    bnds->update();
    spDis->rhsSolve();
    viscousTerm->calViscous();   // 任务 4：粘性项累加
    sourceTerm->calSource();
    #pragma omp parallel for
    for(int i=0;i<n;i++)
    {
        (*cons)[i]=1.0/3.0*tempdata[i]-2.0/3.0*dt*(*rhs)[i]+2.0/3.0*(*cons)[i];
    }
    info->t-=dt/2;
}


/**
 * @brief 使用DTS Euler方法进行时间推进
 * 
 * 实现基于对角化时间步长(Diagonalized Time Stepping)的欧拉时间积分方法，
 * 用于隐式时间推进。
 * 
 * @param dt 时间步长
 */
void BlockSolver::DTS_Euler(real dt)
{
    
    
    int n=cons->size();
    Data tempdata(*cons);
    Data tempRhs(*cons);
    real maxres=0;

    int imStep=0;
    do
    {
        
        auto dtau=calLocalCFL();
        // 将右端项(rhs)重置为零，开始新的隐式迭代计算
        rhs->setZeros();
        // 将当前守恒变量转换为原始变量，用于计算通量和雅可比矩阵
        eqn->consToPrim();
        // 更新边界条件，确保在迭代过程中边界值保持最新
        bnds->update();
        // 计算空间离散项（对流项的数值导数），贡献到右端项中
        spDis->rhsSolve();
        // 计算源项（如粘性、热传导等）并添加到右端项中
        sourceTerm->calSource();
        // cgnsIO.BlockCgnsOutput(block,info);
        // cgnsIO.solCgnsOutput(rhs,info);
        tempRhs.setZeros();
        #pragma omp parallel for
        for(int i=0;i<n;i++)
        {
            real coef=dtau[i]/(dt+dtau[i]);
            tempRhs[i]=(tempdata[i]-(*cons)[i]-dt*(*rhs)[i])*coef;
            (*cons)[i]=(*cons)[i]+tempRhs[i];
        }
        maxres=tempRhs.getLinf(0);
        // 修复：使用传统的iostream方式替代std::format以避免编译错误
        std::cout << "time = " << std::fixed << std::setprecision(4) << info->t 
                  << " dt=" << std::setprecision(10) << dt 
                  << " imStep=" << imStep 
                  << " rhoL2 = " << std::setprecision(4) << maxres << "  \n";
    }while(imStep++<info->maxImplicitStep && maxres>1e-7);


}

/**
 * @brief 执行单步时间推进
 * 
 * 根据设定的时间积分方法执行一步时间推进计算。
 * 
 * @param dt 时间步长
 */
void BlockSolver::solve(real dt)
{
    RK3_SSP(dt);
    info->t+=dt;
}

/**
 * @brief BlockSolver类的析构函数
 * 
 * 释放所有动态分配的内存资源。
 */
BlockSolver::~BlockSolver()
{
    delete info;
    delete block;
    delete initer;
    delete eqn;
    delete bnds;
    delete spDis;
    delete sourceTerm;
    delete viscousTerm;
}

/**
 * @brief 输出网格信息到CGNS文件
 */
void BlockSolver::outputGrid()
{
    cgnsIO.BlockCgnsOutput(block,info);
}

/**
 * @brief 输出守恒变量到CGNS文件
 */
void BlockSolver::outputCons()
{
    cgnsIO.solCgnsOutput(eqn->getCons(),info);
}

/**
 * @brief 输出原始变量到CGNS文件
 * 
 * 首先将守恒变量转换为原始变量，然后输出到CGNS文件。
 */
void BlockSolver::outputPrim()
{
    eqn->consToPrim();
    cgnsIO.solCgnsOutput(eqn->getPrim(),info);
}

/**
 * @brief 基于固定步数的时间步进循环
 * 
 * 按照固定的步数进行时间推进，每隔一定步数输出一次结果。
 */
void BlockSolver::stepsLoop()
{
    for(;info->step<info->endStep;info->step++)
    {
        if(info->step%info->outputInterval==0)
        {
            outputGrid();
            outputPrim();
        }

        real dt=info->dt;
        solve(dt);
        timesteps++;
        // 修复：使用传统的iostream方式替代std::format以避免编译错误
        std::cout << "time = " << std::fixed << std::setprecision(4) << info->t 
                  << "  rhoL2 = " << rhs->getL2(0) << "  \n";
    }
    outputGrid();
    outputPrim();
}

/**
 * @brief 基于CFL条件控制时间步长的求解循环
 * 
 * 此函数实现了基于CFL条件的自适应时间步长求解循环。在每次迭代中：
 * 1. 检查是否需要输出当前结果
 * 2. 计算满足CFL条件的最大时间步长
 * 3. 根据输出间隔调整时间步长
 * 4. 执行一次时间步进求解
 * 5. 更新时间和步数计数器
 * 6. 输出当前求解状态信息
 * 循环持续到达到设定的结束步数为止。
 * 
 * @note 该方法使用可变时间步长，确保数值稳定性
 * @see getTimeIntervalExplicit() 获取满足CFL条件的时间步长
 * @see solve() 执行单步求解
 * @see outputGrid() 输出网格信息
 * @see outputPrim() 输出原始变量
 */
void BlockSolver::stepsLoopCFL()
{
    const real eps=1e-15;                                                     // 定义一个极小值，用于浮点数比较
    while(info->step<info->endStep)                                           // 当前步数小于结束步数时继续循环
    {
        // 当输出时间接近0时，输出网格和流场数据
        if(info->outputT<eps)                                                 // 判断是否需要输出结果（输出时间接近0）
        {
            outputGrid();                                                     // 输出网格信息到文件
            outputPrim();                                                     // 输出原始变量到文件
            if(info->dim==3) appendDiagnostics(info, eqn->getPrim(), block);  // 三维：追加 K/Ω 诊断行（任务 7）
        }

        auto dt = info->viscous ? getTimeIntervalExplicitViscous()
                                : getTimeIntervalExplicit();                  // 两套 dt 路线：粘性档增黏性谱半径项（方案 §4.3）
        // dt=info->dt;
        
        // 如果下一个时间步会超过输出间隔，则调整时间步长以精确匹配输出时刻
        if(info->outputT+dt>info->outputDt)                                   // 判断加上当前步长是否会超过输出时间间隔
        {
            dt=info->outputDt-info->outputT;                                  // 调整时间步长，使其恰好到达输出时刻
            info->outputT=0.0;                                                // 重置输出时间计数器
            info->step++;                                                     // 步数增加1
        }
        else info->outputT+=dt;                                               // 否则更新输出时间计数器

        // 执行求解并计时
        auto start = std::chrono::high_resolution_clock::now();               // 记录求解开始时间
        solve(dt);                                                            // 执行一次时间步长为dt的求解
        auto stop = std::chrono::high_resolution_clock::now();                // 记录求解结束时间
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count(); // 计算求解耗时（毫秒）
        timesss+=duration;                                                    // 累加总求解时间
        timesteps++;                                                          // 时间步数计数器增加

        // 输出当前时间、时间步长和密度Linf范数
        // std::cout<<std::format("time = {:.4f} dt={:.10f}  rhoLinf = {:.4f}  \n",info->t,dt,rhs->getLinf(0)); // 输出当前时间、时间步长和守恒变量第0个分量的Linf范数
        
        // 修复：使用传统的ostringstream方式替代std::format以避免编译错误
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(4) << "time = " << info->t 
            << " dt=" << std::setprecision(10) << dt 
            << "  rhoLinf = " << std::setprecision(4) << rhs->getLinf(0) << "  \n";
        std::string output_str = oss.str();
        std::cout << output_str;

        // // 【王鸿飞】begin辅助工具
        
        // // 将相同信息写入txt文件，使用与CGNS文件相同的命名方式，但不包含时间戳
        // std::string base_filename = info->filename();
        // // 移除时间戳部分（从最后一个" - t="开始的部分）
        // size_t pos = base_filename.rfind(" - t=");
        // if (pos != std::string::npos) {
        //     // 修复：使用构造函数替代直接赋值以避免consteval错误
        //     base_filename = std::string(base_filename.begin(), base_filename.begin() + pos) + " - rholinf.txt";
        // } else {
        //     // 如果没有找到时间戳模式，则直接替换扩展名
        //     pos = base_filename.find(".cgns");
        //     if (pos != std::string::npos) {
        //         base_filename.replace(pos, 5, ".txt");
        //     }
        // }
        
        // std::ofstream log_file(base_filename, std::ios::app);
        // if (log_file.is_open()) {
        //     log_file << output_str;
        //     log_file.close();
        // }
        // // 【王鸿飞】end辅助工具

    }
    // 求解完成后输出最终结果
    outputGrid();                                                             // 求解完成后输出最终网格信息
    outputPrim();                                                             // 求解完成后输出最终原始变量
    if(info->dim==3) appendDiagnostics(info, eqn->getPrim(), block);          // 三维：追加 K/Ω 诊断行（任务 7）
}

/**
 * @brief 测试函数
 * 
 * 执行简单的测试计算并输出结果。
 */
void BlockSolver::Test()
{
    outputGrid();
    RK3_SSP(info->dt);

    //outputCons();
    //outputPrim();
    cgnsIO.solCgnsOutput(rhs,info);

}

/**
 * @brief 计算局部CFL数
 * 
 * 根据当前流动状态和网格信息计算每个网格点的局部CFL数。
 * 
 * @return 包含每个网格点CFL数的向量
 */
std::vector<real> BlockSolver::calLocalCFL()
{
   
    real dt,lambda=0;
    
    Data* prim = eqn->getPrim();
    int n = prim->getN();
    std::vector<real> res(n);
    // 将守恒变量转换为原始变量，以便计算局部CFL数
    eqn->consToPrim();
    // 并行计算每个网格点的局部CFL时间步长
    #pragma omp parallel for reduction(max:lambda)
    for(int i=0;i<n;i++)
    {
        real iLambda;
        real dh=block->getMinDh(i);
        real RT=std::sqrt(GAMMA*(*prim)(i,2)/(*prim)(i,0))+std::abs((*prim)(i,1));
        if(info->dim==2) iLambda=(std::sqrt(GAMMA*(*prim)(i,3)/(*prim)(i,0))+std::abs((*prim)(i,1))+std::abs((*prim)(i,2)))/dh;
        else if(info->dim==1)   iLambda=(std::sqrt(GAMMA*(*prim)(i,2)/(*prim)(i,0))+std::abs((*prim)(i,1)))/dh;
        res[i]=info->implicitCFL/iLambda;
    }
    if(info->eqType!=EULER)
    { std::cout<<"Block Solver error: CFL Loop with not EULER solver\n";}
    return res;
    

}

/**
 * @brief 基于DTS方法的时间步进循环
 * 
 * 使用对角化时间步长(DTS)方法进行时间推进。
 */
void BlockSolver::stepsLoopDTS()
{
    const real eps=1e-15;
    while(info->step<info->endStep)
    {
        if(info->outputT<eps)
        {
            outputGrid();
            outputPrim(); 
        }

        real dt, lambda = 0;
        Data* prim = eqn->getPrim();
        int n = prim->getN();
        // 将守恒变量转换为原始变量，以便计算局部波速
        eqn->consToPrim();
        // 并行计算全局最大特征值（lambda），用于确定稳定的时间步长
        #pragma omp parallel for reduction(max:lambda)
        for(int i=0;i<n;i++)
        {
            real iLambda;
            real dh=block->getMinDh(i);
            real RT=std::sqrt(GAMMA*(*prim)(i,2)/(*prim)(i,0))+std::abs((*prim)(i,1));
            if(info->dim==2) iLambda=(std::sqrt(GAMMA*(*prim)(i,3)/(*prim)(i,0))+std::abs((*prim)(i,1))+std::abs((*prim)(i,2)))/dh;
            else if(info->dim==1)   iLambda=(std::sqrt(GAMMA*(*prim)(i,2)/(*prim)(i,0))+std::abs((*prim)(i,1)))/dh;
            lambda=std::max(lambda,iLambda);
        }
        if(info->eqType!=EULER)
        { std::cout<<"Block Solver error: CFL Loop with not EULER solver\n"; lambda=1;}
        dt=info->CFL/lambda;
        // dt=info->dt;
        if(info->outputT+dt>info->outputDt)
        {
            dt=info->outputDt-info->outputT;
            info->outputT=0.0;
            info->step++;
        }
        else info->outputT+=dt;

        DTS_Euler(dt);
        info->t+=dt;
        //std::cout<<std::format("time = {:.4f} dt={:.10f}  rhoL2 = {:.4f}  \n",info->t,dt,rhs->getL2(0));
    }
    outputGrid();
    outputPrim();
}

/**
 * @brief 计算显式时间积分的时间步长
 * 
 * 根据CFL条件和当前流动状态计算满足稳定条件的最大时间步长。
 * 
 * @return 计算得到的时间步长
 */
real BlockSolver::getTimeIntervalExplicit()
{
    real dt,lambda=0;
    int dim=info->dim;
    if(info->eqType == EULER)
    {
        Data* prim = eqn->getPrim();
        int n = prim->getN();
        // 将守恒变量转换为原始变量，用于计算声速和流动速度
        eqn->consToPrim();
        std::vector<real> lambdas(dim);
        for(int idim=0;idim<dim;idim++)
        {
            real maxLambda=0;
            #pragma omp parallel for reduction(max:maxLambda)
            for(int i=0;i<n;i++)
            {
                real iLambda;
                auto dhs=block->getCellInterval(i);
                real dh=info->constH?info->geth(idim):dhs.at(idim);
                iLambda=(std::sqrt(GAMMA*(*prim)(i,dim+1)/(*prim)(i,0))+std::abs((*prim)(i,idim+1)))/dh;
                maxLambda=std::max(maxLambda,iLambda);
            }
            lambdas[idim]=maxLambda;
        }
        lambda=0;
        for(auto iLambda:lambdas) lambda+=iLambda;

    }
    else
    { std::cout<<"Block Solver error: CFL Loop with not EULER solver\n"; lambda=1;}
    dt=info->CFL/lambda;
    return dt;
}

/**
 * @brief 获取粘性档显式时间间隔（对流＋黏性谱半径之和）
 *
 * 两套 dt 路线（飞哥已定，开发方案 §4.3）：无粘档沿用 getTimeIntervalExplicit（原样）；
 * 粘性档用本函数：λ_d＝对流谱半径＋max(4/3, γ/Pr)·(μ/ρ)/Δx_d²，逐向求和后 dt=CFL/λ。
 */
real BlockSolver::getTimeIntervalExplicitViscous()
{
    real lambda=0;
    int dim=info->dim;
    Data* prim = eqn->getPrim();
    int n = prim->getN();
    eqn->consToPrim();
    real mu=1.0/info->Re;
    real coef=std::max(4.0/3.0, GAMMA/info->Pr);
    std::vector<real> lambdas(dim);
    for(int idim=0;idim<dim;idim++)
    {
        real maxLambda=0;
        #pragma omp parallel for reduction(max:maxLambda)
        for(int i=0;i<n;i++)
        {
            real iLambda;
            auto dhs=block->getCellInterval(i);
            real dh=info->constH?info->geth(idim):dhs.at(idim);
            iLambda=(std::sqrt(GAMMA*(*prim)(i,dim+1)/(*prim)(i,0))+std::abs((*prim)(i,idim+1)))/dh
                    + coef*(mu/(*prim)(i,0))/(dh*dh);
            maxLambda=std::max(maxLambda,iLambda);
        }
        lambdas[idim]=maxLambda;
    }
    for(auto iLambda:lambdas) lambda+=iLambda;
    return info->CFL/lambda;
}