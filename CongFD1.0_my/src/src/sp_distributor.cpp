/**
 * @file sp_distributor.cpp
 * @brief 空间分布器类的实现文件
 */

#include "sp_distributor.hpp"

/**
 * @brief 求解右手端
 * 
 * 对各个维度分别进行空间离散计算，得到右端项。
 * 支持1D/2D/3D问题的并行计算，使用OpenMP进行并行化。
 */
void SpDistributor::rhsSolve()
{
    // 检查维度是否合法（1-3维）
    if (dim <= 0 || dim > 3) {
        std::cout << "SpDistributor error: dim error";
        return;
    }

    // 1维空间离散
    if (dim >= 1) {
        long long timep = 0;
#pragma omp parallel for collapse(2) reduction(+ : timep) schedule(static)
        for (int i = 0; i < iMax[1]; i++)
            for (int j = 0; j < iMax[2]; j++) {
                // std::cout<<std::format("thread num: {} \n",omp_get_thread_num());
                auto oneDBnds = bnds->getOneDBnd(1, i, j);
                auto offsets = calOffset(1, i, j, iMax);
                SpaceDis spDis(iMax[0], prim, rhs, oneDBnds[0], oneDBnds[1], info);

                spDis.setConstNorm({ 1, 0, 0 });
                spDis.setOffset(offsets[0], offsets[1]);
                spDis.setIDim(0);
                spDis.difference();
                timep += spDis.timep;
            }
        timepp += timep;
    }

    // 2维空间离散
    if (dim >= 2) {
        long long timep = 0;
        // 使用OpenMP并行处理，collapse(2)将两个嵌套循环合并为一个并行区域
        // reduction(+ : timep)确保每个线程的时间累加到全局timep变量
        // schedule(static)将迭代均匀分配给各线程
#pragma omp parallel for collapse(2) reduction(+ : timep) schedule(static)
        for (int i = 0; i < iMax[0]; i++)
            for (int j = 0; j < iMax[2]; j++) {
                // 获取当前(i,j)位置上第二维的边界条件
                auto oneDBnds = bnds->getOneDBnd(2, i, j);
                // 计算当前(i,j)位置上的偏移量，用于确定数据访问位置
                auto offsets = calOffset(2, i, j, iMax);
                // 创建空间离散对象，传入网格点数、原始数据、右端项、左右边界条件和信息对象
                SpaceDis spDis(iMax[1], prim, rhs, oneDBnds[0], oneDBnds[1], info);

                // 设置法向量为Y方向单位向量(0,1,0)，用于第二维计算
                spDis.setConstNorm({ 0, 1, 0 });
                // 设置当前处理的维度索引为1（对应第二维）
                spDis.setIDim(1);
                // 设置数据访问的起始偏移量和步长偏移量
                spDis.setOffset(offsets[0], offsets[1]);
                // 执行差分计算，得到该维度对右端项的贡献
                spDis.difference();
                // 累加本次计算所用时间
                timep += spDis.timep;
            }
        // 将本轮计算时间累加到总时间
        timepp += timep;
    }

    if (dim >= 3) {
        long long timep = 0;
#pragma omp parallel for collapse(2) reduction(+ : timep) schedule(static)
        for (int i = 0; i < iMax[0]; i++)
            for (int j = 0; j < iMax[1]; j++) {
                auto oneDBnds = bnds->getOneDBnd(3, i, j);
                auto offsets = calOffset(3, i, j, iMax);
                SpaceDis spDis(iMax[2], prim, rhs, oneDBnds[0], oneDBnds[1], info);

                spDis.setConstNorm({ 0, 0, 1 });
                spDis.setIDim(2);
                spDis.setOffset(offsets[0], offsets[1]);

                spDis.difference();
                timep += spDis.timep;
            }
        timepp += timep;
    }
}