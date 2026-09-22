/**
 * @file block.hpp
 * @brief 定义Block类，用于处理计算网格块的相关数据和操作
 */

#pragma once
#include"data.hpp"
#include <cgnslib.h>

/**
 * @brief Block类，用于管理计算域中的网格块
 * 
 * 该类包含网格顶点坐标、单元中心坐标、网格间距等信息，
 * 并提供了相关的网格操作方法。
 */
class Block
{
    public:
    /**
     * @brief 将网格数据输出到CGNS文件
     */
    void outputCgns();
    
    /**
     * @brief 重载括号运算符，获取指定单元的坐标
     * @param ic 单元索引
     * @param idim 坐标维度索引
     * @return 指定单元在指定维度上的坐标值
     */
    real operator()(int,int);
    
    /**
     * @brief 获取内部单元的最大索引数组
     * @return 内部单元最大索引数组
     */
    std::array<int,3> getICMax();
    
    /**
     * @brief 获取网格节点最大索引数组
     * @return 网格节点最大索引数组
     */
    std::array<int,3> getIMax();
    
    /**
     * @brief 获取指定方向上的最小网格间距
     * @param i 网格节点索引
     * @return 最小网格间距
     */
    real getMinDh(int i);
    
    /**
     * @brief 获取网格维度
     * @return 网格维度 (1/2/3D)
     */
    int getDim();
    
    /**
     * @brief 获取单元中心坐标
     * @param idim 坐标维度索引
     * @return 单元中心坐标向量
     */
    std::vector<real> getCellCoor(int idim);
    
    /**
     * @brief 获取顶点坐标
     * @param idim 坐标维度索引
     * @return 顶点坐标向量
     */
    std::vector<real> getVertexCoor(int idim);
    
    /**
     * @brief 获取单元间距
     * @param i 单元索引
     * @return 单元间距向量
     */
    std::vector<real> getCellInterval(int);


    protected:
    friend class Initializer;
    Data coorVer;               ///< 顶点坐标数据
    Data coorCel;               ///< 单元中心坐标数据
    Data intervalCel;           ///< 单元间距数据
    int nVer;                   ///< 顶点数量
    int nCel;                   ///< 单元数量
    bool inited;                ///< 初始化标志
    int dim;                    ///< 网格维度
    std::array<int,3> iMax;     ///< 网格节点最大索引
    std::array<int,3> icMax;    ///< 内部单元最大索引

};