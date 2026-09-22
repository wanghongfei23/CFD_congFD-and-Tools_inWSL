/**
 * @file cgnsio.hpp
 * @brief 定义CGNS输入输出类，用于网格和解数据的CGNS格式读写
 */

#pragma once
#include "macro.hpp"
#include "cgnslib.h"
#include "block.hpp"
#include "data.hpp"
#include "info.hpp"

/**
 * @brief CGNS输入输出类
 * 
 * 该类提供了将网格数据和求解结果写入CGNS格式文件的功能。
 */
class CgnsIO
{
    public:
    /**
     * @brief 将网格数据输出到CGNS文件
     * @param block 网格块指针
     * @param info 信息对象指针
     */
    void BlockCgnsOutput(Block* block,Info* info);
    
    /**
     * @brief 将解数据输出到CGNS文件
     * @param data 数据指针
     * @param info 信息对象指针
     */
    void solCgnsOutput(Data* data,Info* info);
    
    /**
     * @brief 一维解数据输出
     * @param info 信息对象指针
     */
    void oneDsolOutput(Info* info);
    
    private:
};