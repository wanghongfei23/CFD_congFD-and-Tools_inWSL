/**
 * @file cgnsio.cpp
 * @brief CGNS输入输出类的实现文件
 */

#include "cgnsio.hpp"

/**
 * @brief 将网格数据输出到CGNS文件
 * @param block 网格块指针
 * @param info 信息对象指针
 * 
 * 该函数创建一个CGNS格式的网格文件，包含网格顶点坐标信息。
 * 支持1D/2D/3D结构化网格。
 */
void CgnsIO::BlockCgnsOutput(Block* block,Info* info)
{
   {                                                                     // 获取网格维度
   int dim=block->getDim();                                              // 网格维度(1/2/3D)
   
   auto iMax=block->getIMax();                                           // 获取网格节点最大索引数组
   cgsize_t isize[3][dim]={0};                                           // CGNS网格尺寸数组
   auto filename=info->filename();                                       // 获取输出文件名
   int index_file,icelldim,iphysdim,index_base;                          // CGNS文件相关索引变量
   int index_zone,index_coord;                                           // CGNS区域和坐标索引
   char basename[33],zonename[33];                                       // 基础名称和区域名称

   if (cg_open(filename.c_str(),CG_MODE_WRITE,&index_file)) cg_error_exit(); // 创建并打开CGNS文件
   strcpy(basename,"Base");                                              // 设置基础名称
   icelldim=dim;                                                         // 单元维度
   iphysdim=dim;                                                         // 物理维度
   if (cg_base_write(index_file,basename,icelldim,iphysdim,&index_base)) cg_error_exit(); // 写入基础信息
   strcpy(zonename,"Zone  1");                                           // 设置区域名称
   for (int idim = 0; idim < dim; idim++)                                // 循环设置各维度的尺寸信息
   {
        isize[0][idim]=iMax[idim];                                      // 顶点数
        isize[1][idim]=iMax[idim]-1;                                    // 单元数
        isize[2][idim]=0;                                               // 边界顶点数(结构化网格始终为0)
   }
   

   if (cg_zone_write(index_file,index_base,zonename,*isize,CGNS_ENUMV(Structured),&index_zone)) cg_error_exit(); // 写入区域信息
   for (int idim = 0; idim < dim; idim++)                               // 循环写入各维度坐标
   {
        std::string name=(idim==0)? "CoordinateX":(idim==1)? "CoordinateY":"CoordinateZ"; // 根据维度确定坐标名称
        auto x=block->getVertexCoor(idim);                               // 获取指定维度的顶点坐标
        if (cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble),
        name.c_str(),
        x.data(),&index_coord)) cg_error_exit();                         // 写入坐标数据
    }
    
   
   cg_close(index_file);}                                                // 关闭CGNS文件
}

/**
 * @brief 将解数据输出到CGNS文件
 * @param data 数据指针
 * @param info 信息对象指针
 * 
 * 该函数将计算得到的解数据写入已存在的CGNS文件中。
 */
void CgnsIO::solCgnsOutput(Data* data,Info* info)
{

    int index_file;                                                       // CGNS文件索引
    auto filename=info->filename();                                       // 获取文件名
    if (cg_open(filename.c_str(),CG_MODE_MODIFY,&index_file)) cg_error_exit(); // 以修改模式打开已存在的CGNS文件
   int index_base,index_zone,index_flow,index_field;                      // CGNS各种索引变量
    index_base=1;                                                        // 基础索引(假设只有一个基础)
    index_zone=1;                                                        // 区域索引(假设只有一个区域)
   std::string solname="flowSolution";                                   // 设置流动解名称
    if (cg_sol_write(index_file,index_base,index_zone,solname.c_str(),CGNS_ENUMV(CellCenter),&index_flow)) cg_error_exit(); // 写入解节点信息
   for(int ivar=0;ivar<data->getNVar();ivar++)                           // 循环写入各个变量
   {
      auto varr=data->getIVar(ivar);                                     // 获取第ivar个变量的数据
      if (cg_field_write(index_file,index_base,index_zone,index_flow,
        CGNS_ENUMV(RealDouble),data->varName[ivar].c_str(),varr.data(),&index_field)) cg_error_exit(); // 写入变量场数据
   }

    if (cg_close(index_file)) cg_error_exit();;                          // 关闭CGNS文件
}
