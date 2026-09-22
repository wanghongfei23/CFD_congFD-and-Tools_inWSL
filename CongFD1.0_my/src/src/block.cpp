/**
 * @file block.cpp
 * @brief Block类的实现文件
 */

#include"block.hpp"


/**
 * @brief 将网格数据输出到CGNS文件
 * 
 * 该函数创建一个CGNS格式的网格文件，包含网格顶点坐标信息。
 * 支持1D/2D/3D结构化网格。
 */
void Block::outputCgns()
{
   cgsize_t isize[3][dim];
   int ni,nj,nk,i,j,k;
   int index_file,icelldim,iphysdim,index_base;
   int index_zone,index_coord;
   char basename[33],zonename[33];

/* create gridpoints for simple example: */
   std::vector<real> x;
   x.reserve(nVer);

/* WRITE X, Y, Z GRID POINTS TO CGNS FILE */
/* open CGNS file for write */
   if (cg_open("grid_c.cgns",CG_MODE_WRITE,&index_file)) cg_error_exit();
/* create base (user can give any name) */
   strcpy(basename,"Base");
   icelldim=dim;
   iphysdim=dim;
   cg_base_write(index_file,basename,icelldim,iphysdim,&index_base);
/* define zone name (user can give any name) */
   strcpy(zonename,"Zone  1");
   for (int idim = 0; idim < dim; idim++)
   {
        /* vertex size */
        isize[0][idim]=iMax[idim];
        /* cell size */
        isize[1][idim]=isize[0][idim]-1;
        /* boundary vertex size (always zero for structured grids) */
        isize[2][idim]=0;
   }
   

/* create zone */
   cg_zone_write(index_file,index_base,zonename,*isize,CGNS_ENUMV(Structured),&index_zone);
/* write grid coordinates (user must use SIDS-standard names here) */
    for (int idim = 0; idim < dim; idim++)
    {
        std::string name=(idim==0)? "CoordinateX":(idim==1)? "CoordinateY":"CoordinateZ";
        for (int i=0 ; i<nVer ; i++ ) x[i]=coorVer(i,idim);
        cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble),
        name.c_str(),
        x.data(),&index_coord);
    }
    
   
/* close CGNS file */
   cg_close(index_file);
}


/**
 * @brief 重载括号运算符，获取指定单元的坐标
 * @param ic 单元索引
 * @param idim 坐标维度索引
 * @return 指定单元在指定维度上的坐标值
 */
real Block::operator()(int ic,int idim)
{
    return coorCel(ic,idim);
}


/**
 * @brief 获取内部单元的最大索引数组
 * @return 内部单元最大索引数组
 */
std::array<int,3> Block::getICMax()
{
   return icMax;
}

/**
 * @brief 获取网格节点最大索引数组
 * @return 网格节点最大索引数组
 */
std::array<int,3> Block::getIMax()
{
   return iMax;
}


/**
 * @brief 获取网格维度
 * @return 网格维度 (1/2/3D)
 */
int Block::getDim()
{
   return dim;
}

/**
 * @brief 获取单元中心坐标
 * @param idim 坐标维度索引
 * @return 单元中心坐标向量
 */
std::vector<real> Block::getCellCoor(int idim)
{
   assert(idim<=dim);
   std::vector<real> res;
   res.reserve(nCel);
   for (int i = 0; i < nCel; i++)
   {
      res.push_back(coorCel(i,idim));
   }
   return res;
   
}

/**
 * @brief 获取顶点坐标
 * @param idim 坐标维度索引
 * @return 顶点坐标向量
 */
std::vector<real> Block::getVertexCoor(int idim)
{
   std::vector<real> res;
   res.reserve(nVer);
   int nVerDim=nVer;
   nVerDim/=( dim == 3 ? 1:( dim==2 ? 2 : 4));
   for (int i = 0; i < nVerDim; i++)
   {
      res.push_back(coorVer(i,idim));
   }
   return res;
   
}

/**
 * @brief 获取单元间距
 * @param i 单元索引
 * @return 单元间距向量
 */
std::vector<real> Block::getCellInterval(int i)
{
   std::vector<real> res;
   res.reserve(nVer);
   for (int idim = 0; idim < dim; idim++)
   {
      res.push_back(intervalCel(i,idim));
   }
   return res;
}

/**
 * @brief 获取指定方向上的最小网格间距
 * @param i 网格节点索引
 * @return 最小网格间距
 */
real Block::getMinDh(int i)
{
   real res=coorVer(i+1,0)-coorVer(i,0);
   if(dim>=2) res=std::min(coorVer(i+iMax[0],1)-coorVer(i,1),res);
   if(dim>=3) res=std::min(coorVer(i+iMax[0]*iMax[1],2)-coorVer(i,2),res);
   return res;
}