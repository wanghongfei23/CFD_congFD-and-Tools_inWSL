/**
 * @file eigenSystem.cpp
 * @brief 特征系统类的实现文件
 */

#include"eigenSystem.hpp"

/**
 * @brief 构造函数，使用原始变量和法向量初始化
 * @param prim 原始变量数组 [rho, u, v, p]
 * @param norm_ 法向量数组
 */
eigensystemEuler2D::eigensystemEuler2D(const std::array<real,4> &prim,const std::array<real,3> & norm_)
{
    r=prim[0],u=prim[1],v=prim[2],p=prim[3];
    gamma=GAMMA,ek=(u*u+v*v)/2;
    h=p/r*gamma/(1-gamma);
    c=std::sqrt(gamma*p/r);
    norm=norm_;
    Vn=norm[0]*u+norm[1]*v;
}

/**
 * @brief 构造函数，使用左右两侧原始变量和法向量初始化
 * @param priml 左侧原始变量数组 [rho, u, v, p]
 * @param primr 右侧原始变量数组 [rho, u, v, p]
 * @param norm_ 法向量数组
 */
eigensystemEuler2D::eigensystemEuler2D(const std::array<real,4> &priml,const std::array<real,4> &primr,const std::array<real,3> & norm_)
{
      norm=norm_;
      
      gamma=GAMMA;
      enum{LL,RR};
      real rl=priml[0],ul=priml[1],vl=priml[2],pl=priml[3];
      real rr=primr[0],ur=primr[1],vr=primr[2],pr=primr[3];

      std::array<real,2> H;
      H[LL]=(ul*ul+vl*vl)/2+pl/rl*gamma/(gamma-1);
      H[RR]=(ur*ur+vr*vr)/2+pr/rr*gamma/(gamma-1);
      real coef1=std::sqrt(rl);
      real coef2=std::sqrt(rr);
      real divisor=1.0/(coef1+coef2);

      r=std::sqrt(rl*rr);
      u=(coef1*ul+coef2*ur)*divisor;
      v=(coef1*vl+coef2*vr)*divisor;
      Vn=norm[0]*u+norm[1]*v;
      real ht=(coef1*H[LL]+coef2*H[RR])*divisor;
      ek=(u*u+v*v)/2;
      h=ht-ek;
      c=std::sqrt((gamma-1)*h);
      p=r*h*((gamma-1)/gamma);

      leftEig={(-norm[Y]*u+norm[X]*v),norm[Y]              ,(-norm[X])             ,0          ,
               (h-ek)                ,u                    , v                     ,-1.0       , 
               (Vn/c+ek/h)/2         ,(-norm[X]/c-u/h)/2   ,(-norm[Y]/c-v/h)/2     ,1.0/(2.0*h),
               (-Vn/c+ek/h)/2        ,(norm[X]/c-u/h)/2    ,(norm[Y]/c-v/h)/2      ,1.0/(2.0*h)};
      
      rightEig={0                  ,1.0/h ,1           ,1          ,
                norm[Y]            ,u/h   ,u-norm[X]*c ,u+norm[X]*c,
                -norm[X]           ,v/h   ,v-norm[Y]*c ,v+norm[Y]*c,
                norm[Y]*u-norm[X]*v,ek/h  ,h+ek-Vn*c   ,h+ek+Vn*c  };
}


/**
 * @brief 将原始变量转换为特征变量
 * @param prim 原始变量数组 [rho, u, v, p]
 * @return 特征变量数组
 */
std::array<real,4> eigensystemEuler2D::primToChar(const std::array<real,4> &prim)
{
    real rt=prim[0],ut=prim[1],vt=prim[2],pt=prim[3];
    real ekt=(ut*ut+vt*vt)/2;
    real rut=rt*ut,rvt=rt*vt,ret=pt/(gamma-1)+rt*ekt;

    std::array<real,4> res;

    res[0]=rut* leftEig[1]
          +rvt* leftEig[2]
          +ret* leftEig[3]
          +rt * leftEig[0];

    res[1]=rut* leftEig[5]
          +rvt* leftEig[6]
          +ret* leftEig[7]
          +rt * leftEig[4];

    res[2]=rut* leftEig[9]
          +rvt* leftEig[10]
          +ret* leftEig[11]
          +rt * leftEig[8];

    res[3]=rut* leftEig[13]
          +rvt* leftEig[14]
          +ret* leftEig[15]
          +rt * leftEig[12];
    
    return res;

}

/**
 * @brief 将特征变量转换为原始变量
 * @param chars 特征变量数组
 * @return 原始变量数组 [rho, u, v, p]
 */
std::array<real,4> eigensystemEuler2D::charToPrim(const std::array<real,4> & chars)
{
    real ch1=chars[0],ch2=chars[1],ch3=chars[2],ch4=chars[3],rt,rut,rvt,ret;
    rt    =ch3* rightEig[2]
          +ch4* rightEig[3]
          +ch2* rightEig[1]
          +ch1* rightEig[0];

    rut   =ch3* rightEig[6]
          +ch4* rightEig[7]
          +ch2* rightEig[5]
          +ch1* rightEig[4];

    rvt   =ch3* rightEig[10]
          +ch4* rightEig[11]
          +ch2* rightEig[9]
          +ch1* rightEig[8];

    ret   =ch3* rightEig[14]
          +ch4* rightEig[15]
          +ch2* rightEig[13]
          +ch1* rightEig[12];

    real ut=rut/rt;
    real vt=rvt/rt;
    real rekt=(rut*rut+rvt*rvt)/rt/2;
    real pt=(gamma-1)*(ret-rekt);
    return{rt,ut,vt,pt};
}

eigensystemEuler1D::eigensystemEuler1D(const std::array<real,3> &priml,const std::array<real,3> &primr)
{
      gamma=GAMMA;
      enum{LL,RR};
      real rl=priml[0],ul=priml[1],pl=priml[2];
      real rr=primr[0],ur=primr[1],pr=primr[2];

      std::array<real,2> H;
      H[LL]=(ul*ul)/2+pl/rl*gamma/(gamma-1);
      H[RR]=(ur*ur)/2+pr/rr*gamma/(gamma-1);
      real coef1=std::sqrt(rl);
      real coef2=std::sqrt(rr);
      real divisor=1.0/(std::sqrt(rl)+std::sqrt(rr));

      r=std::sqrt(rl*rr);
      u=(coef1*ul+coef2*ur)*divisor;
      real ht=(coef1*H[LL]+coef2*H[RR])*divisor;
      ek=(u*u)/2;
      h=ht-ek;
      c=std::sqrt((gamma-1)*h);
      p=r*h*((gamma-1)/gamma);

      real gamma_1=gamma-1;
      leftEig={ht+c*(u-c)/gamma_1      , -u-c/gamma_1   , 1.0,
               -2.0*ht+4.0*c*c/gamma_1 , 2.0*u          ,-2.0,
               ht-c*(u+c)/gamma_1      , -u+c/gamma_1   , 1.0};
      real factorEig=0.5*gamma_1/(c*c);


      for(int ii=0;ii<9;ii++) leftEig[ii]*=factorEig;

      rightEig={1.0  , 1.0     , 1.0,
               u-c   , u       , u+c,
               ht-u*c , 0.5*u*u , ht+u*c};
}


std::array<real,3> eigensystemEuler1D::primToChar(const std::array<real,3> &prim)
{
      real rt=prim[0],ut=prim[1],pt=prim[2];
      real ekt=(ut*ut)/2;
      real rut=rt*ut,ret=pt/(gamma-1)+rt*ekt;

      std::array<real,3> res;

      res[0]=rt * leftEig[0]
            +rut* leftEig[1]
            +ret* leftEig[2];

      res[1]=rt * leftEig[3]
            +rut* leftEig[4]
            +ret* leftEig[5];

      res[2]=rt * leftEig[6]
            +rut* leftEig[7]
            +ret* leftEig[8];

    return res;
}

std::array<real,3> eigensystemEuler1D::charToPrim(const std::array<real,3> & chars)
{
    real ch1=chars[0],ch2=chars[1],ch3=chars[2],rt,rut,ret;

    rt    =ch1* rightEig[0]
          +ch2* rightEig[1]
          +ch3* rightEig[2];

    rut   =ch1* rightEig[3]
          +ch2* rightEig[4]
          +ch3* rightEig[5];

    ret   =ch1* rightEig[6]
          +ch2* rightEig[7]
          +ch3* rightEig[8];

    real ut=rut/rt;
    real Et=ret/rt;
    real ekt=(ut*ut)/2;
    real pt=(gamma-1)*(ret-rt*ekt);
    return{rt,ut,pt};
}

/* ================= 三维（任务 4 新增；与 2D 类逐行对应） ================= */

/**
 * @brief 单态构造函数（三维）：直接委托给两态构造（左右取同一状态 → 无 Roe 平均）
 */
eigensystemEuler3D::eigensystemEuler3D(const std::array<real,5> &prim,const std::array<real,3> & norm_)
{
    *this = eigensystemEuler3D(prim, prim, norm_);
}

/**
 * @brief 两态（Roe 平均）构造函数（三维）
 * @param priml 左侧原始变量 [rho, u, v, w, p]
 * @param primr 右侧原始变量 [rho, u, v, w, p]
 * @param norm_ 单位法向量（本码仅轴向）
 */
eigensystemEuler3D::eigensystemEuler3D(const std::array<real,5> &priml,const std::array<real,5> &primr,const std::array<real,3> & norm_)
{
    norm=norm_;
    gamma=GAMMA;
    enum{LL,RR};
    real rl=priml[0],ul=priml[1],vl=priml[2],wl=priml[3],pl=priml[4];
    real rr=primr[0],ur=primr[1],vr=primr[2],wr=primr[3],pr=primr[4];

    // Roe 平均（sqrt(rho) 加权；H 为总焓）
    std::array<real,2> H;
    H[LL]=(ul*ul+vl*vl+wl*wl)/2+pl/rl*gamma/(gamma-1);
    H[RR]=(ur*ur+vr*vr+wr*wr)/2+pr/rr*gamma/(gamma-1);
    real coef1=std::sqrt(rl);
    real coef2=std::sqrt(rr);
    real divisor=1.0/(coef1+coef2);

    r=std::sqrt(rl*rr);
    u=(coef1*ul+coef2*ur)*divisor;
    v=(coef1*vl+coef2*vr)*divisor;
    w=(coef1*wl+coef2*wr)*divisor;
    Vn=norm[0]*u+norm[1]*v+norm[2]*w;
    real ht=(coef1*H[LL]+coef2*H[RR])*divisor;
    ek=(u*u+v*v+w*w)/2;
    h=ht-ek;                          // 比焓 h = H - 动能
    c=std::sqrt((gamma-1)*h);         // 声速
    p=r*h*((gamma-1)/gamma);

    // 轴向法向下的切向基 t1、t2（与法向正交；与 2D 的 t=(ny,-nx) 构造同型）
    std::array<real,3> t1,t2;
    if (norm[0]>=norm[1] && norm[0]>=norm[2]) { t1={0.0,1.0,0.0}; t2={0.0,0.0,1.0}; }
    else if (norm[1]>=norm[2])                { t1={1.0,0.0,0.0}; t2={0.0,0.0,1.0}; }
    else                                      { t1={1.0,0.0,0.0}; t2={0.0,1.0,0.0}; }

    real t1q=t1[0]*u+t1[1]*v+t1[2]*w;
    real t2q=t2[0]*u+t2[1]*v+t2[2]*w;

    // 左特征矩阵（行主序，5×5）：剪切波×2 / 熵波 / 声波-
    leftEig={-t1q, t1[0], t1[1], t1[2], 0.0,
             -t2q, t2[0], t2[1], t2[2], 0.0,
             (h-ek), u, v, w, -1.0,
             (Vn/c+ek/h)/2, (-norm[0]/c-u/h)/2, (-norm[1]/c-v/h)/2, (-norm[2]/c-w/h)/2, 1.0/(2.0*h),
             (-Vn/c+ek/h)/2, (norm[0]/c-u/h)/2, (norm[1]/c-v/h)/2, (norm[2]/c-w/h)/2, 1.0/(2.0*h)};

    // 右特征矩阵（行主序，5×5）：按分量行存储（列为特征向量的对应分量）
    rightEig={0.0, 0.0, 1.0/h, 1.0, 1.0,
              t1[0], t2[0], u/h, u-norm[0]*c, u+norm[0]*c,
              t1[1], t2[1], v/h, v-norm[1]*c, v+norm[1]*c,
              t1[2], t2[2], w/h, w-norm[2]*c, w+norm[2]*c,
              t1q, t2q, ek/h, h+ek-Vn*c, h+ek+Vn*c};
}

/**
 * @brief 将原始变量转换为特征变量（三维）
 *
 * char = L·U（先把原始量组成守恒量向量，再乘左特征矩阵）。
 */
std::array<real,5> eigensystemEuler3D::primToChar(const std::array<real,5> &prim)
{
    real rt=prim[0],ut=prim[1],vt=prim[2],wt=prim[3],pt=prim[4];
    real ekt=(ut*ut+vt*vt+wt*wt)/2;
    real rut=rt*ut,rvt=rt*vt,rwt=rt*wt,ret=pt/(gamma-1)+rt*ekt;

    std::array<real,5> res;

    res[0]=rt * leftEig[0]
          +rut* leftEig[1]
          +rvt* leftEig[2]
          +rwt* leftEig[3]
          +ret* leftEig[4];

    res[1]=rt * leftEig[5]
          +rut* leftEig[6]
          +rvt* leftEig[7]
          +rwt* leftEig[8]
          +ret* leftEig[9];

    res[2]=rt * leftEig[10]
          +rut* leftEig[11]
          +rvt* leftEig[12]
          +rwt* leftEig[13]
          +ret* leftEig[14];

    res[3]=rt * leftEig[15]
          +rut* leftEig[16]
          +rvt* leftEig[17]
          +rwt* leftEig[18]
          +ret* leftEig[19];

    res[4]=rt * leftEig[20]
          +rut* leftEig[21]
          +rvt* leftEig[22]
          +rwt* leftEig[23]
          +ret* leftEig[24];

    return res;
}

/**
 * @brief 将特征变量转换为原始变量（三维）
 *
 * U = R·char（再拆回原始量）。
 */
std::array<real,5> eigensystemEuler3D::charToPrim(const std::array<real,5> & chars)
{
    real ch1=chars[0],ch2=chars[1],ch3=chars[2],ch4=chars[3],ch5=chars[4],rt,rut,rvt,rwt,ret;

    rt    =ch1* rightEig[0]
          +ch2* rightEig[1]
          +ch3* rightEig[2]
          +ch4* rightEig[3]
          +ch5* rightEig[4];

    rut   =ch1* rightEig[5]
          +ch2* rightEig[6]
          +ch3* rightEig[7]
          +ch4* rightEig[8]
          +ch5* rightEig[9];

    rvt   =ch1* rightEig[10]
          +ch2* rightEig[11]
          +ch3* rightEig[12]
          +ch4* rightEig[13]
          +ch5* rightEig[14];

    rwt   =ch1* rightEig[15]
          +ch2* rightEig[16]
          +ch3* rightEig[17]
          +ch4* rightEig[18]
          +ch5* rightEig[19];

    ret   =ch1* rightEig[20]
          +ch2* rightEig[21]
          +ch3* rightEig[22]
          +ch4* rightEig[23]
          +ch5* rightEig[24];

    real ut=rut/rt;
    real vt=rvt/rt;
    real wt=rwt/rt;
    real rekt=(rut*rut+rvt*rvt+rwt*rwt)/rt/2;
    real pt=(gamma-1)*(ret-rekt);
    return{rt,ut,vt,wt,pt};
}