// #include "SpaceDis.hpp"

// void SpaceDis::initBuffer()
// {
//     flux.clear();
//     var.clear();
//     fluxNode.clear();
//     center=0;

//     flux.resize(nVar);
//     var.resize(nPrim);
//     fluxNode.resize(nVar);


//     int nFluxHalf,nFluxNode,nVarNode;
//     if (diffMethod==TRAD2) {nFluxHalf=2;nFluxNode=0;}
//     else if(diffMethod==TRAD6) {nFluxHalf=6;nFluxNode=0;}
//     else if(diffMethod==HDS6) {nFluxHalf=2;nFluxNode=5;}
//     else{
//         std::cout<<"SpaceDis error: unexpected diffmethod type\n";
//     }

//     //确定节点通量需要计算的数量和
//     if(interMethod==FIRSTORDER) {nVarNode=std::max(nFluxNode,2);}
//     else if(interMethod==MUSCL){nVarNode=std::max(nFluxNode,4);}
//     else if(interMethod==WCNS5){nVarNode=std::max(nFluxNode,6);}
//     else if(interMethod==WCNSZ5){nVarNode=std::max(nFluxNode,6);}
//     else{
//         std::cout<<"SpaceDis error: unexpected intermethod type\n";
//     }
//     centerOffset=(nVarNode-1)/2;
//     assert(centerOffset> ((nVarNode-1)/2-0.1));


//     //循环数组空间初始化
//     for(auto iflux:flux) iflux.set_capacity(nFluxHalf);
//     for(auto ivar:var) ivar.set_capacity(nVarNode);
//     for(auto ifnode:fluxNode) ifnode.set_capacity(nFluxNode);

//     //将刚开始的值push进circular buffer
//     for(int i=-centerOffset;i<=centerOffset;i++)
//     {
//         auto variables=(this->*variablesGetter)(i);
//         for(int ivar=0;ivar<nPrim;ivar++)
//         {
//             var[ivar].push_back(variables[ivar]);
//         }
//         if(var[0].full())
//         {
            
//         }
//     }
    
    
// }

// std::vector<real> SpaceDis::getEulerPrim(int i)
// {
//     std::vector<real> res(nPrim);
//     for (int ivar = 0; ivar < nPrim; ivar++) res[ivar]=at(i,ivar);
//     return res;
// }
