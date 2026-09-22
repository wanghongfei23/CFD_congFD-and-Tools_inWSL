#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>

double u0(double x,double u,double t)
{
    return 0.25+0.5*sin(M_PI*(x-u*t));
}
double du0(double x,double u,double t)
{
    return -0.5*M_PI*t*cos(M_PI*(x-u*t));
}
double f(double x,double u,double t)
{
    return u-u0(x,u,t);
}
double df(double x,double u,double t)
{
    return 1.0-du0(x,u,t);
}
int main(int argc, char** argv){
    double tend,dt;
    int n;
    std::cout<<argc<<std::endl;
    if (argc==4) 
    {
        n=std::stoi(argv[1]);
        tend=std::stod(argv[2]);
        dt=std::stod(argv[3]);
        std::cout<<"success running"<<std::endl;
    }
    else 
    {
        n=200,tend=10,dt=0.1;
    }
        std::cout<<"n=  "<<n<<"   t=  "<<tend<<"   dt=  "<<dt<<std::endl;
    double h=2.0/n,eps=1e-8;
    for(double t=0;t<tend+dt/2;t=t+dt)
    {
        std::string fname="BurgersT"+std::to_string(t)+".dat";
        std::fstream file(fname,std::ios::out);
        file<<"Title=\"exact solution of burgers equation at t = "<<t<<"\""<<'\n';
        file<< "Variables=\"x\",\"u\" \n";

        for(int i =0;i<n;i++)
        {
            double u=0.25;
            double x=h/2.0+i*h-1.0;
            double delta=100;
            do
            {
                delta=-f(u,x,t)/df(u,x,t);
                u=u+0.05*delta;
                //std::cout<<times<<' '<<theta<<'\n';
            } while (std::abs(delta)>eps);
            file<<x<<' '<<u<<'\n';
        }
        file.close();
        std::cout<<"t= "<<t<<"  finished"<<'\n';
    }

    


    
    
}
