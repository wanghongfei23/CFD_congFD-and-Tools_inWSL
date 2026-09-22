#include<iostream>
#include<cmath>
#include<vector>
#include<algorithm>
using namespace std;
double max(double a,double b,double c)
{
    return max(max(a,b),c);
}
int main()
{
    double Q[2][3]={1,0.5,2.5,0.125,0.1,0.25};
    double F[5][3];
    double F2[5][3];
    double a[5][3];
    int N=1;
    double dt;
    double dx;
    double gama=1.4;


    double rou_L,rou_R,u_L,u_R,H_L,H_R,p_L,p_R,c_L,c_R;
    double rou,u,H,c,p;
    double a1,a2,a3,a4,a5;
    double drou,dp,du;
    double lamda[3],lamda_L[3],lamda_R[3],abs_lamda[3];
    int i,j,k;
    vector<double> rou_0(N+2);
	vector<double>  u_0(N+2);
	vector<double> p_0(N+2);
	vector<double> m(3);
    double corrector[3][3];//第一个表示不同的变量，第二个表示不同变量的左右
    double error=1e-8;
    double eps;
    //Roe平均
    for(i=0;i<N;i++)
    {
        rou_L=Q[i][0];
        rou_R=Q[i+1][0];
        if(rou_L==0|rou_R==0)
        {
            rou_L=rou_L+error;
            rou_R=rou_R+error;
        }
        drou=rou_R-rou_L;
        rou=std::sqrt(rou_L*rou_R);
        u_L=Q[i][1]/rou_L;
        u_R=Q[i+1][1]/rou_R;
        du=u_R-u_L;
        u=(u_L*std::sqrt(rou_L)+u_R*std::sqrt(rou_R))/(std::sqrt(rou_L)+std::sqrt(rou_R));
        p_L=(gama-1)*(Q[i][2]-rou_L*u_L*u_L*0.5);
        p_R=(gama-1)*(Q[i+1][2]-rou_R*u_R*u_R*0.5);
        dp=p_R-p_L;
        c_L=std::sqrt(gama*p_L/rou_L);
        c_R=std::sqrt(gama*p_R/rou_R);
        H_L=gama*p_L/rou_L/(gama-1)+0.5*u_L*u_L;
        H_R=gama*p_R/rou_R/(gama-1)+0.5*u_R*u_R;
        H=(H_L*std::sqrt(rou_L)+H_R*std::sqrt(rou_R))/(std::sqrt(rou_L)+std::sqrt(rou_R));
        c=(gama-1)*(H-0.5*u*u);
        p=rou*c/gama;
        corrector[0][0]=u_L-std::sqrt(c_L);
        corrector[0][1]=u-std::sqrt(c);
        corrector[0][2]=u_R-std::sqrt(c_R);
        corrector[1][0]=u_L;
        corrector[1][1]=u;
        corrector[1][2]=u_R;
        corrector[2][0]=u_L+std::sqrt(c_L);
        corrector[2][1]=u+std::sqrt(c);
        corrector[2][2]=u_R+std::sqrt(c_R);
        for(j=0;j<3;j++)
        {
        	m[j]=fabs(corrector[j][1]);
		}
        /*for(j=0;j<3;j++)
        {
            eps=4.0*max(0,corrector[j][1]-corrector[j][0],corrector[j][2]-corrector[j][1]);
            if(fabs(corrector[j][1])<fabs(eps))
            {
                m[j]=fabs(pow(corrector[j][1],2)*0.5/eps+0.5*eps);
            }
        }*/
        a1=m[1]*(drou-dp/c);
        a2=m[2]*(dp+rou*std::sqrt(c)*du)/(2*c);
        a3=m[0]*(dp-rou*std::sqrt(c)*du)/(2*c);
        a4=a1+a2+a3;
        a5=std::sqrt(c)*(a2-a3);
        a[i][0]=a4;
        a[i][1]=u*a4+a5;
        a[i][2]=H*a4+u*a5-c*a1/(gama-1);
    }
    for(i=0;i<3;i++)
    {
        a[N+1][i]=a[N][i];
    }
    for(i=0;i<=N+1;i++)
    {
    	if(Q[i][0]==0)
    	{
    		Q[i][0]=Q[i][0]+error;
		}
		rou_0[i]=Q[i][0];
    	u_0[i]=Q[i][1]/Q[i][0];
    	p_0[i]=(gama-1)*(Q[i][2]-0.5*rou_0[i]*u_0[i]*u_0[i]);
	}
    for(i=0;i<N+1;i++)
    {
    F2[i][0]=rou_0[i]*u_0[i];
	F2[i][1]=rou_0[i]*u_0[i]*u_0[i]+p_0[i];
	F2[i][2]=(Q[i][2]+p_0[i])*u_0[i];
    }
     for(i=0;i<N;i++)
    {
	   F[i][0]=0.5*(F2[i][0]+F2[i+1][0])-0.5*a[i][0];
	   F[i][1]=0.5*(F2[i][1]+F2[i+1][1])-0.5*a[i][1];
	   F[i][2]=0.5*(F2[i][2]+F2[i+1][2])-0.5*a[i][2];
    }
    F[N+1][2]=F[N][2];
    F[N+1][1]=F[N][1];
    F[N+1][0]=F[N][0];
    for(i=1;i<=N;i++)
    {
	   for(j=0;j<=2;j++)
	   {
		Q[i][j]=Q[i][j]-dt*(F[i][j]-F[i-1][j])/dx ;	
	   }
    }
}

