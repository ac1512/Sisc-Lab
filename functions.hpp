#ifndef  FUNCTIONS_H_INCLUDED
#define  FUNCTIONS_H_INCLUDED
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

using namespace std;
const double pi = 3.14159;
template <typename T,typename U>
void DMS(T *xL,T *xR,U &a1)  //DMS=Domain Setup
// void DMS(double *xL,double *xR)
{
  // cout<<a1<<"    "<<a1<<endl;
  if(a1==0||a1==1)
   {
        *xL=0.0;
        *xR=1.0;
   }
   else if(a1==3)
   {
     *xL=0.0;
     *xR=120.0;
   }
   else if(a1==4)
    {
      *xL=-15.0;
      *xR=15.0;
    }
   else if(a1==41)
   {
     *xL=-15.0;
     *xR=35.0;
   }
   else if(a1==42)
   {
     *xL=-35.0;
     *xR=15.0;
   }
   else if(a1==5)
   {
     *xL=-1.0;
     *xR=1.0;
   }
   else if(a1==6)
   {
     *xL=0.0;
     *xR=1.0;
   }
};
// //
template <typename T,typename U>
void setbed(T &x,U &A,T &B,T &b) //T is double vector and U is int vector.
{
  double nx=A[2]; //value of expl.
  vector<double> xc;
  xc.resize(nx+2*B[4],0.0); //B[4] contains nbc;
  b.resize(nx+2*B[4],0.0);  //same as size of xc and contains the base profile and its the output of the function..
  for(int i=0;i<(nx+2*B[4]);i++)
  {
    xc[i]=(x[i]+x[i+1])/2;
  }
  if(A[2]==0)
  {
    double bL,bR,xgate;
    bL=-0.1;
    bR=-0.45;
    xgate=0.5;
    for(int i=0;i<(nx+2*B[4]);i++)
    {
      b[i]=bL*(xc[i]<xgate)+bR*(xc[i]>xgate);
    }
  }
  else if(A[2]==1)
  {
    for(int i=0;i<(nx+2*B[4]);i++)
    {
      b[i]=sin(pi*xc[i])*sin(pi*xc[i]);
    }
  }
  else if(A[2]==3)
  {
    double D,delta,x_a;
    D=1;
    delta=0.019;
    x_a=sqrt(4*D/(3*delta))*acosh(sqrt(20.0));
    for(int i=0;i<(nx+2*B[4]);i++)
    {
      if(xc[i]>=2*x_a)
      {
        b[i]=(xc[i]-2*x_a)/19.85;
      }
    }
  }
  else if(A[2]==4||A[2]==41||A[3]==42)
  {
    double alpha=A[10];
    for(int i=0;i<(nx+2*B[4]);i++)
    {
      b[i]=-tan(alpha)*xc[i];
    }
  }
  else if(A[2]==5)
  {
    for(int i=0;i<(nx+2*B[4]);i++)
    {
      b[i]=sin(pi*xc[i])*sin(pi*xc[i]);
    }
  }
  else if(A[2]==6)
  {
    for(int i=0;i<(nx+2*B[4]);i++)
    {
      b[i]=sin(4*pi*xc[i])*(xc[i]<=0.5)+(sin(4*pi*xc[i])-2.0)*(xc[i]>0.5);
    }
  }
}
//yogi
template <typename T>
void SetInitU(parameter &p, T &x, T &b,T &h, T &u){

	int nx=p.nx;
	int j;
	T xc[nx];
	for(j=0; j<nx;j++)
		xc[j]=0.5*(x[j]+x[j+1]);

	if (p.exl==0){
		for(j=0;j<nx;j++){
			h[j]=0.1+0.0*xc[j];
			u[j]=1.5+0.0*xc[j];
		}
	}else if(p.exl==1){
		double th;
		for(j=0;j<nx;j++){
			th=cos(2*pi*xc[j]);
			h[j]=5.0+exp(th);
			u[j]=sin(th)/h[j];
		}
	}else if(p.exl==3){
		if(p.wb==1){		//test well-balancing
			for(j=0;j<nx;j++){
				if (1 - b[j] > 0)
					h[j] = 1 - b[j];
				else
					h[j] = 0;
				u[j]=0.0*h[j];
			}
		}else{
			int D=1;
			double delta=0.019, gamma;
			gamma=sqrt(3*delta/(4*D));
			x_a=2.178272/gamma;
			for(j=0;j<nx;j++){
				if (D + delta / (cosh(gamma*(xc[j] - x_a))*cosh(gamma*(xc[j] - x_a))) - b[j] > 0)
					h[j] = D + delta / (cosh(gamma*(xc[j] - x_a))*cosh(gamma*(xc[j] - x_a))) - b[j];
				else
					h[j] = 0;
				u[j]=sqrt(p.g/D)*delta/ (cosh(gamma*(xc[j] - x_a))*cosh(gamma*(xc[j] - x_a)));
			}
		}
	}else if(p.exl==4||p.exl==41){
		double alpha=p.alpha;
		for(j=0;j<nx;j++){
			h[j]=0;
			if(xc[j]<0)
				h[j]=1+xc[j]*tan(alpha);
			u[j]=0.0;

		}
	}else if(p.exl==42||p.exl==41){
		double alpha=p.alpha;
		for(j=0;j<nx;j++){
			h[j]=0;
			if(xc[j]<0)
				h[j]=1+xc[j]*tan(alpha);
			if(h[j]<0)
				h[j]=1+xc[j]*tan(alpha);
			u[j]=0.0;
		}
	}else if (p.exl==5){
		for(j=0;j<nx;j++){
			h[j]=2+0.0*xc[j]-b[j];
			u[j]=0.0*xc[j];
		}
	}else if (p.exl==6){
		for(j=0;j<nx;j++){
			if(xc>0.5){
				if (-1.5 - b[j] > 0)
					h[j] = -1.5 - b[j];
				else
					h[j] = 0;
			}
			else{
				if (-0.5 - b[j] > 0)
					h[j] = -0.5 - b[j];
				else
					h[j] = 0;
			}
			u[j]=0.0*xc[j];
		}
	}
}
//
//
// template <typename T,typename U>
// void iniCN(T xl,T xr,T *out,U &A,U &B,parameter &p)
// {
//   //out contains the output paramters
//   //inp contains the parameters in this case its A.
//   double nx=A[0];
//   double nbc=B[4];
//   //add get_mesh;
//
//   //
//   vector<double> i;
//   i.resize(nx+1+2*nbc);
//   for(int k=0;k<(nx+1+2*nbc);k++)
//   {
//     i[k]=k+1;
//   }
//   dx=(xr-xl)/nx;
//   vector<double> x_bed;
//   for(int k=0;k<(nx+1+2*nbc);k++)
//   {
//     x_bed[i]=xl+(i[k]-(nbx+1))*dx; //with bdd cells.
//   }
//   //adding setbed
//
//
//   //
//
//   //add setInitU
//
//   //
// }
//
//













#endif
