#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>
#include <math.h>
#include <vector>
#include "functions.hpp"
#include "PROTOTYPE.h"
#include <numeric>
using namespace std;
// const double pi = 3.14159;
// class parameter p;
int main()
{
  double xL=0.0,xR=0.0;
  parameter p;
  vector<int> a;
  vector<double> b,out1,out2,ho,uo;
  a.resize(12,0.0);
  b.resize(6,0.0);
  p.set_param();
  p.get_param(a,b);
//setting up the domain size
  DMS(&xL,&xR,a[2]);
  cout<<"==================================================="<<endl;
  cout<<"           make this code great again              "<<endl;
  cout<<"==================================================="<<endl;

 iniCN(xL,xR,out1,out2,ho,uo,a,b,p);
 cout<<ho.size()<<endl;
 // for(int i=0;i<50;i++)
 // {
 //   cout<<"ho="<<ho[i]<<endl;
 // }
int nbc,nn,itstop;
nbc=a[12];
nn=a[0];
itstop=10;
// while((b[6]<b[0])&&(a[11]<10000000))
// {
//   a[11]=a[11]+1;
// }
//conversion of vectors into pointers
// cout<<out1.size()<<endl;
// cout<<out2.size()<<endl;
// cout<<ho.size()<<endl;
// cout<<uo.size()<<endl;
//
double* cn_x= new double[(int)out1.size()];
double* cn_h= new double[(int)ho.size()];
double* cn_u= new double[(int)uo.size()];
double* cn_b= new double[(int)out2.size()];
for(int i=0;i<out1.size();i++)
     cn_x[i]=out1[i];
for(int i=0;i<ho.size();i++)
     cn_h[i]=ho[i];
for(int i=0;i<uo.size();i++)
     cn_u[i]=uo[i];
for(int i=0;i<out2.size();i++)
    cn_b[i]=out2[i];
//get the timestep.
double dt,Nx=a[1],g=b[5],cfl=b[4];
GetDt(cn_x,cn_h,cn_u,&dt,Nx,g,cfl);
cout<<dt<<endl;

}
