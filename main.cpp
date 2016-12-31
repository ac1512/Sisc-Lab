#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>
#include <math.h>
#include <vector>
#include "functions.hpp"
#include "function2.hpp"
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
int nbc,nn,itstop;
nbc=a[12];
nn=a[0];
itstop=10;
double dt,Nx,g,cfl,itre=0.0,l2;
  l2=0.0;
ofstream norm_out;
vector<double> hresult,uresult;
norm_out.open("l2data.txt",ios::app);
while((b[6]<b[0])&&(a[11]<10000000))
{
  a[11]=a[11]+1;
  cout<<a[11]<<endl;
//out1=cn_x,out2=cn_b,ho=cn_h,uo=cn_u
Nx=a[1];
g=b[5];
cfl=b[4];
 GetDt(out1,ho,uo,&dt,g,cfl);
cout<<"dt="<<setprecision(16)<<fixed<<dt<<endl;
 Run(out1,ho,uo,out2,hresult,uresult,dt,a,b);
 b[6]=b[6]+dt;
 // cout<<"total time="<<b[6]<<endl;
 //copying hresult and u result in cn_u and cn_h
 for(int i=0;i<hresult.size();i++)
   {ho[i]=hresult[i];
       l2=l2+ho[i]*ho[i];}
 for(int i=0;i<uresult.size();i++)
    uo[i]=uresult[i];
   cout<<"l2="<<setprecision(16)<<fixed<<sqrt(l2)<<endl;
   itre=itre+1;
   cout<<"itre="<<itre<<endl;
   norm_out<<itre<<"\t"<<sqrt(l2)<<endl;
   l2=0.0;
}
norm_out.close();
}
