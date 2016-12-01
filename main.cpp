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
 // for(int i=0;i<50;i++)
 // {
 //   cout<<"ho="<<ho[i]<<endl;
 // }
int nbc,nn,itstop;
nbc=a[12];
nn=a[0];
itstop=10;
//out1=cn_x,out2=cn_b,ho=cn_h,uo=cn_u
double dt,Nx=a[1],g=b[5],cfl=b[4];
GetDt(out1,ho,uo,&dt,g,cfl);
cout<<dt<<endl;
cout<<"uo"<<uo.size()<<endl;
vector<double> hresult,uresult;
 Run(out1,ho,uo,out2,hresult,uresult,dt,a,b);
vector<double> d(4,0.0),gg(4,0.0);
}
