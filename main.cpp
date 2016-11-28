#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>
#include <math.h>
#include <vector>
#include "functions.hpp"
#include "PROTOTYPE.h"
using namespace std;
// const double pi = 3.14159;
// class parameter p;
int main()
{
  double xL=0.0,xR=0.0;
  parameter p;
  vector<int> a;
  vector<double> b,out1,out2,ho,uo;
  a.resize(11,0.0);
  b.resize(7,0.0);
  p.set_param();
  p.get_param(a,b);
//setting up the domain size
  DMS(&xL,&xR,a[2]);
  cout<<"==================================================="<<endl;
  cout<<"           make this code great again              "<<endl;
  cout<<"==================================================="<<endl;

 iniCN(xL,xR,out1,out2,ho,uo,a,b,p);





}
