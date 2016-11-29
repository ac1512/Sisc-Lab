#include "PROTOTYPE.h"
// #include "functions.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
using namespace std;
const double pi = 3.14159;
void parameter::set_param(void)
{
  // double pi=3.14159;
  // cout<<"set"<<pi<<endl;
	ne=50;
  // cout<<ne<<endl;
  exl=6;
  draw_frames=1;
  ord=1;
  if(ord==1)
  {
    mth_lim=1;
    mth_rk=1;
  }
  else
  {
    mth_lim=3; //slop_limiter=1(zero),2(mm),3(mc),4(vanleer),5(superbee),6(none)
    mth_rk=2;
  }
  Flux_method=2;
  wb=0;
  bdd=1;
  if (exl ==1)
      {
        tstop     = 0.1;
      }
  else if(exl == 3)
     {
      tstop = 700;
     }
  else if(exl == 4)
      {
        tstop = 2.0;
        alpha = 0;
      }
  else if(exl == 41)
  {
   tstop = 300;
// param.alpha = pi/60;
   alpha = -pi/60;
  //  param.alpha = 0;
  }
  else if(exl == 42)
   {
   tstop = 30;
   alpha = pi/60;
   }
// % param.alpha = -pi/60;
  // % param.alpha = 0;
else if(exl == 5)
   {
      tstop= 1.0;
    }
else if(exl == 6)
   {
     tstop= 2.0;
   }
   IT = 0;
   tme = 0.0;
   eps  = 1.0e-6;
   htol = 1.0e-7;
   hdry = 1.0e-12;
   nbc = 2;
   cfl = 0.5;
   g = 9.812;
}

//function to get paramters for external use
void parameter::get_param(vector<int>& a,vector<double>& b)
{
  a[0]=ne;
  a[1]=nx;
  a[2]=exl;
  a[3]=draw_frames;
  a[4]=ord;
  a[5]=mth_lim;
  a[6]=mth_rk;
  a[7]=Flux_method;
  a[8]=wb;
  a[9]=bdd;
  a[10]=alpha;
  a[11]=IT;
  a[12]=nbc;
  b[0]=tstop;
  b[1]=eps;
  b[2]=htol;
  b[3]=hdry;
  b[4]=cfl;
  b[5]=g;
  b[6]=tme;
}

void parameter::set_paramNx(int total){
  nx = total;
}

//Prototype function for updating values in other parameter class object.
// void parameter::interfere(parameter &p, int val)
// {
//   p.nx = val;
// }
