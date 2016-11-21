#include "PROTOTYPE.h"
// #include "function.h"
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

void parameter::set_param(parameter &p)
{
  double pi=3.14159;
  cout<<"set"<<pi<<endl;
	p.ne=50;
  cout<<p.ne<<endl;
  p.exl=6;
  p.draw_frames=1;
  p.ord=1;
  if(p.ord==1)
  {
    p.mth_lim=1;
    p.mth_rk=1;
  }
  else
  {
    p.mth_lim=3; //slop_limiter=1(zero),2(mm),3(mc),4(vanleer),5(superbee),6(none)
    p.mth_rk=2;
  }
  p.Flux_method=2;
  p.wb=0;
  p.bdd=1;
  if (p.exl ==1)
      {
        p.tstop     = 0.1;
      }
  else if(p.exl == 3)
     {
      p.tstop = 700;
      }
  else if(p.exl == 4)
      {
        p.tstop = 2.0;
        p.alpha = 0;
      }
  else if(p.exl == 41)
  {
   p.tstop = 300;
// param.alpha = pi/60;
   p.alpha = -pi/60;
  //  param.alpha = 0;
  }
  else if(p.exl == 42)
   {
   p.tstop = 30;
   p.alpha = pi/60;
   }
// % param.alpha = -pi/60;
  // % param.alpha = 0;
else if(p.exl == 5)
   {
      p.tstop= 1.0;
    }
else if(p.exl == 6)
   {
     p.tstop= 2.0;
   }
   p.IT = 0;
   p.tme = 0.0;
   p.eps  = 1.0e-6;
   p.htol = 1.0e-7;
   p.hdry = 1.0e-12;
   p.nbc = 2;
   p.cfl = 0.5;
   p.g = 9.812;
};

//function to get paramters for external use
void parameter::get_param(parameter &p,vector<int>& a,vector<double>& b)
{
  a[0]=p.ne;
  a[1]=p.exl;
  a[2]=p.draw_frames;
  a[3]=p.ord;
  a[4]=p.mth_lim;
  a[5]=p.mth_rk;
  a[6]=p.Flux_method;
  a[7]=p.wb;
  a[8]=p.bdd;
  a[9]=p.alpha;
  a[10]=p.IT;
  b[0]=p.tstop;
  b[1]=p.eps;
  b[2]=p.htol;
  b[3]=p.hdry;
  b[4]=p.nbc;
  b[5]=p.cfl;
  b[6]=p.g;
}
