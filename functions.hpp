#ifndef  FUNCTIONS_H_INCLUDED
#define  FUNCTIONS_H_INCLUDED
// #include "functions.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

using namespace std;

template <typename T,typename U>
void DMS(T *xL,T *xR,U &a)  //DMS=Domain Setup
// void DMS(double *xL,double *xR)
{
  cout<<a[0]<<"    "<<a[1]<<endl;
  if(a[1]==0||a[1]==1)
   {
        *xL=0.0;
        *xR=1.0;
   }
   else if(a[1]==3)
   {
     *xL=0.0;
     *xR=120.0;
   }
   else if(a[1]==4)
    {
      *xL=-15.0;
      *xR=15.0;
    }
   else if(a[1]==41)
   {
     *xL=-15.0;
     *xR=35.0;
   }
   else if(a[1]==42)
   {
     *xL=-35.0;
     *xR=15.0;
   }
   else if(a[1]==5)
   {
     *xL=-1.0;
     *xR=1.0;
   }
   else if(a[1]==6)
   {
     *xL=0.0;
     *xR=1.0;
   }
};
#endif
