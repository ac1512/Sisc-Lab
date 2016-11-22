#ifndef  FUNCTIONS_H_INCLUDED
#define  FUNCTIONS_H_INCLUDED
// #include "functions.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "PROTOTYPE.h"

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

template <typename T, typename U>
void GetMesh(parameter &param, T Point, U NCell, T x){
  int i, nitv, itv, total, istart;
  double dx;
  // if NCell is a vector of integers that counts number of discretized grids
  // between physical distances defined in Points, for example:
  // physical coordinates in Points = [0.0, 10.0, 15.0] & NCell = [50, 25]
  // there has to be someway to count the entries in NCell, so the grid domain
  // can be properly discretized.

  nitv = NCell.size();

  //Calculates total discretized grids over the whole domain.
  total = accumulate(NCell.begin(), NCell.end(), 0.0);

  //store the sum of the grid points into param struct
  param.set_paramNx(total);

  // filling in physical coordinates for grid points into x vector
  istart = 0;
  for (itv = 0; itv < nitv; itv++){
    dx = (Point[itv+1] - Point[itv])/NCell[itv];
    for (i = 0; i < NCell[itv]; i++){
      x[i+istart] = Point[itv] + i*dx;
    }
    istart = i;
  }

  //coordinate of last grid point
  x[total] = Point[itv];
};

void GetSound(double g, double *h, double *c, int vec_length){
  int i;

  // g=9.812, b[6] in get_param function
  for (i = 0; i < vec_length; i++){
    c[i] = sqrt(g * h[i]);
  }
};

#endif
