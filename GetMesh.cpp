#include <iostream>
#include <string>

using namespace std;

void GetMesh(parameter &param, double &Point, int &NCell, double &x){
  int i, nitv, itv, total, istart;
  double dx;
  // if NCell is a vector of integers that counts number of discretized grids
  // between physical distances defined in Points, for example:
  // physical coordinates in Points = [0.0, 10.0, 15.0] & NCell = [50, 25]
  // there has to be someway to count the entries in NCell, so the grid domain
  // can be properly discretized.

  nitv = 0; //TODO: placeholder of NCell, tracks length of NCell

  //Calculates total discretized grids over the whole domain.
  total = 0;
  for (i = 0; i < nitv; i++){
    total += NCell[i]; //sum(NCell) in matlab
  }

  //store the sum of the grid points into param struct
  param.setNx(total); //TODO: may require you to add into parameter class
                      //TODO: length of x vector depends on this variable

  // filling in physical coordinates for grid points into x vector
  istart = 0;
  for (itv = 0; itv < nitv; itv++){
    dx = (Point[itv+1] - Point[itv])/NCell[itv];
    for (i = istart; i < NCell[itv]; i++){
      x[i] = Point[itv] + (i-istart)*dx;
    }
    istart = i;
  }

  //coordinate of last grid point
  x[total-1] = Point[itv-1];


  //TODO: need to test run to see if the coordinates at the bounds are properly
  //      accounted for.

}
