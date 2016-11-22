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
  vector<int> A;
  vector<double> B;
  A.resize(11,0.0);
  B.resize(7,0.0);
  p.set_param();
  p.get_param(A,B);
  cout<<A[2]<<endl;
  cout<<pi<<endl;
//setting up the domain size
    DMS(&xL,&xR,A[2]);
    // vector<double> iniCN_out;  //contains the output of the iniCN

  cout<<"==================================================="<<endl;
  cout<<"||           make this code great again             ||"<<endl;
  cout<<"==================================================="<<endl;
 // double nbc;
 // nbc=b[4];






}
