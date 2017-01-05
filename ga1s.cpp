#include<vector>
#include<iostream>
using namespace std;

#include "dco.hpp"
using namespace dco;

template<class T>
void f(const vector<T>& x, T &y) {
  y=0;
  for (size_t i=0; i<x.size(); i++) y=y+x[i]*x[i];
  y=y*y;
}

void fg(const vector<double> &xv, double &yv, vector<double> &g) {
  typedef ga1s<double> DCO_M;
  typedef DCO_M::type DCO_T;
  typedef DCO_M::tape_t DCO_TAPE_T;
  size_t n=xv.size();
  vector<DCO_T> x(n); DCO_T y;
  DCO_M::global_tape=DCO_TAPE_T::create();
  for (size_t i=0;i<n;i++) {
    x[i]=xv[i];
    DCO_M::global_tape->register_variable(x[i]);
  }
  f(x,y);
  DCO_M::global_tape->register_output_variable(y);
  yv=value(y);
  derivative(y)=1;
  // DCO_M::global_tape->interpret_adjoint();
  // for (size_t i=0;i<n;i++) { g[i]=derivative(x[i]); }
  DCO_TAPE_T::remove(DCO_M::global_tape);
}

int main(int argc, char* argv[]) {
  assert(argc==2); cout.precision(15);
  size_t n=atoi(argv[1]);
  vector<double> x(n,0),g(n,0); double y=0;
  for (size_t i=0;i<n;i++) x[i]=cos(static_cast<double>(i));
  fg(x,y,g);
  cout << y << endl;
  for (size_t i=0;i<n;i++) cout << g[i] << endl;
  return 0;
}
