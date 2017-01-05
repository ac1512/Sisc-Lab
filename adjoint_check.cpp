#include<vector>
#include<iostream>
using namespace std;

#include "dco.hpp"
using namespace dco;

template<class T>
void f(const vector<vector<T> >& x, T &y)
{
  y=0;
  int breadth=x[1].size();
  int length=x.size();
  for (size_t i=0; i<length; i++)
   {
     for(int j=0;j<breadth;j++)
     {
       y=y+x[i][j]*x[i][j];
     }
   }
}

template<typename T>
void fset(vector<vector<T> > &x)
{
  int breadth=x[1].size();
  int length=x.size();
  for(int i=1;i<length;i++)
  {
    for(int j=0;j<breadth;j++)
    {
      x[i][j]=0;
      x[i][j]=x[i][j]+sin(x[i-1][j])*x[i-1][j];
    }
  }
}

template<typename T>
void fg(vector<vector<T> > &x, T &y,vector<vector<T> > derv)
{
  typedef ga1s<T> DCO_M;
  typedef typename DCO_M::type DCO_T;
  typedef typename DCO_M::tape_t DCO_TAPE_T;
  int breadth=x[1].size();
  int length=x.size();
  vector<vector<DCO_T> > xdco;
  DCO_T ydco;
  xdco.resize(length,vector<DCO_T>(breadth));
  // cout<<xdco[1].size()<<endl;
  DCO_M::global_tape=DCO_TAPE_T::create();
  for(int j=0;j<breadth;j++)
  {
    xdco[0][j]=x[0][j];
   for(int i=0;i<length;i++)
   {
     DCO_M::global_tape->register_variable(xdco[i][j]);
   }
  }
  fset(xdco);
  f(xdco,ydco);
  DCO_M::global_tape->register_output_variable(ydco);
  y=value(ydco);
  derivative(ydco)=1;
  DCO_M::global_tape->interpret_adjoint();
  for(int i=0;i<length;i++)
  {
  for(int j=0;j<breadth;j++)
  {
   derv[i][j]=derivative(xdco[i][j]);
   cout<<derv[i][j]<<"   ";
  }
  cout<<endl;
}

   DCO_TAPE_T::remove(DCO_M::global_tape);
}
int main()
{
  double y;
  vector<vector<double> > x,derv;
  x.resize(10,vector<double> (10,0.0));
  derv.resize(10,vector<double> (10,0.0));
  for (int j=0;j<10;j++)
  {
    x[0][j]=j+sin(j)*(j+2);
  }
  fg(x,y,derv);
}
