#ifndef  FUNCTIONS_H_INCLUDED
#define  FUNCTIONS_H_INCLUDED
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "PROTOTYPE.h"
#include <numeric>

using namespace std;
const double pi = 3.14159;
template <typename T,typename U>
void DMS(T *xL,T *xR,U &a1)  //DMS=Domain Setup
// void DMS(double *xL,double *xR)
{
  // cout<<a1<<"    "<<a1<<endl;
  if(a1==0||a1==1)
   {
        *xL=0.0;
        *xR=1.0;
   }
   else if(a1==3)
   {
     *xL=0.0;
     *xR=120.0;
   }
   else if(a1==4)
    {
      *xL=-15.0;
      *xR=15.0;
    }
   else if(a1==41)
   {
     *xL=-15.0;
     *xR=35.0;
   }
   else if(a1==42)
   {
     *xL=-35.0;
     *xR=15.0;
   }
   else if(a1==5)
   {
     *xL=-1.0;
     *xR=1.0;
   }
   else if(a1==6)
   {
     *xL=0.0;
     *xR=1.0;
   }
};
//
template <typename T,typename U>
void setbed(T &x,U &A,T &B,T &b) //T is double vector and U is int vector.
{
  double nx=A[1]; //value of expl.
  vector<double> xc;
  xc.resize(nx+2*A[12]+1,0.0); //A[12] contains nbc;
  b.resize(nx+2*A[12]+1,0.0);  //same as size of xc and contains the base profile and its the output of the function..
  for(int i=0;i<(nx+2*A[12]);i++)
  {
    xc[i]=(x[i]+x[i+1])/2;
  }
  if(A[2]==0)
  {
    double bL,bR,xgate;
    bL=-0.1;
    bR=-0.45;
    xgate=0.5;
    for(int i=0;i<(nx+2*A[12]+1);i++)
    {
      b[i]=bL*(xc[i]<xgate)+bR*(xc[i]>xgate);
    }
  }
  else if(A[2]==1)
  {
    for(int i=0;i<(nx+2*A[12]);i++)
    {
      b[i]=sin(pi*xc[i])*sin(pi*xc[i]);
    }
  }
  else if(A[2]==3)
  {
    double D,delta,x_a;
    D=1;
    delta=0.019;
    x_a=sqrt(4*D/(3*delta))*acosh(sqrt(20.0));
    for(int i=0;i<(nx+2*A[12]);i++)
    {
      if(xc[i]>=2*x_a)
      {
        b[i]=(xc[i]-2*x_a)/19.85;
      }
    }
  }
  else if(A[2]==4||A[2]==41||A[3]==42)
  {
    double alpha=A[10];
    for(int i=0;i<(nx+2*A[12]);i++)
    {
      b[i]=-tan(alpha)*xc[i];
    }
  }
  else if(A[2]==5)
  {
    for(int i=0;i<(nx+2*A[12]);i++)
    {
      b[i]=sin(pi*xc[i])*sin(pi*xc[i]);
    }
  }
  else if(A[2]==6)
  {
    for(int i=0;i<(nx+2*A[12]+1);i++)
    {
      b[i]=sin(4*pi*xc[i])*(xc[i]<=0.5)+(sin(4*pi*xc[i])-2.0)*(xc[i]>0.5);
      // cout<<"b="<<b[i]<<endl;
    }
  }
}
// //yogi
template <typename U,typename T>
void SetInitU(U &A,T &B, T &x, T &b,T &h, T &u){

	int nx=A[1];
	int j;
	T xc(nx,0.0);
	for(j=0; j<nx;j++)
		xc[j]=0.5*(x[j]+x[j+1]);

	if (A[2]==0){
		for(j=0;j<nx;j++){
			h[j]=0.1+0.0*xc[j];
			u[j]=1.5+0.0*xc[j];
		}
	}else if(A[2]==1){
		double th;
		for(j=0;j<nx;j++){
			th=cos(2*pi*xc[j]);
			h[j]=5.0+exp(th);
			u[j]=sin(th)/h[j];
		}
	}else if(A[2]==3){
		if(A[8]==1){		//test well-balancing
			for(j=0;j<nx;j++){
				if (1 - b[j] > 0)
					h[j] = 1 - b[j];
				else
					h[j] = 0;
				u[j]=0.0*h[j];
			}
		}else{
			int D=1;
			double delta=0.019, gamma,x_a;
			gamma=sqrt(3*delta/(4*D));
			x_a=2.178272/gamma;
			for(j=0;j<nx;j++){
				if (D + delta / (cosh(gamma*(xc[j] - x_a))*cosh(gamma*(xc[j] - x_a))) - b[j] > 0)
					h[j] = D + delta / (cosh(gamma*(xc[j] - x_a))*cosh(gamma*(xc[j] - x_a))) - b[j];
				else
					h[j] = 0;
				u[j]=sqrt(B[5]/D)*delta/ (cosh(gamma*(xc[j] - x_a))*cosh(gamma*(xc[j] - x_a)));
			}
		}
	}else if(A[2]==4||A[2]==41){
		double alpha=A[10];
		for(j=0;j<nx;j++){
			h[j]=0;
			if(xc[j]<0)
				h[j]=1+xc[j]*tan(alpha);
			u[j]=0.0;

		}
	}else if(A[2]==42||A[2]==41){
		double alpha=A[10];
		for(j=0;j<nx;j++){
			h[j]=0;
			if(xc[j]<0)
				h[j]=1+xc[j]*tan(alpha);
			if(h[j]<0)
				h[j]=1+xc[j]*tan(alpha);
			u[j]=0.0;
		}
	}else if (A[2]==5){
		for(j=0;j<nx;j++){
			h[j]=2+0.0*xc[j]-b[j];
			u[j]=0.0*xc[j];
		}
	}else if (A[2]==6){
		for(j=0;j<nx;j++){
			if(xc[j]>0.5){
				if (-1.5 - b[j] > 0)
					h[j] = -1.5 - b[j];
				else
					h[j] = 0;
			}
			else{
				if (-0.5 - b[j] > 0)
					h[j] = -0.5 - b[j];
				else
					h[j] = 0;
			}
			u[j]=0.0*xc[j];
		}
	}
}
//chew
//V is vector class
//T is double.
//U ins int
template <typename T, typename U, typename V>
void GetMesh(parameter &param, T *Point, U &NCell, V &x){
  int i, nitv, itv, total, istart,sz;
  total=0;
  double dx;
  // if NCell is a vector of integers that counts number of discretized grids
  // between physical distances defined in Points, for example:
  // physical coordinates in Points = [0.0, 10.0, 15.0] & NCell = [50, 25]
  // there has to be someway to count the entries in NCell, so the grid domain
  // can be properly discretized.

  nitv = NCell.size();

  //Calculates total discretized grids over the whole domain.
  total = accumulate(NCell.begin(), NCell.end(), 0.0);
  // sz=NCell.size();
  // for(int i=0;i<sz;i++)
  // {
  //   total=total+NCell[i];
  // }

  //store the sum of the cells (grid points = Cells + 1) into param struct
  param.set_paramNx(total);

  x.resize(total+1,0.0);

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
//
template <typename T,typename U,typename D>
void iniCN(D xl,D xr,T &out,T &b,T &ho,T &uo,U &A,T &B,parameter &p)
{
  //out contains the output
  //inp contains the parameters in this case its A.
  double nx=A[0];
  double nbc=A[12],dx;
  double *Point=new double[2];
  out.resize(nx+1,0.0);
  Point[0]=xl;
  Point[1]=xr;
  vector<double> Ncell;
  Ncell.resize(1,0.0);
  Ncell[0]=A[0];
  ho.resize(nx,0.0);
  uo.resize(nx,0.0);
  //add get_mesh;
GetMesh(p,Point,Ncell,out);
  p.get_param(A,B);
// for(int i=0;i<51;i++)
// {
//   cout<<"ho="<<out[i]<<endl;
// }
// //   //
  vector<double> i;
  i.resize(nx+1+2*nbc);
  for(int k=0;k<(nx+1+2*nbc);k++)
  {
    i[k]=k+1;
  }
  dx=(xr-xl)/nx;
  vector<double> x_bed;
  x_bed.resize(nx+1+2*nbc);
  for(int k=0;k<(nx+1+2*nbc);k++)
  {
    x_bed[k]=xl+(i[k]-(nbc+1))*dx; //with bdd cells.
  }
// //   // adding setbed
setbed(x_bed,A,B,b);
// for(int i=0;i<(nx+2*A[12]);i++)
// {
//   cout<<"ho="<<b[i]<<endl;
// }
// // //
vector<double> b_sec;
b_sec.resize(nx,0.0);
for(int j=0;j<nx;j++)
{
  b_sec[j]=b[nbc+j];
}
// //   //add setInitU
SetInitU(A,B,out,b_sec,ho,uo);

}
//
//

void GetSound(double g, vector<double> &h, vector<double> &c){
  int i, vec_length;

  vec_length = c.size();
  // g=9.812, b[6] in get_param function
  for (i = 0; i < vec_length; i++){
    c[i] = sqrt(g * h[i]);
  }
};

template <typename T>
void GetEigen(T &h, T &u, T &eigL, T &eigR, double g){
  int i, vec_length;

  vec_length = eigL.size();

  vector<double> c (vec_length);

  GetSound(g, h, c);

  for(i = 0; i < vec_length; i++){
    eigL[i] = u[i] - c[i];
    eigR[i] = u[i] + c[i];
  }

}

template <typename T>
void GetDt(T &x, T &h, T &u, double *dt, double g, double cfl){

  double tmp, min;
  int i, size;

  size = x.size() - 1;

  vector<double> dx (size);
  vector<double> lambda (size);
  vector<double> eigL (size);
  vector<double> tmp_calc (size);
  vector<double> eigR (size);


  GetEigen(h, u, eigL, eigR, g);

  min = 1e10;
  for (i = 0; i < size ; i++){
    dx[i] = x[i+1] - x[i];
    tmp = (fabs(eigL[i]) >= fabs(eigR[i]))? fabs(eigL[i]):fabs(eigR[i]);
    lambda[i] = tmp + 1e-16;

    tmp_calc[i] = dx[i]/lambda[i];
    min = (min < tmp_calc[i])? min:tmp_calc[i];
  }

  tmp = min*cfl;
  if (tmp < 1e-10){
    printf("Error: Dt = %e is too small!\n", tmp);
  }

  *dt = tmp;
}

void GetHydroPre(vector<double> &h, double g, vector<double> &p){
  int i, size_h;

  size_h = h.size();

  for (i = 0; i < size_h; i++){
    p[i] = 0.5*g*h[i]*h[i];
  }
}

template <typename T>
void GetFlux( T &h, T &u, T &f1, T &f2, parameter &p) {
	int n = sizeof(h);
	for (int j = 0; j < n;j++) {
		f1[j] = h[j] * u[j];
		f2[j] = h[j] * u[j] * u[j] + 0.5*p.g*h[j] * h[j];
	}
}

template <typename T>
void GetNumFlux(T &hL, T &uL, T &hR, T &uR, T &fh, T &fm, parameter &p){

	int fluxmethod=p.Flux_method;
	int j, n=sizeof(uL);
	double fhL[n], fmL[n], fhR[n], fmR[n],eigL1[n], eigL2[n], eigR1[n], eigR2[n], mL[n], mR[n], sL[n], sR[n],s[n];
	GetFlux(hL,uL, fhL, fmL, p);
	GetFlux(hR, uR, fhR, fmR, p);
	GetEigen(hL,uL,n,eigL1,eigL2,p.g);
	GetEigen(hR, uR, n, eigR1, eigR2, p.g);
	for (int j = 0; j < n; j++) {
		mL[j] = hL[j] * uL[j];
		mR[j] = hR[j] * uR[j];
		sL[j] = (((eigL1 < eigR1) ? eigL1 : eigR1) < pow(10, -8)) ? ((eigL1 < eigR1) ? eigL1 : eigR1):pow(10, -8);
		sR[j] = (((eigL2>eigR2) ? eigL2 : eigR2) > pow(10, -8)) ? ((eigL2 > eigR2) ? eigL2 : eigR2):pow(10, -8);
	}
	if (fluxmethod == 1) {
		for (int j = 0; j < n; j++) {
			s[j] = (sL[j] > sR[j]) ? sL[j] : sR[j];
			fh[j] = 0.5*(fhL + fhR - s[j] * (hR[j] - hL[j]));
			fm[j] = 0.5*(fmL + fmR - s[j] * (mR[j] -mL[j]));
		}
	}else if (fluxmethod == 2) {
		for (int j = 0; j < n; j++) {
			fh[j] = (sR[j]*fhL - sL[j]* fhR + sL[j] *sR[j]* (hR[j] - hL[j]))/(sR[j]-sL[j]);
			fm[j] = (sR[j] * fmL - sL[j] * fmR + sL[j] * sR[j] * (mR[j] - mL[j])) / (sR[j] - sL[j]);
		}
	}else if (fluxmethod == 3) {
		for (int j = 0; j < n; j++) {
			fh[j] = (sR[j] * fhL - sL[j] * fhR - sL[j] * sR[j] * (hL[j] - hR[j])) / (sR[j] - sL[j]);
			fm[j] = (sR[j] * fmL - sL[j] * fmR - sL[j] * sR[j] * (mL[j] - mR[j])) / (sR[j] - sL[j]);
		}
	}
}

template <typename T, typename U>
void SetBDC(U &h0, U &u0, U &b0,
            T &paramvec_a, U &paramvec_b, U &h, U &u, U &b){
  int i, j, b0_length;
  int nx, iL, iR, nbc, *iBL, *iBR;
  double g;

  b0_length = b0.size();

  nx = paramvec_a[1]; //parameter variable of Nx
  nbc = paramvec_a[12]; //parameter variable of nbc
  g = paramvec_b[6]; //parameter variable of g

  for(i = 0; i < b0_length ; i++){
    b[i] = b0[i];
  }

  for (i = 0; i < nx; i++){
    h[i+nbc] = h0[i];
    u[i+nbc] = u0[i];
  }

  iL = nbc;
  iR = iL + nx - 1;

  for (i = 0; i < iL; i++){
    iBL[i] = i; // indices of lhs ghost cells
  }

  for (i = 0 ; i < nbc; i++){ //TODO :Check the indices properly!
    iBR[i] = iR + i; // indices of rhs ghost cells
  }

  if (paramvec_a[2] == 0){
  // left boundary condition
    for (i = 0; i < nbc; i++){
      j = iL - i;
      b[j] = b[iL];
      h[j] = h[iL];
      u[j] = u[iL];
    }
  // right boundary condition
    for (i = 0; i < nbc; i++){
      j = iR + i;
      b[j] = b[iR];
      h[j] = h[iR];
      u[j] = u[iR];
    }
  }
  else if (paramvec_a[2] == 1){
    if (paramvec_a[9] = 1){
      h[1] = h[nx + nbc - 1];
      h[0] = h[nx + nbc - 2];
      u[1] = u[nx-1];
      u[0] = u[nx + nbc - 2];

      h[nx + nbc] = h[nbc];
      h[nx + nbc + 1] = h[nbc + 1];
      u[nx + nbc] = u[nbc];
      u[nx + nbc + 1] = u[nbc + 1];
    }
    else{
      for (i = 0; i < iL ; i++){
        iBL[i] = i; // indices of lhs ghost cells
        h[i] = h[iR];
        u[i] = u[iR];
      }

      for (i = 0 ; i < nbc; i++){ //TODO :Check the indices properly!
        iBR[i] = iR + i; // indices of rhs ghost cells
        h[iBR[i]] = h[iL];
        u[iBR[i]] = u[iL];
      }
    }
  }
  else if (paramvec_a[2] == 3){
    h[1] = 1.0;
    h[0] = 1.0;

    for (i = 0; i < iL; i++){
      u[iBL[i]] = 0.0;
    }
    for (i = 0 ; i < nbc; i++){
      h[iBR[i]] = 0.0 * iBR[i];
      u[iBR[i]] = 0.0 * iBR[i];
    }
  }
  else if (paramvec_a[2] == 4){
    for (i = 0; i < iL; i++){
      h[iBL[i]] = h[iL];
      b[iBL[i]] = b[iL];
    }

    if (h[nbc] >= paramvec_b[2]){
      for (i = 0; i < iL; i++){
        u[iBL[i]] = u[nbc]*h[nbc]/h[nbc-1];
      }
    }
    else{
      for (i = 0; i < iL; i++){
        u[iBL[i]] = 0.0;
      }
    }

    for (i = 0 ; i < nbc; i++){
      h[iBR[i]] = h[iR];
      b[iBR[i]] = b[iR];
    }

    if (h[iR+1] >= paramvec_b[2]){
      for (i = 0 ; i < nbc; i++){
        u[iBR[i]] = u[iR]*h[iR]/h[iR+1];
      }
    }
    else{
      for (i = 0 ; i < nbc; i++){
        u[iBR[i]] = 0.0;
      }
    }
  }
  else if (paramvec_a[2] == 41){
    for (i = 0; i < iL; i++){
      h[iBL[i]] = h[iL];
      u[iBL[i]] = 0.0;
    }
    for (i = 0 ; i < nbc; i++){
      h[iBR[i]] = h[iR];
    }

    if (h[iR+1] >= paramvec_b[2]){
      for (i = 0 ; i < nbc; i++){
        u[iBR[i]] = u[iR]*h[iR]/h[iR+1];
      }
    }
    else{
      for (i = 0 ; i < nbc; i++){
        u[iBR[i]] = 0.0;
      }
    }
  }
  else if (paramvec_a[2] == 42){
    for (i = 0; i < iL; i++){
      h[iBL[i]] = h[iL];
      u[iBL[i]] = 0.0;
    }
    for (i = 0 ; i < nbc; i++){
      h[iBR[i]] = h[iR];
      u[iBR[i]] = u[iR];
    }
  }
  else if (paramvec_a[2] == 5){
    for (i = 0; i < iL ; i++){
      h[iBL[i]] = h[iR];
      u[iBL[i]] = 0.0;
    }
    for (i = 0 ; i < nbc; i++){
      h[iBR[i]] = h[iL];
      u[iBR[i]] = 0.0;
    }
  }
  else if (paramvec_a[2] == 6){
    for (i = 0; i < iL ; i++){
      h[iBL[i]] = h[iL];
      u[iBL[i]] = 0.0;
    }
    for (i = 0 ; i < nbc; i++){
      h[iBR[i]] = h[iR];
      u[iBR[i]] = 0.0;
    }
  }
}

 #endif
