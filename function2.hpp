#ifndef  FUNCTION2_H_INCLUDED
#define  FUNCTION2_H_INCLUDED
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "PROTOTYPE.h"
#include "functions.hpp"

using namespace std;

double limit(double dul, double dur, int mth){
	double a,b, du;
	int sign;
	a=fabs(dul);
	b=fabs(dur);
	if(dul<0 && dur<0)
		sign=-1;
	else if (dul>0 && dur>0)
		sign=1;
	else
		sign=0;

	if(mth==1)
		du=a*0;
	else if(mth==2)
		du=(a<b)?a:b;
	else if(mth==3){
		double th=1.3;
		du=((th*a<th*b)?th*a:th*b)<(0.5*(a+b))?((th*a<th*b)?th*a:th*b):(0.5*(a+b));
	}else if(mth==4){
		double eps=pow(10,-8);
		du=2*a*b/(((a+b)>eps)?(a+b):eps);
	}else if (mth==5){
		double m1=(a<2*b)?a:2*b;
		double m2=(2*a<b)?2*a:b;
		du=(m1>m2)?m1:m2;
	}else
		du=0.5*(a+b);
	du=sign*du;
	return du;
}

template <typename T, typename U>
void RecCH( T &x, T &h, T &u, T &b, T &wL, T &hL, T &uL, T &bL, T &wR, T &hR, T &uR, T &bR, U &a) {
	int nx = a[1], nbc = a[12];
	int j, n=sizeof(h);
	double w[n], dw[nx+nbc+1], db[nx + nbc + 1], du[nx + nbc + 1], dh[nx + nbc + 1];
	for (j = 0; j < n; j++)
		w[j] = h[j] + b[j];
	for (j = nbc - 1; j < (nx + nbc + 1); j++) {
		dw[j] = limit(w[j] - w[j - 1], w[j + 1] - w[j], a[5]);
		db[j] = limit(b[j] - b[j - 1], b[j + 1] - b[j], a[5]);
		du[j] = limit(u[j] - u[j - 1], u[j + 1] - u[j], a[5]);
		dh[j] = dw[j] - db[j];
	}
	for (j = nbc - 1; j < (nx + nbc + 1); j++) {
		if (h[j] - 0.5*dh[j] < 0) {
			dh[j] = 2*h[j];
			db[j] = dw[j] - dh[j];
		}
		hR[j] = h[j] - 0.5*dh[j];
		uR[j] = u[j] - 0.5*du[j];
		wR[j] = w[j] - 0.5*dw[j];
		bR[j] = b[j] - 0.5*db[j];
	}
	for (j = nbc - 1; j < (nx + nbc + 1); j++) {
		if (h[j] +0.5*dh[j] < 0) {
			dh[j] = -2*h[j];
			db[j] = dw[j] - dh[j];
		}
	}
	for (j = nbc; j < (nx + nbc + 1); j++) {
		hL[j] = h[j-1] + 0.5*dh[j - 1];
		uL[j] = u[j - 1] + 0.5*du[j - 1];
		wL[j] = w[j - 1] + 0.5*dw[j - 1];
		bL[j] = b[j - 1] + 0.5*db[j - 1];
	}
}

// template <typename T, typename U>
// void GetRHS(T &x,T &h, T&u, T &b, T &rhsH, T &rhsM, U &a, double g){//g=parameter b[5]
//
// 	int nx=a[1], nbc=a[12];
// 	int j, n=sizeof(u);
//
// 	vector<double> hL(n,0.0), hR(n,0.0), wL(n,0.0), wR(n,0.0), uL(n,0.0), uR(n,0.0), bL(n,0.0),bR(n,0.0);
// 	RecCH(x, h, u, b, wL, hL, uL, bL, wR, hR, uR, bR, a);
//
//
// 	vector<double> bmax(nx + nbc + 1,0.0), bo(nx + nbc + 1,0.0), wmin(nx + nbc + 1,0.0), hLs(nx + nbc + 1,0.0), hRs(nx + nbc + 1,0.0), dpL(nx + nbc + 1,0.0), dpR(nx + nbc + 1,0.0);
// 	for (j = nbc; j < (nx + nbc + 1); j++) {
// 		bmax[j] = (bL[j] > bR[j]) ? bL[j] : bR[j];
// 		wmin[j] = (wL[j] < wR[j]) ? wL[j] : wR[j];
// 		bo[j] = (bmax[j] < wmin[j]) ? bmax[j] : wmin[j];
// 		hLs[j] = ((wL[j] - bo[j]) < hL[j]) ? (wL[j] - bo[j]) : hL[j];
// 		hRs[j] = ((wR[j] - bo[j]) < hR[j]) ? (wR[j] - bo[j]) : hR[j];
// 		dpL[j] = g*(hL[j] + hLs[j])*(bo[j]-bL[j]) / 2;
// 		dpR[j] = g*(hR[j] + hRs[j])*(bo[j] - bR[j]) / 2;
// 	}
// 	vector<double> fh(nx + nbc+1,0.0), fm(nx + nbc + 1,0.0), fmL(nx + nbc + 1,0.0), fmR(nx + nbc + 1,0.0);
// 	GetNumFlux((hLs+nbc), (uL + nbc), (hRs + nbc), (uR+ nbc), (fh + nbc), (fm + nbc), a, g);
// 	for (j = nbc; j < (nx + nbc+1); j++) {
// 		fmL[j] = fm[j] + dpL[j];
// 		fmR[j] = fm[j] + dpR[j];
// 	}
// 	double dx[nx];
// 	for (j = 0; j < nx; j++)
// 		dx[j] = x[j + 1] - x[j];
// 	for (j = 0; j < nx-1; j++) {
// 		rhsH[j] = -1 * (fh[j + nbc + 1] - fh[nbc]) / dx[j];
// 		rhsM[j]= -1 * (fmL[j + nbc + 1] - fmR[nbc]) / dx[j]- g*h[j+nbc]*(bL[j + nbc + 1] - bR[nbc])/dx[j];
// 	}
// }

template<typename T,typename U>
void Run(T &cnx,T &cnh,T &cnu,T &cnb,T &hresult,T &uresult,double &dt,U &A,T &B)
{
  int nx,ord;
  double eps;
  nx=A[0];
  eps=B[1];
  ord=A[6];
  double* hu=new double[(int) cnh.size()]();
  double* h1=new double[(int) cnh.size()]();
  double* hu1=new double[(int) cnh.size()]();
  vector<double> h(cnb.size(),0.0);
  vector<double> u(cnb.size(),0.0);
  vector<double> b(cnb.size(),0.0);
  vector<double> rhsH(nx,0.0);
  vector<double> rhsM(nx,0.0);
  for(int i=0;i<cnh.size();i++)
     hu[i]=cnh[i]*cnu[i];
  int mrk=ord;
  // for(int i=0;i<mrk;i++)
  // {
  //    if((mrk==2)&&(i==1))
  //    {
  //      for(int k=0;k<cnh.size();k++)
  //      {
  //        h1[k]=hu[k];
  //        hu1[k]=hu[k];
  //      }
  //    }
  //    SetBDC(cnh,cnu,cnb,A,B,h,u,b);
  //    GetRHS(cnx,h,u,b,rhsH,rhsM,A,B[5]);
  //    for(int k=0;k<cnh.size();k++)
  //     {
  //       cnh[k]=cnh[k]+dt*rhsH[k];
  //       hu[k]=hu[k]+dt*rhsM[k];
  //     }
  //    if((mrk=2)&&(i==2))
  //    {
  //      for(int k=0;k<cnh.size();k++)
  //      {
  //        cnh[k]=0.5*(cnh[k]+h[k]);
  //        hu[k]=0.5*(hu1[k]+hu[k]);
  //      }
  //    }
  //    for(int k=0;k<cnu.size();k++)
  //    {
  //      cnu[k]=hu[k]/max(h[k],b[1]);
  //    }
	//
      // }
  // hresult.resize(cnh.size(),0.0);
  // uresult.resize(cnh.size(),0.0);
  // for(int k=0;k<cnh.size();k++)
  //     hresult[k]=cnh[k];
  // for(int k=0;k<cnu.size();k++)
  //     uresult[k]=cnu[k];
}
 #endif
