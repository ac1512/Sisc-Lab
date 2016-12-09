#ifndef  FUNCTION2_H_INCLUDED
#define  FUNCTION2_H_INCLUDED
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "PROTOTYPE.h"
#include "functions.hpp"
#include <iomanip>
using namespace std;


double limit(double dul, double dur, int mth){
	double a,b, du;
	int sign;
	a=fabs(dul);
	b=fabs(dur);
	// cout<<"a="<<a<<"b="<<b<<endl;
	if(dul<0 && dur<0)
		sign=-1.0;
	else if (dul>0 && dur>0)
		sign=1.0;
	else
		sign=0;

	if(mth==1)
		du=a*0.0;
	else if(mth==2)
		du=minimum(a,b);
	else if(mth==3){
		double th=1.3;
		du=maximum(minimum(th*a,th*b),0.5*(a+b));
	}else if(mth==4){
		double eps=1.0e-8;
		du=2.0*a*b/(maximum((a+b),eps));
	}else if (mth==5){
		double m1=minimum(a,2.0*b);
		double m2=minimum(2.0*a,b);
		du=maximum(m1,m2);
	}else
		du=0.5*(a+b);
	du=sign*du;
	return du;
}

template <typename T, typename U>
void RecCH( T &x, T &h, T &u, T &b, T &wL, T &hL, T &uL, T &bL, T &wR, T &hR, T &uR, T &bR, U &a) {
	int nx = a[1], nbc = a[12];
	int j, n=h.size();
	double w[n], dw[n], db[n], du[n], dh[n];
	for (j = 0; j < n; j++)
		{w[j] = h[j] + b[j];
			// cout<<"h["<<j<<"]="<<h[j]<<"  "<<b[j]<<endl;
		}
	for (j = nbc - 1; j < (nx + nbc + 1); j++) {
		dw[j] = limit(w[j] - w[j - 1], w[j + 1] - w[j], a[5]);
		db[j] = limit(b[j] - b[j - 1], b[j + 1] - b[j], a[5]);
		du[j] = limit(u[j] - u[j - 1], u[j + 1] - u[j], a[5]);
		dh[j] = dw[j] - db[j];
		 //cout<<"dh="<<dh[j]<<"db="<<db[j]<<"dw "<<dw[j]<<" du"<<du[j]<<endl;
	}
	for (j = nbc - 1; j < (nx + nbc+1); j++) {
		// cout<<"hR="<<hR[j]<<endl;
		if (h[j] - 0.5*dh[j] < 0) {
			// cout<<"j="<<j+nbc-1<<endl;
			dh[j+nbc-1] = 2.0*h[j+nbc-1];
			db[j+nbc-1] = dw[j+nbc-1] - dh[j+nbc-1];
		}
	}
	for (j = nbc - 1; j < (nx + nbc+1); j++) {
		hR[j] = h[j] - 0.5*dh[j];
		uR[j] = u[j] - 0.5*du[j];
		wR[j] = w[j] - 0.5*dw[j];
		bR[j] = b[j] - 0.5*db[j];
		// cout<<"wR["<<j<<"]="<<wR[j]<<endl;
	}
	for (j = nbc - 1; j < (nx + nbc+1); j++) {
		if (h[j] +0.5*dh[j] < 0) {
			dh[j+nbc-1] = -2.0*h[j+nbc-1];
			db[j+nbc-1] = dw[j+nbc-1] - dh[j+nbc-1];
		}
	}
	for (j = nbc; j < (nx + nbc + 1); j++) {
		hL[j] = h[j-1] + 0.5*dh[j - 1];
		uL[j] = u[j - 1] + 0.5*du[j - 1];
		wL[j] = w[j - 1] + 0.5*dw[j - 1];
    bL[j] = b[j - 1] + 0.5*db[j - 1];
		// cout<<"w["<<j<<"]="<<w[j-1]<<endl;
	}
}

 template <typename T, typename U>
 void GetRHS(T &x,T &h, T&u, T &b, T &rhsH, T &rhsM, U &a, double g){//g=parameter b[5]

 	int nx=a[1], nbc=a[12];
 	int j, n=h.size();
 	vector<double> hL(n,0.0), hR(n,0.0), wL(n,0.0), wR(n,0.0), uL(n,0.0), uR(n,0.0), bL(n,0.0),bR(n,0.0);
 	RecCH(x, h, u, b, wL, hL, uL, bL, wR, hR, uR, bR, a);


 	vector<double> bmax(nx + nbc + 1,0.0), bo(nx + nbc + 1,0.0), wmin(nx + nbc + 1,0.0), hLs(nx + nbc + 1,0.0), hRs(nx + nbc + 1,0.0), dpL(nx + nbc + 1,0.0), dpR(nx + nbc + 1,0.0);

 	for (j = nbc; j < (nx + nbc + 1); j++) {
 		bmax[j] = maximum(bL[j],bR[j]);
		wmin[j] = minimum(wL[j],wR[j]);
 		bo[j] = minimum(bmax[j],wmin[j]);
 		hLs[j] = minimum((wL[j] - bo[j]),hL[j]);
		hRs[j] = minimum((wR[j] - bo[j]),hR[j]);
 		dpL[j] = g*(hL[j] + hLs[j])*(bo[j]-bL[j]) / 2.0;
 		dpR[j] = g*(hR[j] + hRs[j])*(bo[j] - bR[j]) / 2.0;
		//cout<<"dpR[]"<<j<<"="<<dpR[j]<<setprecision(20)<<fixed<<"dpL[]"<<j<<"="<<dpL[j]<<endl;
 	}
	vector<double> fh(nx+1,0.0), fm(nx+1,0.0), fmL(nx+1,0.0), fmR(nx+1,0.0);
	vector<double> hLss(nx+1,0.0),uLs(nx+1,0.0),hRss(nx+1,0.0),uRs(nx+1,0.0);
	for(int j=0;j<nx+1;j++){
  		hLss[j]=hLs[nbc+j];
		hRss[j]=hRs[nbc+j];
		uLs[j]=uL[nbc+j];
		uRs[j]=uR[nbc+j];
		// cout<<"hls=["<<j<<"]="<<hLss[j]<<"    fmr=["<<j<<"]="<<hRss[j]<<endl;
	}
	GetNumFlux(hLss, uLs, hRss, uRs, fh, fm, a, g);
	for (j = 0; j < (nx+1); j++) {
		fmL[j] = fm[j] + dpL[j+nbc];
		fmR[j] = fm[j] + dpR[j+nbc];
	}
	vector<double> dx(nx,0.0);
	for (j = 0; j < nx; j++)
		{dx[j] = x[j + 1] - x[j];
			// cout<<"dx=["<<j<<"]="<<dx[j]<<endl;
		}
	for (j = 0; j < nx; j++) {
		rhsH[j] = -1.0* (fh[j+ 1] - fh[j]) / dx[j];
		rhsM[j]= -1.0 * (fmL[j+ 1] - fmR[j]) / dx[j]- g*h[j+nbc]*(bL[j + nbc + 1] - bR[j+nbc])/dx[j];
		// if(j==20)
		//  printf(" rhsH[%d]=%0.16e  rhsM[%d]=%0.16e  fml[20]=%0.16e fmR[20]=%0.16e dx[20]=%0.16e g=%e h[20]=%0.16e bL[20]=%0.16e bR[20]=%0.16e \n",j,rhsH[j],j,rhsM[j],fmL[j+1],fmR[j],dx[20],g,h[j + nbc],bL[j + nbc + 1],bR[j + nbc]);
	}
}


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
	//cout<<"cnbsize="<<cnb.size()<<endl;
  for(int i=0;i<cnh.size();i++)
      {hu[i]=cnh[i]*cnu[i];
	// 			cout<<"hu="<<hu[i]<<endl;
}
  int mrk=ord;
	// cout<<"mrk="<<mrk<<endl;
   for(int i=0;i<mrk;i++)
  {
     if((mrk==2)&&(i==0))
     {
       for(int k=0;k<cnh.size();k++)
       {
         h1[k]=cnh[k];
         hu1[k]=hu[k];
       }
     }
		//  cout<<"cnb.size="<<cnb.size()<<endl;
		//  cout<<"cnh.size="<<cnh.size()<<endl;
     SetBDC(cnh,cnu,cnb,A,B,h,u,b);
		//  int  tmp_size = b.size();
		//  for (int q = 0; q < tmp_size; q++){
		// 	 printf("mrk = %d\t [%d] b = %e\n", mrk, q, b[q]);
		//  }
     GetRHS(cnx,h,u,b,rhsH,rhsM,A,B[5]);
     for(int k=0;k<cnh.size();k++)
      {
        cnh[k]=cnh[k]+dt*rhsH[k];
         hu[k]=hu[k]+dt*rhsM[k];

				// cout<<"rhsm="<<rhsM[k]<<"rhsh="<<rhsH[k]<<endl;
      }
     if((mrk==2)&&(i==1))
     {
       for(int k=0;k<cnh.size();k++)
       {
         cnh[k]=0.5*(cnh[k]+h1[k]);
         hu[k]=0.5*(hu1[k]+hu[k]);
       }
     }
     for(int k=0;k<cnu.size();k++)
     {
       cnu[k]=hu[k]/maximum(cnh[k],eps);
			//  cout<<"h["<<k<<"]="<<cnh[k]<<endl;

     }
	//
   }
  hresult.resize(cnh.size(),0.0);
  uresult.resize(cnu.size(),0.0);
  for(int k=0;k<cnh.size();k++)
      hresult[k]=cnh[k];
  for(int k=0;k<cnu.size();k++)
      uresult[k]=cnu[k];
}
 #endif
