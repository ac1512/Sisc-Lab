#include<iostream>
#include <fstream>
#include<vector>
#include<math.h>

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

template <typename T>
void RecCH( T &x, T &h, T &u, T &b, T &wL, T &hL, T &uL, T &bL, T &wR, T &hR, T &uR, T &bR, int &a) {
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

template <typename T>
void GetRHS(T &x,T &h, T&u, T &b, T &rhsH, T &rhsM, int &a, double &paramb){
	
	int nx=a[1], nbc=a[12];
	int j, n=sizeof(u);

	double hL[n], hR[n], wL[n], wR[n], uL[n], uR[n], bL[n],bR[n];
	RecCH(x, h, u, b, wL, hL, uL, bL, wR, hR, uR, bR, a);


	double bmax[nx + nbc + 1], bo[nx + nbc + 1], wmin[nx + nbc + 1], hLs[nx + nbc + 1], hRs[nx + nbc + 1], dpL[nx + nbc + 1], dpR[nx + nbc + 1];
	for (j = nbc; j < (nx + nbc + 1); j++) {
		bmax[j] = (bL[j] > bR[j]) ? bL[j] : bR[j];
		wmin[j] = (wL[j] < wR[j]) ? wL[j] : wR[j];
		bo[j] = (bmax[j] < wmin[j]) ? bmax[j] : wmin[j];
		hLs[j] = ((wL[j] - bo[j]) < hL[j]) ? (wL[j] - bo[j]) : hL[j];
		hRs[j] = ((wR[j] - bo[j]) < hR[j]) ? (wR[j] - bo[j]) : hR[j];
		dpL[j] = paramb[5]*(hL[j] + hLs[j])*(bo[j]-bL[j]) / 2;
		dpR[j] = paramb[5]*(hR[j] + hRs[j])*(bo[j] - bR[j]) / 2;
	}
	double fh[nx + nbc+1], fm[nx + nbc + 1], fmL[nx + nbc + 1], fmR[nx + nbc + 1];
	GetNumFlux(*(hLs+nbc), *(uL + nbc), *(hRs + nbc), *(uR+ nbc), *(fh + nbc), *(fm + nbc), a,b);
	for (j = nbc; j < (nx + nbc+1); j++) {
		fmL[j] = fm[j] + dpL[j];
		fmR[j] = fm[j] + dpR[j];
	}
	double dx[nx];
	for (j = 0; j < nx; j++)
		dx[j] = x[j + 1] - x[j];
	for (j = 0; j < nx-1; j++) {
		rhsH[j] = -1 * (fh[j + nbc + 1] - fh[nbc]) / dx[j];
		rhsM[j]= -1 * (fmL[j + nbc + 1] - fmR[nbc]) / dx[j]- paramb[5]*h[j+nbc]*(bL[j + nbc + 1] - bR[nbc])/dx[j];
	}
}
