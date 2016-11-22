#include<iostream>
#include <fstream>
#include<vector>
#include<math.h>

using namespace std;
template <typename T>
RecCH( T &x, T &h, T &u, T &b, T &wL, T &hL, T &uL, T &bL, T &wR, T &hR, T &uR, T &bR, parameter &p) {
	int nx = p.nx, nbc = p.nbc;
	int j, n=h.size();
	T w[n], dw[nx+nbc+1], db[nx + nbc + 1], du[nx + nbc + 1], dh[nx + nbc + 1];
	for (j = 0; j < n; j++)
		w[j] = h[j] + b[j];
	for (j = nbc - 1; j < (nx + nbc + 1); j++) {
		dw[j] = limit(w[j] - w[j - 1], w[j + 1] - w[j]);
		db[j] = limit(b[j] - b[j - 1], b[j + 1] - b[j]);
		du[j] = limit(u[j] - u[j - 1], u[j + 1] - u[j]);
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
void GetRHS(T &x,T &h, T&u, T &b, T &rhsH, T &rhsM, parameter &p){
	
	int nx=p.nx, nbc=p.nbc;
	int j, n=u.size();

	T hL[n], hR[n], wL[n], wR[n], uL[n], uR[n], bL[n],bR[n];
	RecCH(p, x, h, u, b, wL, hL, uL, bL, wR, hR, uR, bR);


	T bmax[nx + nbc + 1], bo[nx + nbc + 1], wmin[nx + nbc + 1], hLs[nx + nbc + 1], hRs[nx + nbc + 1], dpL[nx + nbc + 1], dpR[nx + nbc + 1];
	for (j = nbc; j < (nx + nbc + 1); j++) {
		bmax[j] = (bL[j] > bR[j]) ? bL[j] : bR[j];
		wmin[j] = (wL[j] < wR[j]) ? wL[j] : wR[j];
		bo[j] = (bmax[j] < wmin[j]) ? bmax[j] : wmin[j];
		hLs[j] = ((wL[j] - bo[j]) < hL[j]) ? (wL[j] - bo[j]) : hL[j];
		hRs[j] = ((wR[j] - bo[j]) < hR[j]) ? (wR[j] - bo[j]) : hR[j];
		dpL[j] = p.g*(hL[j] + hLs[j])*(bo[j]-bL[j]) / 2;
		dpR[j] = p.g*(hR[j] + hRs[j])*(bo[j] - bR[j]) / 2;
	}
	T fh[nx + nbc+1], fm[nx + nbc + 1], fmL[nx + nbc + 1], fmR[nx + nbc + 1];
	GetNumFlux(*(hLs+nbc), *(uL + nbc), *(hRs + nbc), *(uR+ nbc), *(fh + nbc), *(fm + nbc));
	for (j = nbc; j < (nx + nbc+1); j++) {
		fmL[j] = fm[j] + dpL[j];
		fmR[j] = fm[j] + dpR[j];
	}
	T dx[nx];
	for (j = 0; j < nx; j++)
		dx[j] = x[j + 1] - x[j];
	for (j = 0; j < nx-1; j++) {
		rhsH[j] = -1 * (fh[j + nbc + 1] - fh[nbc]) / dx[j];
		rhsM[j]= -1 * (fmL[j + nbc + 1] - fmR[nbc]) / dx[j]- p.g*h[j+nbc]*(bL[j + nbc + 1] - bR[nbc])/dx[j];
	}
}
