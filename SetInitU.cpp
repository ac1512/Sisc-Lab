#include<iostream>
#include <fstream>
#include<vector>
#include<math.h>

using namespace std;

#define pi 3.141592653

template <typename T>
void SetInitU(parameter &p, T &x, T &b,T &h, T &u){
	
	int nx=p.nx;
	int j;
	T xc[nx];
	for(j=0; j<nx;j++)
		xc[j]=0.5*(x[j]+x[j+1]);

	if (p.exl==0){
		for(j=0;j<nx;j++){
			h[j]=0.1+0.0*xc[j];
			u[j]=1.5+0.0*xc[j];
		}
	}else if(p.exl==1){
		double th;
		for(j=0;j<nx;j++){
			th=cos(2*pi*xc[j]);
			h[j]=5.0+exp(th);
			u[j]=sin(th)/h[j];
		}
	}else if(p.exl==3){
		if(p.wb==1){		//test well-balancing
			for(j=0;j<nx;j++){
				if (1 - b[j] > 0)
					h[j] = 1 - b[j];
				else
					h[j] = 0;
				u[j]=0.0*h[j];
			}
		}else{
			int D=1; 
			double delta=0.019, gamma;
			gamma=sqrt(3*delta/(4*D));
			x_a=a2.178272/gamma;
			for(j=0;j<nx;j++){
				if (D + delta / (cosh(gamma*(xc[j] - x_a))*cosh(gamma*(xc[j] - x_a))) - b[j] > 0)
					h[j] = D + delta / (cosh(gamma*(xc[j] - x_a))*cosh(gamma*(xc[j] - x_a))) - b[j];
				else
					h[j] = 0;
				u[j]=sqrt(p.g/D)*delta/ (cosh(gamma*(xc[j] - x_a))*cosh(gamma*(xc[j] - x_a)));
			}
		}
	}else if(p.exl==4||p.exl==41){
		double alpha=p.alpha;
		for(j=0;j<nx;j++){
			h[j]=0;
			if(xc[j]<0)
				h[j]=1+xc[j]*tan(alpha);
			u[j]=0.0;
			
		}
	}else if(p.exl==42||p.exl==41){
		double alpha=p.alpha;
		for(j=0;j<nx;j++){
			h[j]=0;
			if(xc[j]<0)
				h[j]=1+xc[j]*tan(alpha);
			if(h[j]<0)
				h[j]=1+xc[j]*tan(alpha);				
			u[j]=0.0;
		}
	}else if (p.exl==5){
		for(j=0;j<nx;j++){
			h[j]=2+0.0*xc[j]-b[j];
			u[j]=0.0*xc[j];
		}
	}else if (p.exl==6){
		for(j=0;j<nx;j++){
			if(xc>0.5){
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
