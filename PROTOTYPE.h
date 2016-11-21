#ifndef  _PROTOTYPE_H_INCLUDED
#define  _PROTOTYPE_H_INCLUDED

#include <vector>
using namespace std;
class parameter
{
    public:
        void set_param(parameter &p);
        void get_param(parameter &p,vector<int>& a,vector<double>& b);
     private:
    	int ne; //number of cells.
    	int exl;
    	int draw_frames;
    	int ord;   // holds the order of the numerical scheme.
    	int mth_lim; //limiter
    	int mth_rk;  //what is thsi mth_rk?
    	int Flux_method;
    	int wb;
    	int bdd;
    	double tstop;
    	int alpha,IT;
    	double tme;
        double eps;
        double htol,hdry,nbc,cfl,g;

};
#endif
