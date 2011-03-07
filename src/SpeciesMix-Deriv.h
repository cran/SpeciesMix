#include<R.h>
#include<Rinternals.h>
#include<Rmath.h>
#include<R_ext/Applic.h>
#undef length

#include <vector>
#include<algorithm>
#include <iostream>

#define MAT_RF(i,j,nx) i+nx*j
#define MAT_3D(i,j,k,nx,ny)  i + nx*j + k*(nx*ny)
using std::vector;         // use vector as abbreviation for std::vector


class Optimise_data{
public:
  Optimise_data();
  ~Optimise_data();
  
  double *y, *X, *tau;
  int *ID;
  int S, G,  Xr, Xc, lpar, ly;
  vector<int> StartEndPos;
  
  void SetVars(double *ty, double *tX, int *tID, int tS, int tG, int tXr, int tXc, int tlpar, int tlobs, double *ttau);

  // holders for derivitive calculation values

  double logl;
  vector<double> sum_f_species; // log( f(yi,Bi)) S*G long
  vector<double> deriv_f_B;  //  d(log f(yi,Bi)) / d( Bi) G*Xc*S long
  vector<double> parpi; // vector containing calculated pi's G long
  vector<double> species_l_contrib; // vector of each species likelihood contribution // S long

};
extern "C"  SEXP SpeciesMix(SEXP R_pars, SEXP R_y, SEXP R_X, SEXP R_ID,SEXP R_tau, SEXP R_gradient);

extern "C"  SEXP Calculate_Gradient(SEXP R_pars, SEXP R_y, SEXP R_X, SEXP R_ID,SEXP R_tau, SEXP R_gradient);

double optimise_function(int n, double *pars, void *ex);
void gradient_function(int n, double *pars, double *gr, void *ex );

vector<double> calc_logl(const vector<double> &pars, Optimise_data &data);
double like_function(vector <double> &estpi, vector < double > &coef, const double *y, const double *X, int Xr, int Xc, int start, int end, double *tau, int s, vector<double> &sum_f_species, vector<double> &deriv_f_B);
double link_function(double p, int link_type);
double inv_link_function(double lpre, int link_type);

void additive_logistic(vector< double> &x,int inv);
