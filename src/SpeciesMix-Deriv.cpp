#include"SpeciesMix-Deriv.h"

extern "C" 
{ 
  SEXP SpeciesMix(SEXP R_pars, SEXP R_y, SEXP R_X, SEXP R_ID,SEXP R_tau, SEXP R_gradient){
    // y is response
    //X is design matrix
    // ID is vector length(y) with species names for each observation
    // ID_names is vector length S of unique ID
    // tau is easy pass out of matrix of tau
    // estpi is easy pass out of pi
    double *pars=NULL, *y=NULL, *X=NULL, *estpi=NULL, *tau=NULL, *pi=NULL, *r_pointer=NULL, *gradient=NULL;
    int *ID=NULL;
    double logl=0;
    int S,G, Xr, Xc, lpar,s,i ,ly;
    SEXP dimensions, R_logl;
    double abstol,reltol;
    int fncount, grcount, ifail,trace,nREPORT;
  
  
    Optimise_data data;

    //R_pars has length = length(pi) + length( G* (number of covariates + intercept))
    // ordered so that first G elements are pi and remaining elements are coefficents. all coefficents from same group are together
    lpar = LENGTH(R_pars);
    ly = LENGTH(R_y);
    //    vector<double> prams( lpar );
    vector< double > params(lpar);
    pars=REAL(R_pars);
    y=REAL(R_y);
    X=REAL(R_X);
    tau=REAL(R_tau);
    ID=INTEGER(R_ID);
    gradient=REAL(R_gradient);

 
    for(i=0;i<lpar;i++){ params.at(i) = pars[i]; }



    dimensions=getAttrib(R_X, R_DimSymbol);
    Xr=INTEGER(dimensions)[0];
    Xc=INTEGER(dimensions)[1];

    dimensions=getAttrib(R_tau, R_DimSymbol);
    S=INTEGER(dimensions)[0];
    G=INTEGER(dimensions)[1];

    
    data.SetVars(y, X, ID, S, G,  Xr, Xc,  lpar, ly, tau);
    //    Rprintf("%d,%d,%d,%d,%d\n",ly,Xr,Xc,S,G);
    vector< double > pll( 1, 0);
    vector< double > grd(lpar,0);
   
     // set up variables for optimisation
    abstol = 1e-8;
    reltol = 1e-8;
    nREPORT = 5;//how often to report
    fncount=0;
    grcount=0;
    ifail=0;
    trace=1;
    vector <int> mask (lpar,1); 
    vector < double > logl_out(lpar,0);

        vmmin(lpar, pars, &logl_out.at(0), optimise_function, gradient_function, 1000, trace, &mask.at(0),  abstol,  reltol,  nREPORT, &data, &fncount, &grcount, &ifail);
    
    std::cout << "ifail = " << ifail << "\n" ;
    logl=1;
    logl = optimise_function(lpar, pars, &data);

    //std::cout << "Optimise = " << logl << "\n";
    gradient_function(lpar, pars, &grd.at(0),&data);
    // std::cout << "Gradient  = ";
    for(i=0;i<lpar;i++){
      //  std::cout << i<< " : " << grd.at(i) << ", ";
      gradient[i] = grd.at(i);
    }

    // std::cout << "\n";


    R_logl = allocVector(REALSXP,1);
    r_pointer = REAL(R_logl);
    *r_pointer = logl;
    return(R_logl);
  }

}

double optimise_function(int n, double *pars, void *ex){

  Optimise_data *data = (Optimise_data *) ex;
  //Optimise_data data =  * (Optimise_data *) ex;
  
  vector<double> logl(1,0);
  int i;
  vector<double> x (n,0);
  
  for(i=0;i<n;i++) x.at(i) = pars[i];

  logl = calc_logl(x,*data);

  //logl = data->F.Forward(0,x);

  return(0-logl.at(0));

}

void gradient_function(int n, double *pars, double *gr, void *ex ){
  Optimise_data *data = (Optimise_data *) ex;
  //Optimise_data data = *(Optimise_data *) ex;
  vector<double> ad_g(n,0);
  vector<double> x (n,0);
  int i,G,g,j,Xc,s,S;
  double add_log_trans=0;
  vector<double> logl(1,0);

  G=data->G;
  S=data->S;
  vector<double> pi_mat_deriv(G*(G-1),0);
  vector<double> dl_dpi(G,0);

  Xc=data->Xc;
  for(i=0;i<n;i++) x.at(i) = pars[i];
 
  logl = calc_logl(x,*data);

  for(g=0;g<(G-1);g++){ 
      add_log_trans+=exp(x.at(g)); //add up transformed pi's
  }
  add_log_trans+=1;
 
  for(g=0;g<G;g++){ //GO through all pi's

      for(j=0;j<Xc;j++){
	for(s=0;s<S;s++){
	  ad_g.at((G-1)+ MAT_RF(g,j,G)) += exp( -1*data->species_l_contrib.at(s) + log(data->parpi.at(g)) + data->sum_f_species.at(MAT_RF(g,s,G))) *data->deriv_f_B.at(MAT_3D(g,j,s,G,Xc));
	  if(j==0) dl_dpi.at(g) += exp(  -1*data->species_l_contrib.at(s)+ data->sum_f_species.at(MAT_RF(g,s,G)));
	}
      }

      for(i=0;i<(G-1);i++){ // go through eta's
	if(g<(G-1)){
	  if(i==g){
	    pi_mat_deriv.at(MAT_RF(i,g,(G-1))) = exp(x.at(i))/add_log_trans - exp(2*x.at(i))/(add_log_trans*add_log_trans);// diag
	    pi_mat_deriv.at(MAT_RF(i,(G-1),(G-1))) += pi_mat_deriv.at(MAT_RF(i,g,(G-1)));
	  }else{
	    pi_mat_deriv.at(MAT_RF(i,g,(G-1))) = -exp(x.at(i))*exp(x.at(g)) / (add_log_trans*add_log_trans); //off-diag
	    pi_mat_deriv.at(MAT_RF(i,(G-1),(G-1))) += pi_mat_deriv.at(MAT_RF(i,g,(G-1)));
	  }
	}
      }

  }
  for(i=0;i<(G-1);i++) pi_mat_deriv.at(MAT_RF(i,(G-1),(G-1))) *= -1;



    for(i=0;i<(G-1);i++){
      for(g=0;g<G;g++){
	ad_g.at(i) += dl_dpi.at(g)* pi_mat_deriv.at(MAT_RF(i,g,(G-1)));
      }
    }

  for(i=0;i<n;i++) gr[i] = 0-ad_g.at(i); 

}


 vector <double> calc_logl(const vector<double> &pars, Optimise_data &data){
  int G,S,Xr,Xc;
  G=data.G;
  S=data.S;
  Xr=data.Xr;
  Xc=data.Xc;
  vector< double > estpi(G-1,0); //vector to hold pi's
  vector< double > coef(Xc*G,0); //vector to hold all coefficents
  //vector< double > logl( S * G ,0); //output log likelihood
  vector< double > logl(1 ,0), tlog(1 ,0); //output log likelihood
  //  double sumlogl=0;
  int i,s,g,j;
  int start, end;
  double tmp;
  
  for(i=0;i<(G-1);i++){
    estpi.at(i) = pars.at(i); //G-1 values for pi
  }

  additive_logistic(estpi,1); // additive logistic transfor on pi;

  for(i=0;i<G;i++){
    data.parpi.at(i) = estpi.at(i);
    for(s=0;s<S;s++){ 
      data.sum_f_species.at(MAT_RF(i,s,G))=0;
      for(j=0;j<Xc;j++) data.deriv_f_B.at(MAT_3D(i,j,s,G,Xc)) =0;
    }
  }
 
 
  for(i=(G-1);i<data.lpar;i++){coef.at(i-(G-1))= pars.at(i);}

  logl.at(0) = 0;

  for(s=0;s<S;s++){
      start = data.StartEndPos.at(s*2);
      end = data.StartEndPos.at(s*2+1);

      tlog.at(0) = like_function(estpi, coef,data.y,data.X,Xr,Xc,start,end, data.tau, s ,data.sum_f_species,data.deriv_f_B);
      logl.at(0)+= tlog.at(0);
      data.species_l_contrib.at(s) = tlog.at(0);
  }
  return(logl);

}


double like_function(vector< double > &estpi, vector < double > &coef, const double *y, const double *X, int Xr, int Xc, int start, int end, double *tau, int s, vector<double> &sum_f_species, vector<double> &deriv_f_B){
  int len,i,j,G,g;
  len=end-start+1;
  //vector< AD<double > > p(len,0);
  vector< double  > p(1,0);
  double lpre=0, eps=0,glogl=0;
  G = estpi.size();
  vector< double  > sump(G,0);

  for(g=0;g<G;g++){
    for(i=0;i<len;i++){
      lpre=0;
      for(j=0;j<Xc;j++){ 
	lpre +=  X[MAT_RF(i,j,Xr)] * coef[MAT_RF(g,j,G)];
      }
      p.at(0)=inv_link_function(lpre,0);  // mu for each archetype at site k
 
      for(j=0;j<Xc;j++) {deriv_f_B.at(MAT_3D(g,j,s,G,Xc)) +=   (y[start+i] - p.at(0)) * X[MAT_RF(i,j,Xr)] ;} 

      if(y[start+i]==0) p.at(0) = 1-p.at(0);

      sump.at(g) += log(p.at(0));
  
    }
    if(g==1) eps = sump.at(g);
    if(sump.at(g) > eps) eps = sump.at(g);
    
  }

  for(g=0;g<G;g++){
    sum_f_species.at(MAT_RF(g,s,G)) = sump.at(g);
    glogl+= estpi.at(g)*exp(sump.at(g) - eps);
  }
  // code for taus can go here
  glogl = log(glogl) + eps;
 
  return(glogl);

}

double link_function(double p, int link_type){
  // using logit link function
  if(link_type==0) return(log(p)-log(1-p));
}
double inv_link_function(double lpre, int link_type){
  // using logit link function
  if(link_type==0) return(exp(lpre)/(1+exp(lpre)));
}

void additive_logistic(vector< double > &x,int inv){
  int i; 
  // inv == 1 gives transfornmation
  //inv = 0 give inverse transformatio
  if(inv==0){
    for(i=0;i<x.size();i++) x.at(i) = log(x.at(i)/x.back());
    return;
  }
    
  vector< double > xt (x.size(),0);
  double sumx=0, sumxt=0;

  for(i=0;i<x.size();i++){
    xt.at(i) = exp(x.at(i));
    sumx+=xt.at(i);
  }
  for(i=0;i<x.size();i++){
    x.at(i) = xt.at(i)/(1+sumx);
    sumxt+=x.at(i);
  }
  x.push_back(1-sumxt);

}

Optimise_data::Optimise_data(){}
Optimise_data::~Optimise_data(){}

void Optimise_data::SetVars(double *ty, double *tX, int *tID, int tS, int tG, int tXr, int tXc, int tlpar, int tly, double *ttau){
  int i, s, j;
  tau = ttau;
  y=ty;
  X=tX;
  ID=tID;
  S=tS;
  G=tG;
  Xr=tXr;
  Xc=tXc;
  lpar=tlpar;
  ly = tly;

  //vector<int> StartEndPos(S*2);

  StartEndPos.push_back(0);
  //s=1;
  for(i=0;i<G;i++){ // vectors of length G
    parpi.push_back(0);
    for(s=0;s<S;s++){ 
      sum_f_species.push_back(0);
      for(j=0;j<Xc;j++) deriv_f_B.push_back(0);
    }
  }

  for(s=0;s<S;s++) species_l_contrib.push_back(0);

  s=ID[0];

  for(i=0;i<ly;i++) // find start & end positions of each species data
    if(ID[i]!=s){
      StartEndPos.push_back(i-1);
      //s++; // index to next species
      //std::cout << s << "," << ID[i] << "," << i <<"\n" ; 
      s=ID[i];
      StartEndPos.push_back(i); // next species start pos
    }
  
  StartEndPos.push_back(ly-1);  //add final position
 
}




extern "C" 
{ 
  SEXP Calculate_Gradient(SEXP R_pars, SEXP R_y, SEXP R_X, SEXP R_ID,SEXP R_tau, SEXP R_gradient){
    // y is response
    //X is design matrix
    // ID is vector length(y) with species names for each observation
    // ID_names is vector length S of unique ID
    // tau is easy pass out of matrix of tau
    // estpi is easy pass out of pi
    double *pars=NULL, *y=NULL, *X=NULL, *estpi=NULL, *tau=NULL, *pi=NULL, *r_pointer=NULL, *gradient=NULL;
    int *ID=NULL;
    double logl=0;
    int S,G, Xr, Xc, lpar,s,i ,ly;
    SEXP dimensions, R_logl;
    double abstol,reltol;
    int fncount, grcount, ifail,trace,nREPORT;
  
  
    Optimise_data data;

    //R_pars has length = length(pi) + length( G* (number of covariates + intercept))
    // ordered so that first G elements are pi and remaining elements are coefficents. all coefficents from same group are together
    lpar = LENGTH(R_pars);
    ly = LENGTH(R_y);
    //    vector<double> prams( lpar );
    vector< double > params(lpar);
    pars=REAL(R_pars);
    y=REAL(R_y);
    X=REAL(R_X);
    tau=REAL(R_tau);
    ID=INTEGER(R_ID);
    gradient=REAL(R_gradient);

     for(i=0;i<lpar;i++){ params.at(i) = pars[i]; }



    dimensions=getAttrib(R_X, R_DimSymbol);
    Xr=INTEGER(dimensions)[0];
    Xc=INTEGER(dimensions)[1];

    dimensions=getAttrib(R_tau, R_DimSymbol);
    S=INTEGER(dimensions)[0];
    G=INTEGER(dimensions)[1];

    
    data.SetVars(y, X, ID, S, G,  Xr, Xc,  lpar, ly, tau);
    //    Rprintf("%d,%d,%d,%d,%d\n",ly,Xr,Xc,S,G);
    vector< double > pll( 1, 0);
    vector< double > grd(lpar,0);
   
     // set up variables for optimisation
    abstol = 1e-8;
    reltol = 1e-8;
    nREPORT = 5;//how often to report
    fncount=0;
    grcount=0;
    ifail=0;
    trace=1;
    vector <int> mask (lpar,1); 
    vector < double > logl_out(lpar,0);

    logl=1;
    logl = optimise_function(lpar, pars, &data);

    //std::cout << "Optimise = " << logl << "\n";
    gradient_function(lpar, pars, &grd.at(0),&data);
    // std::cout << "Gradient  = ";
    for(i=0;i<lpar;i++){
      //  std::cout << i<< " : " << grd.at(i) << ", ";
      gradient[i] = grd.at(i);
    }

    // std::cout << "\n";


    R_logl = allocVector(REALSXP,1);
    r_pointer = REAL(R_logl);
    *r_pointer = logl;
    return(R_logl);
  }

}
