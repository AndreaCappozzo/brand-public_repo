#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// General functions

// [[Rcpp::export]]
double LogSumExp(arma::colvec logX){
 double a = max(logX);
return(  a + log(accu( exp( logX-a ) )));
}


const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
double dmvnrm_arma(arma::mat x,  
                      arma::rowvec mean,  
                      arma::mat sigma, 
                      bool logd = false) { 
  //int n = x.n_rows;
  int xdim = x.n_cols;
  double out;
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  
    arma::vec z = rooti * arma::trans( x - mean) ;    
    out      = constants - 0.5 * arma::sum(z%z) + rootisum;     

  if (logd == false) {
    out = exp(out);
  }
  return(out);
}





// [[Rcpp::export]]
arma::colvec UPD_Sticks_Beta_cpp(arma::mat AB, 
                             int L, 
                             double alphaDP){
  
  arma::colvec beta_lab           = AB.col(1);
  arma::uvec   inds               = find(beta_lab>0);
  arma::colvec beta_lab_newgroups = beta_lab.elem(inds);
  
  arma::colvec Sticks(L);
  
  if(beta_lab_newgroups.n_rows==0){
    
    for(int l=0; l<(L-1); l++){
      Sticks[l] = R::rbeta(1 , 
                           alphaDP  );  }
  }else{
    for(int l=0; l<(L-1); l++){
      Sticks[l] = R::rbeta(1 + accu( beta_lab_newgroups == (l+1)), 
                           alphaDP + accu(beta_lab_newgroups > (l+1)) );
  }
  }
  
  
  
  Sticks[L-1]=1.;
  return Sticks; 
  }


// [[Rcpp::export]] 
arma::colvec StickBreaker_cpp(arma::colvec V) {
  int N = V.size();
  arma::colvec pi(N), mV2(N);
  arma::colvec mV = 1-V;
  mV2   = exp(arma::cumsum(log(mV)));
  mV2   = arma::shift(mV2,+1); 
  mV2(0) = 1;
  return(V%mV2);
}



// [[Rcpp::export]] 
List Update_Theta_cpp(arma::mat Y, 
                      arma::mat AB, 
                      arma::colvec m_H, double k_H, double v_H, arma::mat S_H, 
                      int L, int p){
  
  
  arma::mat    mu_new(p,L); 
  arma::cube   Sig_new(p,p,L);
  arma::colvec beta_lab = AB.col(1);
  
  for(int l=0; l<L; l++){
    arma::uvec ind_l = find( beta_lab == (l+1) ) ;
    arma::mat  Y_l   = Y.rows(ind_l);
    int n_l          = (Y_l.n_rows);
    
    double k_ = n_l + k_H;
    double v_ = n_l + v_H;
    
    arma::colvec m_(p);
    arma::mat    S_(p,p);
    // Rcout << "ok1 \n"; 
    
    
    
    if(n_l>1){
      //Rcout << Y_l;
      arma::colvec my = mean(Y_l,0).t();
      //Rcout << my;
      m_ = (k_H * m_H + my * n_l) / k_;
      // Rcout << "ok2 \n"; 

      S_ =  S_H + cov(Y_l) * (n_l-1) + 
                  (k_H * n_l)/(k_H+n_l) * (my-m_H)*((my-m_H).t());
      //Rcout << (my-m_H)*((my-m_H).t());
      // Rcout << S_;
      // Rcout << "ok3 \n"; 
    }else if(n_l == 1 ){  
  
    arma::colvec my = Y_l.t();
    m_ = (k_H * m_H + my * n_l) / k_;
    S_ =  S_H + (k_H* n_l)/(k_H+n_l) * (my-m_H)*((my-m_H).t());
    
    }else{
      m_ = m_H;
      S_ = S_H;
    }
    

    // Rcout << "ok4 \n"; 
    
    Sig_new.slice(l) = arma::iwishrnd( (S_), v_);
    // Rcout << "ok5 \n"; 
    
    mu_new.col(l)    = arma::mvnrnd(m_, Sig_new.slice(l) /k_);
    // Rcout << "ok6 \n"; 
    
  }
  
  return List::create(
    _["mu"]      = mu_new,
    _["Sigma"]   = Sig_new);
}













// [[Rcpp::export]] 
List Update_Theta_cpp_TRAIN(arma::mat Y, 
                      arma::mat AB, 
                      arma::mat MU_g,  double k_g, 
                      arma::cube SIGMA_g, double v_g, 
                      int G, int p){
  
  
  arma::mat    Xg(p,G); 
  arma::cube   S2g(p,p,G);
  arma::colvec alpha_lab = AB.col(0);
  
  for(int g=0; g<G; g++){
    arma::uvec ind_g = find( alpha_lab == (g+1) ) ;
    arma::mat  Y_g   = Y.rows(ind_g);
    int n_g          = (Y_g.n_rows);
    
    double k_ = n_g + k_g;
    double v_ = n_g + v_g;
    
    arma::colvec m_(p);
    arma::mat    S_(p,p);
    // Rcout << "ok1 \n"; 
    
    arma::colvec   m_g = MU_g.col(g);
    arma::mat      S_g = SIGMA_g.slice(g);
    
    if(n_g>1){
      //Rcout << Y_l;
      arma::colvec my = mean(Y_g,0).t();
      //Rcout << my;
      m_ = (k_g * m_g + my * n_g) / k_;
      // Rcout << "ok2 \n"; 
      
      S_ =  S_g + cov(Y_g) * (n_g-1) + 
        (k_g * n_g)/(k_g+n_g) * (my-m_g)*((my-m_g).t());
      //Rcout << (my-m_H)*((my-m_H).t());
      // Rcout << S_;
      // Rcout << "ok3 \n"; 
    }else if(n_g == 1 ){  
      
      arma::colvec my = Y_g.t();
      m_ = (k_g * m_g + my * n_g) / k_;
      S_ =  S_g + (k_g* n_g)/(k_g+n_g) * (my-m_g)*((my-m_g).t());
      
    }else{
      m_ = m_g;
      S_ = S_g;
    }
    
    
    // Rcout << "ok4 \n"; 
    
    S2g.slice(g) = arma::iwishrnd( (S_), v_);
    // Rcout << "ok5 \n"; 
    
    Xg.col(g)    = arma::mvnrnd(m_, S2g.slice(g) /k_);
    // Rcout << "ok6 \n"; 
    
  }
  
  return List::create(
    _["mu"]      = Xg,
    _["Sigma"]   = S2g);
}


// ISH-JAMES sampler

// [[Rcpp::export]]
arma::mat Upd_alphabeta_cpp(arma::mat Y, 
                            arma::mat mu, arma::cube Sigma,
                            arma::colvec pidir,
                            arma::colvec omega, 
                            int J, int n,  int L,
                            arma::colvec poss_lab,
                            arma::mat x_bar,
                            arma::cube S2_j){
  arma::mat res_AB(n,2);
  res_AB.fill(0);
  double pi0 = pidir(J);
  pidir.shed_row(J); // removed last element
  //Rcout << "ok1";
 for(int ii=0; ii<n; ii++) {
   arma::colvec likold(J);
   arma::colvec liknew(L);
   for(int j=0; j<J; j++){
     likold(j) = dmvnrm_arma(Y.row(ii), x_bar.col(j).t(), S2_j.slice(j), true);
   }
   for(int l=0; l<L; l++){
     liknew(l) = dmvnrm_arma(Y.row(ii), mu.col(l).t() ,Sigma.slice(l),true);
   }
   //Rcout << "ok2";
   
   arma::colvec oldprob = log(pidir) + likold;
   arma::colvec newprob = log(pi0) + log(omega) + liknew;
                            
   arma::colvec r = join_cols(oldprob, newprob);
                r = exp( r - max(r));
   int ind = RcppArmadillo::sample(poss_lab, 1, TRUE, r)[0];
   //Rcout << "ok3";
   
   if (ind <= J) {
     res_AB(ii,0) = ind;
   } else{
     res_AB(ii,1) = ind - J;
   }
 }
 return(res_AB);
 }








// SLICE sampler

// [[Rcpp::export]]
arma::mat Upd_alphabeta_cpp_SLICE(
                            arma::mat Y, 
                            arma::mat mu, arma::cube Sigma,
                            arma::colvec pidir,
                            arma::colvec Ui, arma::colvec xi,
                            arma::colvec omega, 
                            int J, int n,  int L,
                            arma::colvec poss_lab,
                            arma::mat x_bar,
                            arma::cube S2_j){
  arma::mat res_AB(n,2);
  res_AB.fill(0);
  double pi0 = pidir(J);
  pidir.shed_row(J); // removed last element
  //Rcout << "ok1";
  arma::colvec likold(J);
  arma::colvec liknew(L);
  
  for(int ii=0; ii<n; ii++) {
    for(int j=0; j<J; j++){
      likold(j) = dmvnrm_arma(Y.row(ii), x_bar.col(j).t(), S2_j.slice(j), true);
    }
    for(int l=0; l<L; l++){
      if(Ui[ii] < xi[l]){
      liknew(l) =  dmvnrm_arma(Y.row(ii), mu.col(l).t() ,Sigma.slice(l),true) - log(xi[l]);  
      }else{
      liknew(l) = - arma::datum::inf;  
      }
    }

    arma::colvec oldprob = log(pidir) + likold;
    arma::colvec newprob = log(pi0) + log(omega) + liknew ;
    
    //Rcout << "$"<< oldprob << "\n";
    //Rcout << "%"  << newprob << "\n";
    
    arma::colvec r = join_cols(oldprob, newprob);
    //Rcout << oldprob << "--1\n" << newprob << "--1\n" << r << "--1\n";
    r = exp( r - LogSumExp(r) );
    
    //if(ii==0){Rcout << r << "---\n";}
    
    int ind = RcppArmadillo::sample(poss_lab, 1, TRUE, r)[0];
    //Rcout << "ok3";
    
    if (ind <= J) {
      res_AB(ii,0) = ind;
    } else{
      res_AB(ii,1) = ind - J;
    }
  }
  return(res_AB);
}







// [[Rcpp::export]]
arma::colvec Upd_Zeta_cpp(
    arma::mat Y, 
    arma::mat mu_all, 
    arma::cube Sigma_all,
    arma::colvec Uitilde, 
    arma::colvec xitilde,
    arma::colvec log_omegatilde, 
    int J, int n,  int L,
    arma::colvec poss_lab){
  int toto = L+J;
  arma::colvec res_Z(n);  res_Z.fill(0);
  arma::colvec likel(toto);
  
  for(int ii=0; ii<n; ii++) {
    
    for(int h=0; h<(toto); h++){
      if(Uitilde[ii] < xitilde[h]){
        likel(h) =  dmvnrm_arma(Y.row(ii), mu_all.col(h).t() ,Sigma_all.slice(h),true) - log(xitilde[h]);  
      }else{
        likel(h) = - arma::datum::inf;
      }
    }
    //Rcout << ii << "\n +++ \n"<< likel << "\n -- \n";
    //Rcout << "ok1";
    likel = log_omegatilde + likel ;
    arma::colvec r = exp( likel - LogSumExp(likel) );
    //Rcout << ""<< r <<"";
    res_Z[ii] = RcppArmadillo::sample(poss_lab, 1, TRUE, r)[0]; // the error happens here
    //Rcout << "ok3";
  }
  return(res_Z);
}



// [[Rcpp::export]]
arma::colvec sasa(arma::colvec X, arma::colvec poss_lab, int N){
return(RcppArmadillo::sample(poss_lab, N, TRUE, X)); // the error happens here
}