#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
double LogSumExp(arma::colvec logX){
  double a = max(logX);
  return(  a + log(accu( exp( logX-a ) )));
}

 // [[Rcpp::export]]
 arma::colvec trial(arma::colvec X, arma::colvec M, arma::colvec S){
 return(log(normpdf(X, M, sqrt(S))));
 }
 

const double log2pi = std::log(2.0 * M_PI);

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
  mV2   =arma::shift(cumprod(mV),+1); 
  mV2(0) = 1;
  return(V%mV2);
}



// [[Rcpp::export]]
arma::mat Upd_alphabeta_t_cpp(arma::mat Y, 
                              arma::mat effe_train,
                              arma::mat effe_test,
                              arma::mat sigma_train,
                              arma::mat sigma_test,
                              arma::colvec pidir,
                              arma::colvec omega, 
                              int J, int n,  int L, int t,
                              arma::colvec poss_lab){
  arma::mat res_AB(n,2);
  res_AB.fill(0);
  double pi0 = pidir(J);
  pidir.shed_row(J); // removed last element
  //Rcout << "ok1";
  arma::colvec likold(J);
  arma::colvec liknew(L);
  
  for(int  ii=0; ii<n; ii++) {
    for(int j=0; j<J ; j++){
      likold(j) = arma::accu(log(normpdf(Y.col(ii), effe_train.col(j), sqrt(sigma_train.col(j) ) )));
    }
    for(int l=0; l<L; l++){
      
      liknew(l) = arma::accu(log(normpdf(Y.col(ii), effe_test.col(l),  sqrt(sigma_test.col(l))  )));
      
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


// [[Rcpp::export]] 
arma::colvec DN(arma::colvec X, arma::colvec M, arma::colvec S){
  return(arma::normpdf(X,M,S));
}



// [[Rcpp::export]]
arma::mat Upd_alphabeta_t_cpp_SLICE(arma::mat Y, 
                                    arma::mat effe_train,
                                    arma::mat effe_test,
                                    arma::mat sigma_train,
                                    arma::mat sigma_test,
                                    arma::colvec pidir,
                                    arma::colvec omega, 
                                    arma::colvec Ui,
                                    arma::colvec xi,
                                    int J, int n,  int L,
                                    arma::colvec poss_lab){
  arma::mat res_AB(n,2);
  res_AB.fill(0);
  double pi0 = pidir(J);
  pidir.shed_row(J); // removed last element
  //Rcout << "ok1";
  arma::colvec likold(J);
  arma::colvec liknew(L);
  
  for(int ii=0; ii<n; ii++) {
    for(int j=0; j<J; j++){
      
      likold(j) = arma::accu(log(normpdf( Y.col(ii), 
                                 effe_train.col(j), 
                                 sqrt(sigma_train.col(j)))
      ));
      
    }
    for(int l=0; l<L; l++){
      if(Ui[ii] < xi[l]){
        liknew(l) = arma::accu(log(normpdf(Y.col(ii), effe_test.col(l), sqrt(sigma_test.col(l)))));  
      }else{
        liknew(l) = -arma::datum::inf;  
      }
    }
    
    arma::colvec oldprob = log(pidir) + likold;
    arma::colvec newprob = log(pi0)   + log(omega) + liknew;
    
    //  Rcout<< oldprob;
    //  Rcout<< newprob;  
    arma::colvec r = join_cols(oldprob, newprob);
    r = exp( r - LogSumExp(r) );
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


/*
 // [[Rcpp::export]] 
 arma::mat Update_effe_test_cpp(arma::mat Y, arma::mat alphabeta, 
 int L, int t, arma::colvec tau_k, arma::colvec R, 
 arma::colvec sigma_test,arma::colvec ak){
 
 arma::mat effe_test_new(t,L);
 arma::colvec beta_lab = alphabeta.col(1);
 arma::colvec k_(t), m_(t), Sy(t);
 
 for(int l=0; l<L; l++){
 
 arma::uvec ind_l = find( beta_lab == (l+1) ) ;
 arma::mat  Y_l   = Y.cols(ind_l);
 
 int n_l          = (Y_l.n_cols); 
 
 if(n_l>1){
 
 Sy  = arma::sum(Y_l,1);  
 k_  = n_l / sigma_test[l] + 1/(tau_k[l]*R);
 m_  = (Sy/sigma_test[l]+ ak[l]/tau_k[l]) % pow(k_,-1);
 
 effe_test_new.col(l)  = m_ + as<arma::vec>(rnorm(t))%arma::sqrt(pow(k_,-1));
 
 
 }else if( n_l==1){
 
 Sy  = (Y_l);  
 k_  = n_l / sigma_test[l] + 1/(tau_k[l]*R);
 m_  = (Sy/sigma_test[l]+ ak[l]/tau_k[l]) % pow(k_,-1);
 
 effe_test_new.col(l)  = m_ + as<arma::vec>(rnorm(t))%arma::sqrt(pow(k_,-1));
 
 
 }else{
 effe_test_new.col(l) = ak[l] + as<arma::vec>(rnorm(t)) % arma::sqrt(tau_k[l]*R);
 }
 
 }
 return(effe_test_new);
 }
 */

// [[Rcpp::export]] 
arma::mat Update_effe_test_t_cpp(arma::mat Y, arma::mat alphabeta, 
                                 int L, int t, arma::colvec tau_k, arma::colvec R, 
                                 arma::mat sigma_test,arma::colvec ak){
  
  arma::mat effe_test_new(t,L);
  arma::colvec beta_lab = alphabeta.col(1);
  arma::colvec k_(t), m_(t), Sy(t);
  

  
  for(int l=0; l<L; l++){
    
    arma::uvec ind_l = find( beta_lab == (l+1) ) ;
    arma::mat  Y_l   = Y.cols(ind_l);
    
    arma::colvec tau_k_times_R = (tau_k[l]*R);
    
    int n_l          = (Y_l.n_cols); 
    
    if(n_l>1){
      
      Sy  = arma::sum(Y_l,1);  
      k_  = n_l * pow(sigma_test.col(l),-1) + 1/tau_k_times_R;
      //m_  = (Sy % pow(sigma_test.col(l),-1) + ak[l]/tau_k[l]) % pow(k_,-1);
      
      // post cap's review
      m_                    = (Sy % pow(sigma_test.col(l),-1) + ak[l]/tau_k_times_R) % pow(k_,-1);
      ////////////////////
      effe_test_new.col(l)  = m_ + as<arma::vec>(rnorm(t))%arma::sqrt(pow(k_,-1));
      
      
    }else if( n_l==1){
      
      Sy  = (Y_l);  
      k_  = n_l * pow(sigma_test.col(l),-1) + 1/tau_k_times_R;
      //m_  = (Sy % pow(sigma_test.col(l),-1)+ ak[l]/tau_k[l]) % pow(k_,-1);
      m_  = (Sy % pow(sigma_test.col(l),-1)+ ak[l]/tau_k_times_R) % pow(k_,-1);
      effe_test_new.col(l)  = m_ + as<arma::vec>(rnorm(t))%arma::sqrt(pow(k_,-1));
      
      
    }else{
      effe_test_new.col(l) = ak[l] + as<arma::vec>(rnorm(t)) % arma::sqrt(tau_k_times_R);
    }
    
  }
  return(effe_test_new);
}



// [[Rcpp::export]] 
arma::mat Update_effe_train_t_cpp(arma::mat Y, arma::mat alphabeta, 
                                 int G, int t,  
                                 arma::mat f_bar_g,
                                 arma::mat sigma_train,
                                 double KAPPAG){
  
  arma::mat    effe_train_new(t,G);
  arma::colvec alpha_lab = alphabeta.col(0);
  arma::colvec k_(t), m_(t), Sy(t);
  
  
  
  for(int g=0; g<G; g++){
    
    arma::uvec ind_g = find( alpha_lab == (g+1) ) ;
    arma::mat  Y_g   = Y.cols(ind_g);
    
    
    int n_g          = (Y_g.n_cols); 
    
    if(n_g>1){
      
      Sy  = arma::sum(Y_g,1);  
      k_  = n_g * pow( sigma_train.col(g), -1 ) + 1 / KAPPAG;
      
      m_  = (Sy % pow(sigma_train.col(g),-1) + f_bar_g.col(g)/KAPPAG) % pow(k_,-1);
      ////////////////////
      effe_train_new.col(g)  = m_ + as<arma::vec>(rnorm(t))%arma::sqrt(pow(k_,-1));
      
      
    }else if( n_g==1){
      
      Sy  = (Y_g);  
      k_  = n_g * pow(sigma_train.col(g),-1) +  1 / KAPPAG;;
      m_  = (Sy % pow(sigma_train.col(g),-1) + f_bar_g.col(g)/KAPPAG) % pow(k_,-1);
      effe_train_new.col(g)  = m_ + as<arma::vec>(rnorm(t)) % arma::sqrt( pow( k_,- 1));
      
      
    }else{
      effe_train_new.col(g) = f_bar_g.col(g) + as<arma::vec>(rnorm(t)) * pow(KAPPAG,.5);
    }
    
  }
  return(effe_train_new);
}


// [[Rcpp::export]] 
arma::mat Update_sigma_test_t_cpp(arma::mat Y, arma::mat alphabeta, 
                                  arma::mat effe_test,
                                  int L, int t, 
                                  double asig, double bsig){
  
  arma::mat sigma_test_new(t,L);
  arma::colvec beta_lab = alphabeta.col(1);
  
  for(int l=0; l<L; l++){
    
    arma::uvec ind_l = find( beta_lab == (l+1) ) ;
    arma::mat  Y_l   = Y.cols(ind_l);
    
    int n_l          = (Y_l.n_cols);
    
    
    if(n_l>1){
      
      arma::mat fake_FL = arma::repelem(effe_test.col(l),1,n_l);
      arma::mat diff2 = pow((Y_l - fake_FL),2);
      
      for( int small_t=0; small_t<t; small_t++){
        double B = accu( diff2.row(small_t) );
        double bpost  = .5  * B + bsig;
        double apost  = n_l * 1/2 + asig;
        sigma_test_new(small_t,l) = 1/ rgamma(1,apost,1/(bpost))[0];
      }
      
    }else if( n_l==1){
      
      //  Rcout << "2.0";
      
      arma::colvec diff2 = pow( (Y_l- effe_test.col(l)) ,2 ) ;
      //  Rcout << "2.1";
      
      for( int small_t=0; small_t<t; small_t++){
        
        double bpost  = .5  *  diff2(small_t) + bsig;
        //      Rcout << "2.2";
        
        double apost  = 1/2 + asig;
        sigma_test_new(small_t,l) = 1/ rgamma(1,apost,1/(bpost))[0];
      }
      //   Rcout << "2";
      
    }else{
      for( int small_t=0; small_t<t; small_t++){
        sigma_test_new(small_t,l) = 1/ rgamma(1,asig,1/(bsig))[0];
      }
      //   Rcout << "3";
    }
  }
  return(sigma_test_new);
}



// [[Rcpp::export]] 
arma::mat Update_sigma_train_t_cpp(arma::mat Y, arma::mat alphabeta, 
                                   arma::mat effe_train,
                                   int G, int t, 
                                   arma::mat a_priorG, arma::mat b_priorG){
  
  arma::mat sigma_train_new(t,G);
  arma::colvec alpha_lab = alphabeta.col(0);
  
  for(int g=0; g<G; g++){
    
    arma::uvec ind_g = find( alpha_lab == (g+1) ) ;
    arma::mat  Y_g   = Y.cols(ind_g);
    
    int n_g          = (Y_g.n_cols);
    
    if(n_g>1){
      
      arma::mat fake_FG = arma::repelem(effe_train.col(g),1,n_g);
      arma::mat diff2 = pow((Y_g - fake_FG),2);
      
      for( int small_t=0; small_t<t; small_t++){
        double B = accu( diff2.row(small_t) );
        double bpost  = .5  * B   + b_priorG(small_t,g);
        double apost  = n_g * 1/2 + a_priorG(small_t,g);
        sigma_train_new(small_t,g) = 1/ rgamma(1,apost,1/(bpost))[0];
      }
      
    }else if( n_g==1){
      
      
      arma::colvec diff2 = pow( (Y_g - effe_train.col(g)) ,2 ) ;
      
      //  Rcout << "2.1";
      
      for( int small_t=0; small_t<t; small_t++){
      
      
      double B = accu( diff2.row(small_t) );
        double bpost  = .5  * B   + b_priorG(small_t,g);
        double apost  =  1/2 + a_priorG(small_t,g);
        sigma_train_new(small_t,g) = 1/ rgamma(1,apost,1/(bpost))[0];
        
      }
      //   Rcout << "2";
      
    }else{
      for( int small_t=0; small_t<t; small_t++){
        sigma_train_new(small_t,g) = 1/ rgamma(1, a_priorG(small_t,g),1/( b_priorG(small_t,g)))[0];
      }
      //   Rcout << "3";
    }
  }
  return(sigma_train_new);
}





// [[Rcpp::export]] 
arma::colvec Update_tau_l(double a_tau, double b_tau,
                         arma::mat effe_test, int L, int t, 
                         arma::colvec OOR, arma::colvec a_l){
  arma::colvec new_tau_l(L);
  double apost = a_tau + t*.5;
  int oldL = effe_test.n_cols;
  
  if(oldL>L){
    
    for(int l=0;l<L;l++){
      
      double bpost = 1/(accu( pow( (effe_test.col(l)-a_l[l]),2 ) % OOR )*.5 + b_tau);   
      new_tau_l[l] =      1/ rgamma(1,apost, 
                                   bpost)[0];
    }
  }else{
    for(int l=0;l<oldL;l++){
      
      double bpost = 1/( accu( pow( (effe_test.col(l)-a_l[l]),2 ) % OOR )*.5 + b_tau);   
      new_tau_l[l] =      1/ rgamma(1,apost, 
                                   bpost)[0];
    }
    for(int l=oldL;l<L;l++){
      new_tau_l[l] =      1/ rgamma(1,a_tau, 
                                   b_tau)[0];
    }
  }
  return(new_tau_l);
}

// [[Rcpp::export]] 
arma::colvec Update_a_l(arma::mat effe_test, int L,
                       arma::colvec tau_l,
                       arma::colvec R, double s, int t, 
                       arma::colvec OneOverR){
  
  arma::colvec a_l_new(L);
  int oldL = effe_test.n_cols;
  if(oldL>L){
    
    for(int l=0;l<L;l++){
      double SumY        = accu(effe_test.col(l) % OneOverR);
      //Rcout << SumY;
      double k_          = accu(OneOverR)/(tau_l[l]) + 1/(s);
      //Rcout << k_;
      double m_          = ( SumY / (tau_l[l]) )/k_;
      a_l_new[l]   =   m_ + (rnorm(1))[0] * pow(1/k_ ,.5); 
    }
  }else{
    for(int l=0;l<oldL;l++){
      double SumY        = accu(effe_test.col(l)%OneOverR);
      //Rcout << SumY;
      double k_          = accu(OneOverR)/(tau_l[l]) + 1/(s);
      //Rcout << k_;
      double m_          = ( SumY/(tau_l[l]) )/k_;
      a_l_new[l]   =   m_ + (rnorm(1))[0] * pow(1/k_ ,.5); 
    }
    for(int l=oldL;l<L;l++){
      a_l_new[l]   =   0 + (rnorm(1))[0] * pow(s,.5); 
    }
    
  }
  return(a_l_new);
}











// [[Rcpp::export]]
arma::mat Upd_ZETA_t_SLICE(arma::mat Y, 
                                    arma::mat effe_ALL,
                                    arma::mat sigma_ALL,
                                    arma::colvec log_omegatilde,
                                    arma::colvec Uitilde,
                                    arma::colvec xitilde,
                                    int J, int n,  int L, 
                                    arma::colvec poss_lab){
  int toto = L+J;
  arma::colvec res_Z(n);  res_Z.fill(0);
  arma::colvec likel(toto);
  
  for(int ii=0; ii<n; ii++) {

    for(int l=0; l<toto; l++){
      if(Uitilde[ii] < xitilde[l]){
        likel(l) =  arma::accu(log(normpdf(Y.col(ii), effe_ALL.col(l), sqrt(sigma_ALL.col(l))))) - log(xitilde[l]);  
      }else{
        likel(l) =  -arma::datum::inf;  
      }
    }
    
    likel = log_omegatilde + likel;
    
    arma::colvec r = exp( likel - LogSumExp(likel) );
    res_Z[ii] = RcppArmadillo::sample(poss_lab, 1, TRUE, r)[0];
    
  }
  
  return(res_Z);
}
