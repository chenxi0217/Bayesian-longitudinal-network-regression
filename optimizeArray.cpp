#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector optimizeArray(int n, int num_visit_maximum, int num_region, int num_subnetwork,
                            NumericMatrix c_latent_sim, NumericVector mit_qr, double ns) {
  
  // Initialize a flattened 4D array (as a 1D NumericVector)
  NumericVector A(n * num_visit_maximum * num_region * num_region);
  
  // Helper function to access elements in a flattened 4D array stored in a 1D vector
  auto array_index = [&](int u, int r, int p, int q) {
    return u + n * (r + num_visit_maximum * (p + num_region * q));
  };
  
  // Loop over the regions, visits, and subnetworks
  for (int p = 0; p < num_region; ++p) {
    for (int q = 0; q < num_region; ++q) {
      for (int u = 0; u < n; ++u) {
        for (int r = 0; r < num_visit_maximum; ++r) {
          for (int s = 0; s < num_subnetwork; ++s) {
            for (int t = 0; t < num_subnetwork; ++t) {
              
              // Only compute if conditions are met
              if (c_latent_sim(p, s) == 1 && c_latent_sim(q, t) == 1) {
                int min_l_p = std::min(s + 1, t + 1);  // Adjust for 1-based indexing
                int max_l_p = std::max(s + 1, t + 1);
                
                int cal = (min_l_p - 1) * num_subnetwork + max_l_p - (min_l_p * (min_l_p - 1)) / 2;
                
                A[array_index(u, r, p, q)] = R::rnorm(mit_qr[u +r* n+ num_visit_maximum*n * (cal-1)], ns);
              }
            }
          }
        }
      }
    }
  }
  
  return A;
}


// [[Rcpp::export]]
List optimize_est(int n, int num_visit_maximum, int num_region, int num_subnetwork, NumericMatrix g_c, NumericVector g_mik_qr, NumericVector g_delta_f) {
  
  // Initialize a flattened 4D array (as a 1D NumericVector)
  NumericVector A_est(n * num_visit_maximum * num_region * num_region);
  NumericVector delta_est(num_region* num_region *3);

  // Helper function to access elements in a flattened 4D array stored in a 1D vector
  auto array_index = [&](int u, int r, int p, int q) {
    return u + n * (r + num_visit_maximum * (p + num_region * q));
  };
  
  // Loop over the regions, visits, and subnetworks
  for (int p = 0; p < num_region; ++p) {
    for (int q = 0; q < num_region; ++q) {
      for (int u = 0; u < n; ++u) {
        for (int r = 0; r < num_visit_maximum; ++r) {
          for (int s = 0; s < num_subnetwork; ++s) {
            for (int t = 0; t < num_subnetwork; ++t) {
              // Only compute if conditions are met
              if (g_c(p, s) == 1 && g_c(q, t) == 1) {
                int min_l_p = std::min(s+1, t +1);  // Adjust for 1-based indexing
                int max_l_p = std::max(s +1, t +1);
                
                int cal = (min_l_p - 1) * num_subnetwork + max_l_p - (min_l_p * (min_l_p - 1)) / 2;
                
                A_est[array_index(u, r, p, q)] = g_mik_qr[u +r* n+ num_visit_maximum*n * (cal-1)];
               for (int m = 0; m < 3; ++m) {
 		              delta_est[p+q*num_region+num_region*num_region*m]=g_delta_f((cal-1),m);
		            }
              }
            }
          }
        }
      }
    }
  }
  
  return List::create(Named("A_est") = A_est,
                      Named("delta_est") = delta_est);
}


// [[Rcpp::export]]
NumericMatrix c_latent_new( int num_region, int num_subnetwork, NumericMatrix max_l_p,NumericMatrix min_l_p,
                            NumericMatrix k_cal, NumericVector A,NumericVector mik_qr, NumericMatrix sigma2qr_matrix,
                            NumericVector paic, int num_visit_maximum, int n,NumericMatrix c_latent_old) {
  
  // Initialize 
  NumericMatrix  new_c_latent(num_region,num_subnetwork);
  NumericMatrix paic_sp_pos(num_region,num_subnetwork);
  double matrix_sum;
  // Helper function to access elements in a flattened 4D array stored in a 1D vector
  auto array_index = [&](int a, int b, int s, int l) {
    return a + n * (b + num_visit_maximum * (s + num_region * l));
  };
  
    
  for (int s=0; s< num_region;s++){
    for (int p=0;p <num_subnetwork;p++){
      double cloglik=0;
      for (int l=0;l< num_region;l++){
        for(int a=0;a< n;a++){
          for(int b=0;b< num_visit_maximum;b++){
           matrix_sum= R::dnorm(A[array_index(a,b,s,l)],mik_qr[a+b*n+k_cal(p,l)*num_visit_maximum*n],
                                 sqrt(sigma2qr_matrix(min_l_p(p,l),max_l_p(p,l))),true);
            cloglik+= matrix_sum;
          }
         }
      }
      paic_sp_pos(s,p)= cloglik+ c_latent_old(s,p)*log(paic[p]);
    }
    paic_sp_pos.row(s)= paic_sp_pos.row(s)- max(paic_sp_pos.row(s));
    paic_sp_pos.row(s)= exp(paic_sp_pos.row(s))/std::accumulate(exp(paic_sp_pos.row(s)).begin(), exp(paic_sp_pos.row(s)).end(), 0.0);
    
  }
  return  paic_sp_pos;
}



// [[Rcpp::export]]
NumericVector bic(int n, int num_visit_maximum, int num_region, int num_subnetwork, NumericMatrix g_c, NumericVector g_mik_qr, NumericVector A, NumericMatrix g_sigma2qr) {
  
  // Initialize a flattened 4D array (as a 1D NumericVector)
  NumericVector bic(n * num_visit_maximum * num_region * num_region);
  NumericVector dic(n * num_visit_maximum * num_region * num_region);

  // Helper function to access elements in a flattened 4D array stored in a 1D vector
  auto array_index = [&](int u, int k, int p, int q) {
    return u + n * (k + num_visit_maximum * (p + num_region * q));
  };
  
  // Loop over the regions, visits, and subnetworks
  for (int p = 0; p < num_region; ++p) {
    for (int q = 0; q < num_region; ++q) {
      for (int u = 0; u < n; ++u) {
        for (int k = 0; k < num_visit_maximum; ++k) {
          for (int s = 0; s < num_subnetwork; ++s) {
            for (int t = 0; t < num_subnetwork; ++t) {
              // Only compute if conditions are met
              if (g_c(p, s) == 1 && g_c(q, t) == 1) {
                int min_l_p = std::min(s+1, t +1);  // Adjust for 1-based indexing
                int max_l_p = std::max(s +1, t +1);
                
                int cal = (min_l_p - 1) * num_subnetwork + max_l_p - (min_l_p * (min_l_p - 1)) / 2;
                
                bic[array_index(u, k, p, q)] = R::dnorm(A[array_index(u, k, p, q)],g_mik_qr[u +k* n+ num_visit_maximum*n * (cal-1)], sqrt(g_sigma2qr(s,t)),true);
              }
            }
          }
        }
      }
    }
  }
  
  return bic;
  
}
