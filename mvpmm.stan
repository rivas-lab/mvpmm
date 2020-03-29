// Manuel A. Rivas, 6.6.2017
// Stanford University 
// Rivas Lab 
// Model to estimate covariance parameters 

data {
    int<lower=0> N; // number of loci
    int<lower=1> M; // number of phenotypes
    matrix[N, M] B; // observed effect sizes
    matrix[N, M] SE; // standard errors
    int<lower=1> K; // number of mixture components
    
}

transformed data{
    vector[M] zeros;
    vector[K] ones; // mass for dirichlet
    vector[M] SE_vec[N];
    
    zeros = rep_vector(0, M);
    ones = rep_vector(1, K);
    ones[1] = 1;
    // fill in the matrix of standard errors
    for (n in 1:N) {
        SE_vec[n] = to_vector(SE[n]);
    }
}

parameters {
    simplex[K] pi; // mixing proportions
    cholesky_factor_corr[M] L_Omega;
    cholesky_factor_corr[M] L_Theta;
    vector<lower=0>[M] tau;
}

transformed parameters{
    matrix[M, M] Sigma;
    matrix[M, M] Sigmas[K];
    Sigma = diag_pre_multiply(tau, L_Omega);
    Sigmas[1] = diag_matrix(rep_vector(0,M));	
    Sigmas[2] = Sigma;
}

model {
    vector[K] ps; // contributions of each
    tau ~ cauchy(0, 2.5);
    L_Omega ~ lkj_corr_cholesky(2.0);
    L_Theta ~ lkj_corr_cholesky(2.0);
    pi ~ dirichlet(ones);
    for (n in 1:N){
       // two components
       for (k in 1:K){
           ps[k] = log(pi[k]) + multi_normal_cholesky_lpdf(B[n] | zeros,  Sigmas[k] + diag_pre_multiply(SE_vec[n], L_Theta));
       }
       target += log_sum_exp(ps);
     }
 
}

generated quantities {
    matrix[M,M] Omegacor;
    matrix[M,M] Thetacor;
    // matrix[N,K] p;
    Omegacor = multiply_lower_tri_self_transpose(L_Omega);
    Thetacor = multiply_lower_tri_self_transpose(L_Theta);
}

