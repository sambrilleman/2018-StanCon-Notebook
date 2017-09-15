// Note that this file contains extracts from the rstanarm package.
// Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
// Copyright (C) 2016, 2017 Sam Brilleman
  
functions {
  
  /** 
  * Create group-specific block-diagonal Cholesky factor, see section 2 of
  * https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
  * @param len_theta_L An integer indicating the length of returned vector, 
  *   which lme4 denotes as m
  * @param p An integer array with the number variables on the LHS of each |
  * @param dispersion Scalar standard deviation of the errors, calles sigma by lme4
  * @param tau Vector of scale parameters whose squares are proportional to the 
  *   traces of the relative covariance matrices of the group-specific terms
  * @param scale Vector of prior scales that are multiplied by elements of tau
  * @param zeta Vector of positive parameters that are normalized into simplexes
  *   and multiplied by the trace of the covariance matrix to produce variances
  * @param rho Vector of radii in the onion method for creating Cholesky factors
  * @param z_T Vector used in the onion method for creating Cholesky factors
  * @return A vector that corresponds to theta in lme4
  */
    vector make_theta_L(int len_theta_L, int[] p, real dispersion,
                        vector tau, vector scale, vector zeta,
                        vector rho, vector z_T) {
      vector[len_theta_L] theta_L;
      int zeta_mark = 1;
      int rho_mark = 1;
      int z_T_mark = 1;
      int theta_L_mark = 1;
      
      // each of these is a diagonal block of the implicit Cholesky factor
      for (i in 1:size(p)) { 
        int nc = p[i];
        if (nc == 1) { // "block" is just a standard deviation
          theta_L[theta_L_mark] = tau[i] * scale[i] * dispersion;
          // unlike lme4, theta[theta_L_mark] includes the dispersion term in it
          theta_L_mark = theta_L_mark + 1;
        }
        else { // block is lower-triangular               
          matrix[nc,nc] T_i; 
          real std_dev;
          real T21;
          real trace_T_i = square(tau[i] * scale[i] * dispersion) * nc;
          vector[nc] pi = segment(zeta, zeta_mark, nc); // gamma(zeta | shape, 1)
          pi = pi / sum(pi);                            // thus dirichlet(pi | shape)
          
          // unlike lme4, T_i includes the dispersion term in it
          zeta_mark = zeta_mark + nc;
          std_dev = sqrt(pi[1] * trace_T_i);
          T_i[1,1] = std_dev;
          
          // Put a correlation into T_i[2,1] and scale by std_dev
          std_dev = sqrt(pi[2] * trace_T_i);
          T21 = 2.0 * rho[rho_mark] - 1.0;
          rho_mark = rho_mark + 1;
          T_i[2,2] = std_dev * sqrt(1.0 - square(T21));
          T_i[2,1] = std_dev * T21;
          
          for (r in 2:(nc - 1)) { // scaled onion method to fill T_i
            int rp1 = r + 1;
            vector[r] T_row = segment(z_T, z_T_mark, r);
            real scale_factor = sqrt(rho[rho_mark] / dot_self(T_row)) * std_dev;
            z_T_mark = z_T_mark + r;
            std_dev = sqrt(pi[rp1] * trace_T_i);
            for(c in 1:r) T_i[rp1,c] = T_row[c] * scale_factor;
            T_i[rp1,rp1] = sqrt(1.0 - rho[rho_mark]) * std_dev;
            rho_mark = rho_mark + 1;
          }
          
          // now vech T_i
          for (c in 1:nc) for (r in c:nc) {
            theta_L[theta_L_mark] = T_i[r,c];
            theta_L_mark = theta_L_mark + 1;
          }
        }
      }
      return theta_L;
    }
  
  /** 
  * Create group-specific coefficients, see section 2 of
  * https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
  *
  * @param z_b Vector whose elements are iid normal(0,sigma) a priori
  * @param theta Vector with covariance parameters as defined in lme4
  * @param p An integer array with the number variables on the LHS of each |
  * @param l An integer array with the number of levels for the factor(s) on 
  *   the RHS of each |
  * @return A vector of group-specific coefficients
  */
    vector make_b(vector z_b, vector theta_L, int[] p, int[] l) {
      vector[rows(z_b)] b;
      int b_mark = 1;
      int theta_L_mark = 1;
      for (i in 1:size(p)) {
        int nc = p[i];
        if (nc == 1) {
          real theta_L_start = theta_L[theta_L_mark];
          for (s in b_mark:(b_mark + l[i] - 1)) 
            b[s] = theta_L_start * z_b[s];
          b_mark = b_mark + l[i];
          theta_L_mark = theta_L_mark + 1;
        }
        else {
          matrix[nc,nc] T_i = rep_matrix(0, nc, nc);
          for (c in 1:nc) {
            T_i[c,c] = theta_L[theta_L_mark];
            theta_L_mark = theta_L_mark + 1;
            for(r in (c+1):nc) {
              T_i[r,c] = theta_L[theta_L_mark];
              theta_L_mark = theta_L_mark + 1;
            }
          }
          for (j in 1:l[i]) {
            vector[nc] temp = T_i * segment(z_b, b_mark, nc);
            b_mark = b_mark - 1;
            for (s in 1:nc) b[b_mark + s] = temp[s];
            b_mark = b_mark + nc + 1;
          }
        }
      }
      return b;
    }
  
  /** 
  * Prior on group-specific parameters
  *
  * @param z_b A vector of primitive coefficients
  * @param z_T A vector of primitives for the unit vectors in the onion method
  * @param rho A vector radii for the onion method
  * @param zeta A vector of primitives for the simplexes
  * @param tau A vector of scale parameters
  * @param regularization A real array of LKJ hyperparameters
  * @param delta A real array of concentration paramters
  * @param shape A vector of shape parameters
  * @param t An integer indicating the number of group-specific terms
  * @param p An integer array with the number variables on the LHS of each |
  * @return nothing
  */
    void decov_lp(vector z_b, vector z_T, vector rho, vector zeta, vector tau,
                  real[] regularization, real[] delta, vector shape,
                  int t, int[] p) {
      int pos_reg = 1;
      int pos_rho = 1;
      target += normal_lpdf(z_b | 0, 1);
      target += normal_lpdf(z_T | 0, 1);
      for (i in 1:t) if (p[i] > 1) {
        vector[p[i] - 1] shape1;
        vector[p[i] - 1] shape2;
        real nu = regularization[pos_reg] + 0.5 * (p[i] - 2);
        pos_reg = pos_reg + 1;
        shape1[1] = nu;
        shape2[1] = nu;
        for (j in 2:(p[i]-1)) {
          nu = nu - 0.5;
          shape1[j] = 0.5 * j;
          shape2[j] = nu;
        }
        target += beta_lpdf(rho[pos_rho:(pos_rho + p[i] - 2)] | shape1, shape2);
        pos_rho = pos_rho + p[i] - 1;
      }
      target += gamma_lpdf(zeta | delta, 1);
      target += gamma_lpdf(tau  | shape, 1);
    }
}
data {
  
  //----- Longitudinal submodels

  // NB these design matrices are evaluated AT the observation times.
  // Also, the response vector, design matrix, predictor meanns, etc
  // are combined across all M longitudinal submodels.
  int<lower=0> N;     // num. of obs. across all long. submodels
  int<lower=1> M;     // num. of long. submodels
  int<lower=0> NM[M]; // num. of obs. in each long. submodel
  int<lower=0,upper=N> idx[M,2]; // indices of first and last obs. for each submodel
  int<lower=0> y_K;   // number of predictors across all long. submodels
  int<lower=0> KM[M]; // num. of predictors in each long. submodel
  int<lower=0,upper=y_K> idx_K[M,2]; // index of first/last beta for each submodel
  vector[N] y;        // outcome vector
  vector[y_K] xbar;   // predictor means
  matrix[N,y_K] X;    // centered predictor matrix in the dense case

  // glmer stuff, see table 3 of
  // https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
  int<lower=0> t;               // num. terms (maybe 0) with a | in the glmer formula
  int<lower=1> p[t];            // num. variables on the LHS of each |
  int<lower=1> l[t];            // num. levels for the factor(s) on the RHS of each |
  int<lower=0> q;               // conceptually equals \sum_{i=1}^t p_i \times l_i
  int<lower=0> len_theta_L;     // length of the theta_L vector
  
  // design matrices for group-specific terms
  int<lower=0> num_non_zero;         // number of non-zero elements in the Z matrix
  vector[num_non_zero] w;            // non-zero elements in the implicit Z matrix
  int<lower=0> v[num_non_zero];      // column indices for w
  int<lower=0> u[t > 0 ? N + 1 : 0]; // where the non-zeros start in each row

  //----- Event submodel
  
  // data for calculating event submodel linear predictor in GK quadrature
  // NB these design matrices are evaluated AT the event time and
  // the (unstandardised) quadrature points
  int<lower=0> e_K;           // num. of predictors in event submodel
  int<lower=0> a_K;           // num. of association parameters
  int<lower=0> Npat;          // num. individuals
  int<lower=0> quadnodes;     // num. of nodes for Gauss-Kronrod quadrature 
  int<lower=0> Npat_times_quadnodes; 
  int<lower=0> nrow_e_Xq;     // num. rows in event submodel predictor matrix
  vector[nrow_e_Xq] e_times;  // event times and quadrature points
  matrix[nrow_e_Xq,e_K] e_Xq; // predictor matrix (event submodel)
  vector[e_K] e_xbar;         // predictor means (event submodel)
  int<lower=0> basehaz_df;    // df for B-splines baseline hazard
  matrix[nrow_e_Xq,basehaz_df] 
    basehaz_X; // design matrix (basis terms) for baseline hazard
  vector[nrow_e_Xq] 
    e_d; // event indicator followed by dummy indicators for quadpoints
  vector[Npat_times_quadnodes] 
    quadweight; // GK quadrature weights with (b-a)/2 scaling 

  // data for calculating long. submodel linear predictor in GK quadrature
  // NB these design matrices are evaluated AT the event time and
  // the (unstandardised) quadrature points
  int<lower=0> nrow_y_Xq[M]; // num. rows in long. predictor matrix
  int<lower=0> idx_q[M,2];   // indices of first and last rows in eta at quadpoints
  matrix[sum(nrow_y_Xq),y_K] y_Xq; // predictor matrix (long submodel) at quadpoints, centred     
  int<lower=0> nnz_Zq;        // number of non-zero elements in the Z matrix (at quadpoints)
  vector[nnz_Zq] w_Zq;    // non-zero elements in the implicit Z matrix (at quadpoints)
  int<lower=0> v_Zq[nnz_Zq]; // column indices for w (at quadpoints)
  int<lower=0> u_Zq[(sum(nrow_y_Xq)+1)]; // where the non-zeros start in each row (at quadpoints)

  //----- Hyperparameters for prior distributions
  
  // means for priors
  vector[M]            y_prior_mean_for_intercept;
  vector[y_K]          y_prior_mean;
  vector[e_K]          e_prior_mean;
  vector[a_K]          a_prior_mean;
  vector<lower=0>[M]   y_prior_mean_for_aux;
  vector[basehaz_df]   e_prior_mean_for_aux;
  
  // scale for priors
  vector<lower=0>[M]   y_prior_scale_for_intercept;
  vector<lower=0>[y_K] y_prior_scale;
  vector<lower=0>[e_K] e_prior_scale;
  vector<lower=0>[a_K] a_prior_scale;
  vector<lower=0>[M]   y_prior_scale_for_aux;
  vector<lower=0>[basehaz_df] e_prior_scale_for_aux;
  
  // hyperparameters for glmer stuff
  vector<lower=0>[t] shape; 
  vector<lower=0>[t] scale;
  int<lower=0> len_concentration;
  real<lower=0> concentration[len_concentration];
  int<lower=0> len_regularization;
  real<lower=0> regularization[len_regularization];
}
transformed data {
  int<lower=0> len_z_T = 0;
  int<lower=0> len_var_group = sum(p) * (t > 0);
  int<lower=0> len_rho = sum(p) - t;
  int<lower=1> pos = 1;
  real<lower=0> delta[len_concentration];  

  len_z_T = 0;
  len_var_group = sum(p) * (t > 0);
  len_rho = sum(p) - t;
  pos = 1;
  for (i in 1:t) {
    if (p[i] > 1) {
      for (j in 1:p[i]) {
        delta[pos] = concentration[j];
        pos = pos + 1;
      }
    }
    for (j in 3:p[i]) len_z_T = len_z_T + p[i] - 1;
  }
}
parameters {
  vector[M] gamma;      // intercepts in long. submodels
  vector[y_K] y_z_beta; // primitive coefs in long. submodels
  vector[e_K] e_z_beta; // primitive coefs in event submodel (log hazard ratios)
  vector[a_K] a_z_beta; // primitive assoc params (log hazard ratios)
  vector<lower=0>[M] y_aux_unscaled; // unscaled residual error SDs 
  vector[basehaz_df] e_aux_unscaled; // unscaled coefs for baseline hazard      

  // decomposition of group-specific parameters
  vector[q] z_b;
  vector[len_z_T] z_T;
  vector<lower=0,upper=1>[len_rho] rho;
  vector<lower=0>[len_concentration] zeta;
  vector<lower=0>[t] tau;
}
transformed parameters {
  vector[y_K] y_beta; // combined coefs vector for all long. submodels
  vector[e_K] e_beta; // coefs in event submodel (log hazard ratios)
  vector[a_K] a_beta; // association params in event submodel (log hazard ratios)
  vector[M] y_aux;  // residual error SDs for long. submodels
  vector[basehaz_df] e_aux;  // weibull shape
  vector[q] b;   // group-specific params (ordered by submodel, 1:M)
  vector[len_theta_L] theta_L;
  
  // coefficients and association parameters 
  y_beta = y_z_beta .* y_prior_scale + y_prior_mean;
  e_beta = e_z_beta .* e_prior_scale + e_prior_mean;
  a_beta = a_z_beta .* a_prior_scale + a_prior_mean;
  
  // auxiliary parameters
  y_aux = y_aux_unscaled .* y_prior_scale_for_aux + y_prior_mean_for_aux;
  e_aux = e_aux_unscaled .* e_prior_scale_for_aux + e_prior_mean_for_aux;
 
  // group-specific parameters
  theta_L = make_theta_L(len_theta_L, p, 1.0, tau, scale, zeta, rho, z_T);
  b = make_b(z_b, theta_L, p, l);
}
model {

  //---- Log-lik for longitudinal submodels

  {
    vector[N] y_eta; // linear predictor for long. submodels 
    
    // evaluate linear predictor for all long. submodels
    y_eta = X * y_beta + csr_matrix_times_vector(N, q, w, v, u, b);
  
    // evaluate linear predictor for just submodel m and accumulate log-lik
    for (m in 1:M) {
      vector[NM[m]] y_eta_m; 
      y_eta_m = y_eta[idx[m,1]:idx[m,2]] + gamma[m];
      target += normal_lpdf(y[idx[m,1]:idx[m,2]] | y_eta_m, y_aux[m]);
    }	  
  }
  
  //----- Log-lik for event submodel (Gauss-Kronrod quadrature)
  
  {
    vector[sum(nrow_y_Xq)] y_eta_q; 
    vector[nrow_e_Xq] e_eta_q; 
    vector[nrow_e_Xq] log_basehaz;
    vector[nrow_e_Xq] ll_haz_q;
    vector[Npat] ll_haz_eventtime;
    vector[Npat_times_quadnodes] ll_surv_eventtime;
    real ll_event;
    
    // Event submodel: linear predictor at event time and quadrature points
    e_eta_q = e_Xq * e_beta;
    
    // Long. submodel: linear predictor at event time and quadrature points
    y_eta_q = y_Xq * y_beta + 
      csr_matrix_times_vector(sum(nrow_y_Xq), q, w_Zq, v_Zq, u_Zq, b);
  
    // Event submodel: add on contribution from association structure to
    // the linear predictor at event time and quadrature points
    for (m in 1:M) {
      vector[nrow_y_Xq[m]] y_eta_q_m;
      y_eta_q_m = y_eta_q[idx_q[m,1]:idx_q[m,2]] + gamma[m];
      e_eta_q = e_eta_q + a_beta[m] * y_eta_q_m;
    }	
    
    // Log baseline hazard at event time and quadrature points
    log_basehaz = basehaz_X * e_aux;	
    
    // Log hazard at event time and quadrature points
    ll_haz_q = e_d .* (log_basehaz + e_eta_q);
    
    // Log hazard contribution to the likelihood
    ll_haz_eventtime = segment(ll_haz_q, 1, Npat);
    
    // Log survival contribution to the likelihood, obtained by summing 
    // over the quadrature points to get the approximate integral
    // (NB quadweight already incorporates the (b-a)/2 scaling such that the
    // integral is evaluated over limits (a,b) rather than (-1,+1))
    ll_surv_eventtime = quadweight .* 
      exp(segment(ll_haz_q, Npat + 1, Npat_times_quadnodes));        
    
    // Log likelihood for event submodel
    ll_event = sum(ll_haz_eventtime) - sum(ll_surv_eventtime);
    target += ll_event;  
  }    
    
  //----- Log-priors
    
    // intercepts
    for (m in 1:M)
      target += normal_lpdf(gamma[m] | 
        y_prior_mean_for_intercept, y_prior_scale_for_intercept);    

    // coefficients    
    target += normal_lpdf(y_z_beta | 0, 1);
    target += normal_lpdf(e_z_beta | 0, 1);
    target += normal_lpdf(a_z_beta | 0, 1);
    
    // auxiliary parameters
    target += normal_lpdf(y_aux_unscaled | 0, 1);
    target += normal_lpdf(e_aux_unscaled | 0, 1);

    // group-specific parameters
    decov_lp(z_b, z_T, rho, zeta, tau, regularization, delta, shape, t, p);
}
generated quantities {
  real alpha[M]; // intercept for long. submodels
  
  // evaluate intercept after correcting for centered predictors
  for (m in 1:M) {
    int K1 = idx_K[m,1]; 
    int K2 = idx_K[m,2]; 
    alpha[m] = gamma[m] - dot_product(xbar[K1:K2], y_beta[K1:K2]);
  }
}
