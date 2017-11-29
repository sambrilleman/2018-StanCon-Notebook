// Note that this file contains extracts from the rstanarm package.
// Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
// Copyright (C) 2016, 2017 Sam Brilleman
  
functions {

  /** 
  * Evaluate the linear predictor for the glmer submodel
  *
  * @param X Design matrix for fe
  * @param Z Design matrix for re, for a single grouping factor
  * @param Z_id Group indexing for Z
  * @param gamma The intercept parameter
  * @param beta Vector of population level parameters
  * @param bMat Matrix of group level params
  * @param shift Number of columns in bMat
  *   that correpond to group level params from prior glmer submodels
  * @return A vector containing the linear predictor for the glmer submodel
  */  
  vector evaluate_eta(matrix X, vector[] Z, int[] Z_id, real[] gamma, 
	                    vector beta, matrix bMat, int shift) {
    int N = rows(X);    // num rows in design matrix
    int K = rows(beta); // num predictors
    int p = size(Z);    // num group level params
    vector[N] eta;
    
    if (K > 0) eta = X * beta;
    else eta = rep_vector(0.0, N);
    
    for (k in 1:p)
      for (n in 1:N)
        eta[n] = eta[n] + (bMat[Z_id[n], k + shift]) * Z[k,n];
    
    return eta;
  }	
	
}
data {
  
  //----- Longitudinal submodels

	// number of long. submodels
	// NB this is fixed equal to 2 for this simplified jm.stan 
	// file. See the jm.stan file in rstanarm for more general 
	// code that allows for 1, 2 or 3 longitudinal submodels.
  int<lower=2,upper=2> M; 
	
  // population level dimensions
  int<lower=0> yNobs[2]; // num observations
  int<lower=0> yNeta[2]; // required length of eta
  int<lower=0> yK[2]; // num predictors

	// population level data
  // NB these design matrices are evaluated AT the observation times.
  vector[yNobs[1]] yReal1; // response vectors
  vector[yNobs[2]] yReal2;  
  matrix[yNeta[1],yK[1]] yX1; // fe design matrix
  matrix[yNeta[2],yK[2]] yX2; 
  vector[yK[1]] yXbar1; // predictor means
  vector[yK[2]] yXbar2;
	
  // group level dimensions
	int<lower=0> bN1; // num groups
	int<lower=0> bK1; // total num params
	int<lower=0> bK1_len[2]; // num params in each submodel
	int<lower=0> bK1_idx[2,2]; // beg/end index for group params

  // group level data
  vector[bK1_len[1] > 0 ? yNeta[1] : 0] y1_Z1[bK1_len[1]]; // re design matrix
  vector[bK1_len[2] > 0 ? yNeta[2] : 0] y2_Z1[bK1_len[2]];
  int<lower=0> y1_Z1_id[bK1_len[1] > 0 ? yNeta[1] : 0]; // group indexing for y1_Z1
  int<lower=0> y2_Z1_id[bK1_len[2] > 0 ? yNeta[2] : 0]; // group indexing for y2_Z1
  
  //----- Event submodel
  
  // data for calculating event submodel linear predictor in GK quadrature
  // NB these design matrices are evaluated AT the event time and
  // the (unstandardised) quadrature points
  real norm_const;            // constant shift for log baseline hazard
  int<lower=0> e_K;           // num. of predictors in event submodel
  int<lower=0> a_K;           // num. of association parameters
  int<lower=0> Npat;          // num. individuals
  int<lower=0> Nevents;       // num. events (ie. not censored)  
	int<lower=0> qnodes;        // num. of nodes for GK quadrature 
  int<lower=0> Npat_times_qnodes; 
  int<lower=0> nrow_e_Xq;     // num. rows in event submodel predictor matrix
  vector[nrow_e_Xq] e_times;  // event times and quadrature points
  matrix[nrow_e_Xq,e_K] e_Xq; // predictor matrix (event submodel)
  vector[e_K] e_xbar;         // predictor means (event submodel)
  int<lower=0> basehaz_df;    // df for B-splines baseline hazard
  matrix[nrow_e_Xq,basehaz_df] basehaz_X; // design matrix (basis terms) for baseline hazard
  vector[Npat_times_qnodes] qwts; // GK quadrature weights with (b-a)/2 scaling 

  // data for calculating long. submodel linear predictor in GK quadrature
  // NB these design matrices are evaluated AT the event time and
  // the (unstandardised) quadrature points
  int<lower=0> nrow_y_Xq[2]; // num. rows in long. predictor matrix at quadpoints
  matrix[nrow_y_Xq[1],yK[1]] y1_xq_eta; // fe design matrix at quadpoints
  matrix[nrow_y_Xq[2],yK[2]] y2_xq_eta; 
  vector[nrow_y_Xq[1]] y1_z1q_eta[bK1_len[1]]; // re design matrix at quadpoints
  vector[nrow_y_Xq[2]] y2_z1q_eta[bK1_len[2]];
  int<lower=0> y1_z1q_id_eta[nrow_y_Xq[1]]; // group indexing for re design matrix
  int<lower=0> y2_z1q_id_eta[nrow_y_Xq[2]]; 

  //----- Hyperparameters for prior distributions
  
  // means for priors
  // coefficients
  vector[M]            y_prior_mean_for_intercept;
  vector[yK[1]]        y_prior_mean1;
  vector[yK[2]]        y_prior_mean2;
  vector[e_K]          e_prior_mean;
  vector[a_K]          a_prior_mean;
  vector<lower=0>[M]   y_prior_mean_for_aux;
  vector[basehaz_df]   e_prior_mean_for_aux;
  
  // scale for priors
  vector<lower=0>[M]     y_prior_scale_for_intercept;
  vector<lower=0>[yK[1]] y_prior_scale1;
  vector<lower=0>[yK[2]] y_prior_scale2;
  vector<lower=0>[e_K]   e_prior_scale;
  vector<lower=0>[a_K]   a_prior_scale;
  vector<lower=0>[M]     y_prior_scale_for_aux;
  vector<lower=0>[basehaz_df] e_prior_scale_for_aux;
  
  // lkj prior stuff
  vector<lower=0>[bK1] b1_prior_scale;
  vector<lower=0>[bK1] b1_prior_df;
  real<lower=0> b1_prior_regularization;
}
transformed data {
  // indexing used to extract lower tri of RE covariance matrix
  int bCov1_idx[bK1 + choose(bK1, 2)];
  if (bK1 > 0) 
    bCov1_idx = lower_tri_indices(bK1);
}
parameters {
  real yGamma1;                 // intercepts in long. submodels
	real yGamma2;
  vector[yK[1]] z_yBeta1;       // primitive coefs in long. submodels
	vector[yK[2]] z_yBeta2;
  real<lower=0> yAux1_unscaled; // unscaled residual error SDs 
  real<lower=0> yAux2_unscaled; 
  vector[e_K] e_z_beta; // primitive coefs in event submodel (log hazard ratios)
  vector[a_K] a_z_beta; // primitive assoc params (log hazard ratios)
  vector[basehaz_df] e_aux_unscaled; // unscaled coefs for baseline hazard      

  // group level params   
  vector<lower=0>[bK1] bSd1; // group level sds  
  matrix[bK1,bN1] z_bMat1;   // unscaled group level params 
  cholesky_factor_corr[bK1 > 1 ? bK1 : 0] 
	  bCholesky1;              // cholesky factor of corr matrix
}
transformed parameters {
	vector[yK[1]] yBeta1;     // population level params for long. submodels
  vector[yK[2]] yBeta2;
	real yAux1[has_aux[1]];   // residual error SDs for long. submodels
  real yAux2[has_aux[2]];  
	vector[e_K] e_beta;       // coefs in event submodel (log hazard ratios)
  vector[a_K] a_beta;       // assoc params in event submodel (log hazard ratios) 
	vector[basehaz_df] e_aux; // b-spline coefs for baseline hazard
	matrix[bN1,bK1] bMat1;    // group level params
  
  // coefs for long. submodels
	yBeta1 = z_yBeta1 .* y_prior_scale1 + y_prior_mean1;
	yBeta2 = z_yBeta2 .* y_prior_scale2 + y_prior_mean2;

	// coefs for event submodel (incl. association parameters)
	e_beta = e_z_beta .* e_prior_scale + e_prior_mean;
  a_beta = a_z_beta .* a_prior_scale + a_prior_mean;
  
  // residual error SDs for long. submodels
	yAux1[1] = yAux1_unscaled[1] * y_prior_scale_for_aux[1] + y_prior_mean_for_aux[1];
	yAux2[1] = yAux2_unscaled[1] * y_prior_scale_for_aux[2] + y_prior_mean_for_aux[2];

	// b-spline coefs for baseline hazard
	e_aux = e_aux_unscaled .* e_prior_scale_for_aux + e_prior_mean_for_aux;
 
  // group level params
	if (bK1 == 1) 
		bMat1 = (bSd1[1] * z_bMat1)'; 
	else if (bK1 > 1) 
		bMat1 = (diag_pre_multiply(bSd1, bCholesky1) * z_bMat1)';
}
model {

  //---- Log-lik for longitudinal submodels

  {
	  // declare linear predictors
		vector[yNeta[1]] yEta1; 
		vector[yNeta[2]] yEta2;

    // evaluate linear predictor for each long. submodel
		yEta1 = evaluate_eta(yX1, y1_Z1, y1_Z1_id, yGamma1, yBeta1, bMat1, 0);
    yEta2 = evaluate_eta(yX2, y2_Z1, y2_Z1_id, yGamma2, yBeta2, bMat1, bK1_len[1]);
		
    // increment the target with the log-lik
    target += normal_lpdf(yReal1 | yEta1, yAux1[1]);
    target += normal_lpdf(yReal2 | yEta2, yAux2[1]);
  }
  
  //----- Log-lik for event submodel (Gauss-Kronrod quadrature)
  
  {
		vector[nrow_y_Xq[1]] yEta1_q; 
		vector[nrow_y_Xq[2]] yEta2_q;
    vector[nrow_e_Xq] e_eta_q; 
		vector[nrow_e_Xq] log_basehaz;  // log baseline hazard AT event time and quadrature points
		vector[nrow_e_Xq] log_haz_q;    // log hazard AT event time and quadrature points
		vector[Nevents] log_haz_etimes; // log hazard AT the event time only
		vector[Npat_times_qnodes] log_haz_qtimes; // log hazard AT the quadrature points
    
    // Event submodel: linear predictor at event time and quadrature points
    e_eta_q = e_Xq * e_beta;
    
    // Long. submodel: linear predictor at event time and quadrature points
		yEta1_q = evaluate_eta(y1_xq_eta, y1_z1q_eta, y1_z1q_id_eta, yGamma1, yBeta1, bMat1, 0);
    yEta2_q = evaluate_eta(y2_xq_eta, y2_z1q_eta, y2_z1q_id_eta, yGamma2, yBeta2, bMat1, bK1_len[1]);
  
    // Event submodel: add on contribution from association structure to
    // the linear predictor at event time and quadrature points
    e_eta_q = e_eta_q + a_beta[1] * yEta1_q + a_beta[2] * yEta1_q;
    
    // Log baseline hazard at event time and quadrature points
    log_basehaz = norm_const + basehaz_X * e_aux;	
    
    // Log hazard at event time and quadrature points
    log_haz_q = log_basehaz + e_eta_q;
  
	  // Log hazard at event times only
	  log_haz_etimes = head(log_haz_q, Nevents);
  
		// Log hazard at quadrature points only
		log_haz_qtimes = tail(log_haz_q, Npat_times_qnodes);
   
    // Log survival contribution to the likelihood, obtained by summing 
    // over the quadrature points to get the approximate integral
    // (NB quadweight already incorporates the (b-a)/2 scaling such that the
    // integral is evaluated over limits (a,b) rather than (-1,+1))
    ll_surv_eventtime = quadweight .* 
      exp(segment(ll_haz_q, Npat + 1, Npat_times_qnodes));        
    
    // Log likelihood for event submodel
		// NB The first term is the log hazard contribution to the log  
		// likelihood for the event submodel. The second term is the log  
		// survival contribution to the log likelihood for the event submodel.  
		// The latter is obtained by summing over the quadrature points to get 
    // the approximate integral (i.e. cumulative hazard). Note that the
    // 'qwts' vector already incorporates (b-a)/2 scaling such that the
    // integral is evaluated over limits (a,b) rather than (-1,+1), where
		// 'a' is baseline, i.e. time 0, and 'b' is the event or censoring
    // time for the individual.
    target += sum(log_haz_etimes) - dot_product(qwts, exp(log_haz_qtimes));  
  }    
    
  //----- Log-priors
    
    // intercepts for long. submodels
    target += normal_lpdf(yGamma1[1] | 
      y_prior_mean_for_intercept[1], y_prior_scale_for_intercept[1]);    
    target += normal_lpdf(yGamma2[1] | 
      y_prior_mean_for_intercept[2], y_prior_scale_for_intercept[2]);    

    // coefficients for long. submodels   
    target += normal_lpdf(z_yBeta1 | 0, 1);
    target += normal_lpdf(z_yBeta2 | 0, 1);
		
		// coefficients for event submodel
    target += normal_lpdf(e_z_beta | 0, 1);
    target += normal_lpdf(a_z_beta | 0, 1);
    
    // residual error SDs for long. submodels
    target += normal_lpdf(yAux1_unscaled[1] | 0, 1);
    target += normal_lpdf(yAux2_unscaled[1] | 0, 1);
		
		// b-spline coefs for baseline hazard
    target += normal_lpdf(e_aux_unscaled | 0, 1);

    // group level terms
      // sds
      target += student_t_lpdf(bSd1 | b1_prior_df, 0, b1_prior_scale);
      // primitive coefs
      target += normal_lpdf(to_vector(z_bMat1) | 0, 1); 
      // corr matrix
      if (bK1 > 1) 
        target += lkj_corr_cholesky_lpdf(bCholesky1 | b1_prior_regularization);
}
generated quantities {
	real yAlpha1; // transformed intercepts for long. submodels
  real yAlpha2;
	real e_alpha; // transformed intercept for event submodel 
  vector[size(bCov1_idx)] bCov1; // var-cov for REs
    
  // Transformed intercepts for long. submodels
	yAlpha1[1] = yGamma1[1] - dot_product(yXbar1, yBeta1);
  yAlpha2[1] = yGamma2[1] - dot_product(yXbar2, yBeta2);
	
  // Transformed intercept for event submodel 
	// NB norm_const is a constant shift in log baseline hazard,
	// used so that the log baseline hazard is equal to the mean
	// log incidence rate, and not equal to zero, when all 
	// coefficients in the event submodel (i.e. the log hazard 
	// ratios and the baseline hazard b-spline coefs) are set 
	// equal to zero.
  e_alpha = norm_const - dot_product(e_xbar, e_beta);
	
	// Transform variance-covariance matrix for REs
	if (bK1 == 1)
    bCov1[1] = bSd1[1] * bSd1[1];
  else
  	bCov1 = to_vector(quad_form_diag(
  	  multiply_lower_tri_self_transpose(bCholesky1), bSd1))[bCov1_idx];
}
