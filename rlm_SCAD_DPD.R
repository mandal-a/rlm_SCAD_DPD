rlm_SCAD_DPD <- function(X=NULL, Y=NULL, data=NULL, Y_index=NULL, X_test=NULL,
                         Y_test=NULL, newdata=NULL, validation_prop=NULL, 
                         trim=0.05, gamma=3.7, SIS_prop = 0.25, SIS_rho.min=NULL, 
                         alpha=seq(0, 0.2, length.out=21), alpha_robust=0.15, 
                         verbose = TRUE, tol=1e-4, maxIter=100, n.lambda=100,
                         seed=NULL){
  
  # Robust linear (regression) model (rlm) fitting using the SCAD-penalized 
  # density power divergent (DPD) estimator for high-dimensional data
  #
  # Linear Regression Model:
  #     Y_i = beta_0 + \sum_{j=1}^p beta_j * X_ij + e_i,   i = 1, 2, ..., n
  
  # Inputs:
  # X: input matrix, of dimension n x p; each row is an observation vector.
  # Y: response variable.
  # data: data frame containing the full data.
  # Y_index: (scalar) the column index or (character) column name for Y in data,
  #          other columns in the data will be in X.
  #          Either (Y, X) or (data, Y_index) must be provided.
  #          (Y, X) will be in (vector, matrix) format 
  #          (data, Y_index) will be in (data frame, scalar, or character) format
  # Y_test, X_test: test data
  # newdata: new test data
  #          Either (Y_test, X_test) or newdata may be provided (optional)
  # validation_prop: proportion of test data for validation method.
  #         If Y_test, X_test, and validation_prop are all NULL, the training mean
  #         absolute deviation (MAD) and trimmed root mean prediction 
  #         error (RMPE)  are calculated.
  #         If newdata is provided, validation_prop is ignored.
  # trim: trimming proportion for the trimmed RMPE 
  # gamma: SCAD parameter
  # SIS_prop: proportion (q) of variables relative to the sample size (n) to 
  #             be used (i.e, n*q) after sure independence screening (SIS) 
  # SIS_rho.min: threshold for correlation coefficient for SIS screening 
  # alpha: (vector) set of DPD parameters
  # alpha_robust: robust alpha for estimating sigma_hat for the RCp (needs tuning)
  # n.lambda: maximum number of lambda values (penalty parameter)
  # verbose: logical value for printing important information
  # tol, maxIter: tolerance limit and maximum number of steps for convergence 
  #               of an MDPDE
  # seed: sets seed if CV is used
  
  # Outputs:
  # A list containing summary and output_beta
  # summary: Summary table for the robust AIC and Cp. The optimum DPD
  #             parameter is selected adaptively using the H-score. 
  # output_beta: A (sparse) matrix containing estimated regression 
  #              coefficient (beta) based on the robust AIC and Cp.
  
  # Simulation Example:
  if (FALSE){
    # Training Data
    library(MASS)
    n = 250; p = 500; sigma = 2; rho = 0.8
    Sigma = rho^abs(outer(1:p, 1:p, "-"))
    X = mvrnorm(n, rep(0, p), Sigma)
    epsilon = rnorm(n, sd = sigma)
    epsilon[sample(n, 12)] = rnorm(12, mean=10) # outliers
    beta = (-1)^sample(1:2, p, replace = TRUE) * runif(p, 1, 2)
    beta[sample(p, p*0.98)] = 0
    beta = c(runif(1, -5, 5), beta) # adds intercept
    Y = beta[1] + X %*% beta[-1] + epsilon
    
    # Penalized DPD Fit
    DPD_fit = rlm_SCAD_DPD(Y=Y, X=X)
    DPD_fit$summary
    
    non_zero_beta = cbind(beta, DPD_fit$beta)
    colnames(non_zero_beta) = c("True beta", "Robust AIC", "Robust Cp")
    (non_zero_beta = non_zero_beta[rowSums(abs(non_zero_beta)) != 0,])
  }
  
  # Version: 11.05.2025
  # Abhijit Mandal, University of Texas at El Paso, USA (Email: amandal@utep.edu)
  # Reference: Mandal, A. and Ghosh, S. Robust Variable Selection Criteria for 
  # the Penalized Regression. Journal of Multivariate Analysis (under review). 
  # ____________________________________________________________________________
  library(MASS) # for ginv
  library(expm) # for bdiag 
  library(Matrix) # for Matrix (sparse representation of beta estimates)
  
  # Global variables for printing important information
  verbose <<- verbose 
  
  if (!is.null(validation_prop)){
    if (is.null(seed)) seed = round(10000*runif(1))
    set.seed(seed); seed <<- seed
    if (verbose)cat(sprintf("Random seed = %d.\n", seed))
  }
  # ____________________________________________________________________________
  #                       Check data
  # ____________________________________________________________________________
  if (!is.null(data) & !is.null(Y_index)){
    if (is.character(Y_index)) Y_index = which(colnames(data) == Y_index)
    Y = as.vector(data[,Y_index])
    X = data.matrix(data[,-Y_index])
    
    if (!is.null(newdata)) {
      data_var_names = colnames(data)
      newdata_var_names = colnames(newdata)
      if (all(newdata_var_names %in% data_var_names)){
        newdata = newdata[,data_var_names]
        if (is.character(Y_index)) Y_index = which(colnames(data) == Y_index)
      } else stop("Variables in data and newdata do not match.")
      
      Y_test = as.vector(newdata[,Y_index])
      X_test = data.matrix(newdata[,-Y_index])
    }
  }
  
  if (is.vector(X)){
    stop("Only one predictor variable: use OLS or RLM.\n")
  } else {X=as.matrix(X); Y = as.vector(Y); n0=nrow(X);  p = ncol(X) }
  if (verbose) cat(sprintf("Sample size: %d, covariates: %d\n", n0, p))
  
  if (length(Y) != n0) stop("Sample size must be same for X and Y.\n")
  # ____________________________________________________________________________
  #                       Test data
  # ____________________________________________________________________________
  if (!is.null(X_test)){
    X_test = as.matrix(X_test); Y_test = as.vector(Y_test)
    if (nrow(X_test) != length(Y_test)) Y_test = NULL
  } else if (is.null(validation_prop)) {
    X_test = X; Y_test = Y
    if (verbose) 
      cat(sprintf("\nCalculating trimmed RMPE and MAD in the training data.\n"))
  } else if (validation_prop < 1) {
    n_test = round(n0 * validation_prop) #test data size
    n = n0 - n_test #training data size
    
    # The design matrix
    test_index = sample(n0, size=n_test)
    X_test = X[test_index, ]
    Y_test = Y[test_index]
    
    X = X[-test_index, ]
    Y = Y[-test_index]
  } else stop("Test data X_test is not valid.")
  
  if (verbose) cat(sprintf("\nTraining data n=%i, test data m=%i.\n", 
                           nrow(X), nrow(X_test)))
  
  X0 = X; X_test0 = X_test #stores original data
  # ____________________________________________________________________________
  if (is.null(SIS_prop) & is.null(SIS_rho.min)){
    if (n < p) {
      stop("n < p: Use eigher non-null SIS_rho.min or SIS_prop for SIS.")
    } else {
      if (verbose) 
        cat(sprintf("\nNo sure independent screening (SIS) is done.\n"))
    }
  }
  
  if (n * SIS_prop >= p){
    SIS_prop = NULL
    cat(sprintf("\nn * SIS_prop >= p: No sure independent screening (SIS) is done.\n"))
  }
  
  n_alpha = length(alpha)
  use_global = TRUE
  
  if (verbose) cat(sprintf("\nNumber of lambda: %d.\n\n", n.lambda))
  # ____________________________________________________________________________
  #                Initialization for output objects
  # ____________________________________________________________________________
  #Information criteria
  IC_type=c("AIC", "Cp")
  n.IC = length(IC_type)
  
  opt_alpha = opt_alpha_lambda =  opt_alpha_RMPE_trim = opt_alpha_MAD = 
    opt_alpha_dim_reduction = opt_alpha_sigma = rep(NA, n.IC)
  names(opt_alpha) = names(opt_alpha_lambda) = names(opt_alpha_RMPE_trim) = 
    names(opt_alpha_MAD) = names(opt_alpha_dim_reduction) = 
    names(opt_alpha_sigma) = IC_type
  
  df_summary = matrix(NA, n_alpha, 5)
  
  #Proportion of DPD failed throughout the lambda path
  DPD_fail = rep(NA, n_alpha)
  names(DPD_fail) = rownames(df_summary) = alpha
  
  # dimension reduction 
  lambda_opt_DPD = dim_reduction_DPD = sigma_DPD = RMPE_trim_DPD = 
    MAD_DPD = H_Score_DPD = matrix(NA, n_alpha, n.IC)
  
  rownames(lambda_opt_DPD) = rownames(dim_reduction_DPD) = 
    rownames(sigma_DPD) = rownames(RMPE_trim_DPD) = rownames(MAD_DPD) = 
    rownames(H_Score_DPD) = alpha
  
  colnames(lambda_opt_DPD) = colnames(dim_reduction_DPD) = 
    colnames(sigma_DPD) = colnames(RMPE_trim_DPD) = colnames(MAD_DPD) = 
    colnames(H_Score_DPD) = IC_type
  
  var_names = colnames(X)
  if (is.null(var_names)) var_names = paste0("X", 1:p)
  beta_DPD = array(0, dim = c(n_alpha, p+1, n.IC), 
                   dimnames = list(alpha, c("intercept", var_names), IC_type))
  
  Create_global(X=X)
  # ____________________________________________________________________________
  #                    MDPDE
  # ____________________________________________________________________________
  RSIS_index = Robust_SIS(X=X, Y=Y, SIS_prop=SIS_prop, 
                          SIS_rho.min=SIS_rho.min, alpha=alpha_robust, 
                          beta_init=NULL, sigma_init=NULL, is.rank=FALSE,
                          tol = tol, maxIter=maxIter, verbose=verbose)
  #update global variables
  if (length(RSIS_index) != p) Create_global(X[, RSIS_index]) 
  
  robust_est = MDPDE_SCAD(X=X[, RSIS_index], Y=Y, alpha=alpha_robust, lambda = 0,
                          use_global=use_global)
  #Robust unbiased estimator of RCp
  sigma_hat_robust_unbiased = robust_est$sigma
  if (!robust_est$converged){
    Huber_fit = rlm(Y ~ X[, RSIS_index], maxit = 100)
    sigma_hat_robust_unbiased = Huber_fit$s
  }
  robust_DPD_converged = robust_est$converged
  if (verbose & !robust_est$converged) sprintf("Robust DPD did not converge.\n")
  
  # ____________________________________________________________________________
  for (j in 1:n_alpha){
    if (verbose) 
      cat(sprintf('  alpha index %i/%i\n', j, n_alpha))
    
    RSIS_index = Robust_SIS(X=X0, Y=Y, SIS_prop=SIS_prop, 
                            SIS_rho.min=SIS_rho.min, alpha=alpha[j], 
                            beta_init=NULL, sigma_init=NULL, is.rank=FALSE,
                            tol = tol, maxIter=maxIter, verbose=verbose)
    if (length(RSIS_index) != p){
      X = X0[, RSIS_index]; X_test = X_test0[, RSIS_index]
      Create_global(X) #update global variables
    }
    
    lambda_vec = SetupLambda(X=X, Y=Y, gamma=gamma,
                             alpha=alpha[j], beta_init=NULL,
                             sigma_init=NULL, n.lambda=n.lambda,
                             lambda.min=NULL, lambda.max=NULL,
                             use_global=use_global,
                             tol=tol, maxIter=maxIter, verbose=verbose)
    
    DPD_fit = Opt_penalty(X=X, Y=Y, alpha=alpha[j], 
                          sigma_unbiased=sigma_hat_robust_unbiased,
                          lambda_vec=lambda_vec, gamma=gamma, 
                          use_global=use_global, IC_type=IC_type,
                          tol=tol, maxIter=maxIter, verbose=verbose)
    
    #Proportion of DPD failed throughout the lambda path
    DPD_fail[j] = mean(!DPD_fit$summary_tab[,"converged"]) 
    df_summary[j,] = DPD_fit$df_summary
    # __________________________________________________________________
    #                       Outputs: DPD
    # __________________________________________________________________
    lambda_opt_DPD[j,] = DPD_fit$opt_lambda
    
    # Trimmed RMPE (for training data, to be done for test data)
    temp = predict_y(X=X, beta_hat_mat=DPD_fit$opt_MDPDE_beta,
                     Y_test=Y_test, X_test=X_test, trim=trim)
    RMPE_trim_DPD[j, ] = temp$RMPE_trim
    MAD_DPD[j, ] = apply(temp$y_hat_mat, 2, function(x) mean(abs(x - Y_test))) 
    
    dim_reduction_DPD[j,] = 100*colMeans(!DPD_fit$X_col_selected)
    
    beta_DPD[j, c(1, 1+RSIS_index), ] = DPD_fit$opt_MDPDE_beta 
    
    sigma_DPD[j, ] = DPD_fit$opt_MDPDE_sigma 
    
    H_Score_DPD[j,] = DPD_fit$H_Score
    
    
  } #end of for (j in 1:n_alpha){ 
  # __________________________________________________________________
  #                       Optimum alpha
  # __________________________________________________________________
  summary = matrix(NA, 6, 2)
  Trimmed_RMPE = paste0(round(trim*100), "% Trimmed RMPE")
  rownames(summary) = c("Opt. DPD alpha", "Opt. Penalty lambda", Trimmed_RMPE, 
                        "MAD", "Dimension Reduction", "Model SD (sigma)")
  colnames(summary) = paste("Robust", IC_type)
  
  opt_alpha_index = apply(H_Score_DPD, MARGIN = 2, which.min) 
  summary["Opt. DPD alpha", ] = alpha[opt_alpha_index]
  summary["Opt. Penalty lambda", ] = diag(lambda_opt_DPD[opt_alpha_index,])
  summary[Trimmed_RMPE, ] = diag(RMPE_trim_DPD[opt_alpha_index,])
  summary["MAD", ] = diag(MAD_DPD[opt_alpha_index,])
  summary["Dimension Reduction", ] = diag(dim_reduction_DPD[opt_alpha_index,])
  summary["Model SD (sigma)", ] = diag(sigma_DPD[opt_alpha_index,])
  
  beta = beta_DPD[1, ,]
  for (i in 1:n.IC)
    beta[, i] = beta_DPD[opt_alpha_index[i], , i]
  
  beta = Matrix(beta, sparse = TRUE)
  
  rm(std_X_global, X_global0, XX_inv_std_global, beta_global, group_std_global, 
     sigma2_global, envir = globalenv())
  
  # check DPD_fail and df_summary
  
  return(list(summary = summary, beta=beta))
}

# ______________________________________________________________________________
#                             Sub-functions
#                     ______________________________
#
# ______________________________________________________________________________
# ______________________________________________________________________________
#                      Optimum lambda (penalty)
# ______________________________________________________________________________
Opt_penalty = function(X, Y, sigma_unbiased, alpha, lambda_vec, gamma, 
                       use_global=FALSE, 
                       IC_type = c("AIC", "BIC", "EBIC", "Cp"),
                       tol=1e-4, maxIter=100, verbose=TRUE){
  # Selects the optimum lambda using "AIC", "BIC", "EBIC", "Cp"
  
  n = length(Y)
  p = ncol(X)
  
  n.lambda = length(lambda_vec)
  n.IC = length(IC_type)
  summary_tab = matrix(NA, n.lambda, 3 + n.IC)
  colnames(summary_tab) = c("lambda", "df", "converged", IC_type)
  summary_tab[,1] = lambda_vec
  
  # entire lambda path
  for (i in 1:n.lambda){
    
    MDPDE_val = MDPDE_SCAD(X=X, Y=Y, alpha=alpha, lambda=lambda_vec[i], 
                           gamma=gamma, use_global=use_global,
                           tol = tol, maxIter=maxIter)
    
    summary_tab[i,"df"] = MDPDE_val$df
    # summary_tab[i,"p"] = MDPDE_val$p.df
    summary_tab[i,"converged"] = MDPDE_val$converged
    
    IC_val = IC_fun(n=n, p=p, sigma_unbiased=sigma_unbiased, 
                    beta_hat=MDPDE_val$beta, sigma_hat=MDPDE_val$sigma,
                    X=X, Y=Y, alpha=alpha, IC_type=IC_type)
    
    summary_tab[i, IC_type] = IC_val
  }
  
  # printing degrees of freedom
  df = summary_tab[,"df"]
  diff_df = diff(df)
  max_diff_index = which.max(diff_df)
  if (verbose) cat(sprintf("\tdf [%d, %d, %d]: max jump %d [%d, %d]\n", 
                           min(df), max(df[-n.lambda]), df[n.lambda], 
                           max(diff_df), df[max_diff_index],
                           df[max_diff_index+1]))
  df_summary = c(min(df), max(df[-n.lambda]), max(diff_df), df[max_diff_index],
                 df[max_diff_index+1])
  names(df_summary) = c("min", "max", "max_diff", "Lower", "Upper")
  
  # Optimum estimates
  opt_lambda = opt_MDPDE_sigma = H_Score  = rep(NA, n.IC)
  opt_MDPDE_beta = matrix(NA, length(MDPDE_val$beta), n.IC)
  MDPDE_y_hat = matrix(NA, n, n.IC)
  X_col_selected = matrix(NA, p, n.IC)
  
  names(opt_lambda) = names(opt_MDPDE_sigma) = names(H_Score) = 
    colnames(opt_MDPDE_beta) = colnames(X_col_selected) = IC_type
  
  # Calculating the optimum ones
  for (i in 1:n.IC){
    opt_lambda[i] = lambda_vec[which.min(summary_tab[,IC_type[i]])]
    
    MDPDE_val = MDPDE_SCAD(X=X, Y=Y, alpha=alpha, lambda=opt_lambda[i], 
                           gamma=gamma, use_global=use_global,
                           tol = tol, maxIter=maxIter)
    
    opt_MDPDE_beta[,i] = MDPDE_val$beta
    opt_MDPDE_sigma[i] = MDPDE_val$sigma
    X_col_selected[,i] = MDPDE_val$X_col_selected
    
    H_Score[i] = H_Score_fun(beta=MDPDE_val$beta, sigma=MDPDE_val$sigma,
                             X=MDPDE_val$X_global, Y=Y, alpha=alpha)
  }
  
  return(list(summary_tab=summary_tab, opt_lambda=opt_lambda, 
              df_summary=df_summary, opt_MDPDE_beta=opt_MDPDE_beta, 
              opt_MDPDE_sigma=opt_MDPDE_sigma, X_col_selected=X_col_selected, 
              H_Score=H_Score))
}
# ______________________________________________________________________________
#                        Generalized AIC
# ______________________________________________________________________________
RAIC_fun <- function(beta, sigma, X, Y, alpha){
  # Generalized AIC using DPD 
  #
  # beta, sigma: estimators of the regression coefficients and sigma
  # alpha: DPD parameter
  # X, Y: data (no column of 1 for intercept in X)
  
  # density
  f_theta_Y = dnorm(Y - cbind(1, X) %*% beta, 0, sigma)
  
  #first term
  if (alpha==0)
    first_term = - sum(log(f_theta_Y[f_theta_Y>1e-200])) else
      first_term = - ((1+alpha)/alpha) * sum(f_theta_Y^alpha)
    
    #number of non-zero coefficients in beta
    df = sum(beta != 0)
    
    #xi_alpha
    xi_alpha = (2*pi)^(-alpha/2) * sigma^(-alpha-2) * (1+alpha)^(-3/2)
    xi_2alpha = (2*pi)^(-alpha) * sigma^(-2*alpha-2) * (1+2*alpha)^(-3/2)
    #eta_alpha
    eta_alpha = (1/4) * (2*pi)^(-alpha/2) * sigma^(-alpha-4) * 
      (2+alpha^2) * (1+alpha)^(-5/2)
    #eta_2alpha
    eta_2alpha = (1/4) * (2*pi)^(-alpha) * sigma^(-2*alpha-4) * 
      (2+4*alpha^2) * (1+2*alpha)^(-5/2)
    
    second_term1 = 2 * df * xi_2alpha/xi_alpha 
    second_term2 = (eta_2alpha - alpha^2 * xi_alpha^2/4)/eta_alpha
    
    RAIC = first_term + second_term1 + second_term2
    return(RAIC)
    
}
# ______________________________________________________________________________
#                       Information Criteria
# ______________________________________________________________________________
IC_fun = function(n, p, sigma_unbiased, beta_hat, sigma_hat,
                  X, Y, alpha, IC_type=c("AIC", "Cp")){
  # Information Criteria using DPD
  
  # n: sample size
  # p: number of covariates 
  # beta_hat, sigma_hat: the DPD estimator using alpha and lambda
  # sigma_unbiased: an unbiased estimator of sigma (preferably using a 
  #         robust alpha with lambda=0) 
  
  #non-zero indices in beta
  df = sum(beta_hat != 0) 
  log_RSS = log(n*sigma_hat^2)
  
  AIC = RAIC_fun(beta=beta_hat, sigma=sigma_hat, X=X, Y=Y, alpha=alpha)
  Cp = n * (sigma_hat/sigma_unbiased)^2 - n + 2*df
  
  IC_vals = c(AIC,  Cp)
  names(IC_vals) = c("AIC", "Cp")
  
  return(IC_vals[IC_type])
}
# ______________________________________________________________________________
#                 H-score for selecting alpha
# ______________________________________________________________________________
H_Score_fun = function(beta, sigma, X, Y, alpha){
  # H-score function for the DPD
  # beta, sigma: estimators of the regression coefficients and sigma
  # X, Y: data
  # alpha: DPD parameter
  
  #gives a large value for a negative sigma
  if (sigma<=0) stop("sigma must be positive.") #return(1e+100)
  
  # density
  Y_Xbeta = Y - cbind(1,X) %*% beta
  f_theta_Y = dnorm(Y_Xbeta/sigma)
  
  D_logL_y = - ((1+alpha)/sigma^2) * Y_Xbeta * f_theta_Y^alpha
  D2_logL_y = ((1+alpha)/sigma^4) * 
    ( alpha * Y_Xbeta^2 - sigma^2 ) * f_theta_Y^alpha
  
  H_score = sum(2 * D2_logL_y + D_logL_y^2)
  
  if (is.na(H_score) | !is.finite(H_score) ){
    stop("wrong H-score.")
  }
  return(H_score)
}
# ______________________________________________________________________________
#             Robust or Rank sure independence screening (RSIS)
# ______________________________________________________________________________
Robust_SIS = function(X, Y, SIS_prop=0.2, SIS_rho.min=NULL, alpha=0.2, 
                      beta_init=NULL, sigma_init=NULL, is.rank=FALSE, 
                      tol = 1e-4, maxIter=100, verbose=TRUE){  
  # Robust or Rank sure independence screening (RSIS)
  
  # X, Y: data
  # SIS_prop: proportion (p) of variables relative to the sample size (n) to 
  #             be used (i.e, n*p) after SIS screening 
  # SIS_rho.min: threshold for correlation coefficient for SIS screening 
  # alpha: DPD parameter
  # beta_init, sigma_init: initial values for MDPDE
  # is.rank: (logical) if true, screening is based on rank SIS from Li et al (2011)
  # tol: tolerance limit for the iteration
  # maxIter: maximum number of iteration
  
  n = length(Y)
  p = ncol(X)
  
  if (is.null(SIS_prop) & is.null(SIS_rho.min)){
    if (n < p) {
      stop("n < p: Use eigher non-null SIS_rho.min or SIS_prop for SIS.")
    } else return(1:p)
  }
  
  if (!is.null(SIS_prop)){
    if (SIS_prop <= 0 | SIS_prop > 1) stop("SIS_prop must be in (0,1].")
    
    # reduced dimension 
    d = round(n * SIS_prop)
    if (d >= p) {
      if (verbose) cat(sprintf("\tNo SIS: d >= p.\n"))
      return(1:p)
    }
  }
  
  if (!is.null(SIS_rho.min)){
    if (SIS_rho.min < 0 | SIS_rho.min > 1) stop("SIS_rho.min must be in [0,1].")
    
    if (SIS_rho.min == 0) {
      if (n < p) {
        stop("n < p: Use non-zero SIS_rho.min or SIS_prop for SIS.")
      } else return(1:p)
    }
  }
  
  beta_hat = rep(NA, p)
  
  for (i in 1:p){
    if (is.rank){
      beta_hat[i] = cor(X[, i], Y, method = "kendall")
    } else
      beta_hat[i] = MDPDE_lm(X=scale(X[, i]), Y=Y, alpha=alpha, 
                             beta_init=beta_init[i+1], sigma_init=sigma_init, 
                             tol = tol, maxIter=maxIter, verbose=verbose)$beta[2]
  }
  
  RSIS_var = NULL
  if (!is.null(SIS_rho.min)){
    if (alpha == 0) sigma_hat = sd(Y) else sigma_hat = mad(Y)
    RSIS_var = which(abs(beta_hat) * sigma_hat > SIS_rho.min)
    
    #checking if further dim. red. is needed
    if (length(RSIS_var) > n) RSIS_var = NULL 
    if (!is.null(SIS_prop)){
      if (length(RSIS_var) > d) RSIS_var = NULL 
    }
  }  
  
  if (is.null(RSIS_var)){  
    if (is.null(SIS_prop)) stop("n < p: Reduce SIS_rho.min for SIS")
    RSIS_var = sort(abs(beta_hat), decreasing = TRUE, index.return = TRUE)$ix[1:d]
  }
  
  if (length(RSIS_var) == 0){
    stop("No varaiable is selected.")
  }
  
  return(RSIS_var)
}
# ______________________________________________________________________________
#                 MDPDE for linear model (unpenalized)
# ______________________________________________________________________________
MDPDE_lm = function(X, Y, alpha=0.2, beta_init=NULL, sigma_init=NULL,
                    tol = 1e-4, maxIter=100, verbose=TRUE){  
  # MDPDE of beta and sigma in normal regression model
  
  # X, Y: data
  # alpha: DPD parameter
  #  beta_init, sigma_init: initial values
  # tol: tolerance limit for the iteration
  # maxIter: maximum number of iteration
  
  n = length(Y)
  X = cbind(1, X)
  #initial values
  if (is.null(beta_init) | is.null(sigma_init)){
    beta = ginv(t(X) %*% X) %*% t(X) %*% Y
    sigma2 = mean((Y - X %*% beta)^2)
  } else {
    beta = beta_init
    sigma2 = sigma_init^2
  }
  
  if (alpha==0) return(list(beta=beta, sigma=sqrt(sigma2), converged=TRUE))
  
  for (iter in 1:maxIter){
    
    beta_old = beta; sigma2_old = sigma2
    
    # density (beta in orthogonal scale)
    f_theta_Y = as.vector(dnorm(Y - X %*% beta, 0, sqrt(sigma2)))
    f_alpha = f_theta_Y^alpha
    D_alpha = diag(f_alpha)
    
    beta = ginv(t(X) %*% D_alpha %*% X) %*% t(X) %*% D_alpha %*% Y
    
    sigma2_1st = sum((Y - X %*% beta)^2 * f_alpha)
    sigma2_2nd = (n * alpha)/( (2*pi)^(alpha/2) * sigma2^(alpha/2 - 1) *
                                 (1 + alpha)^(3/2) )
    sigma2 = (sigma2_1st + sigma2_2nd)/(sum(f_alpha))
    
    if (sigma2 <= 0) stop("sigma must be positive.") 
    if (max(c(abs(beta - beta_old), abs(sigma2 - sigma2_old))) < tol) break
    
  }
  
  converged = ifelse(iter == maxIter, FALSE, TRUE)
  if (!converged & verbose) 
    cat(sprintf("MDPDE did not converge for alpha=%g.\n", alpha))
  
  return(list(beta=beta, sigma=sqrt(sigma2), converged=converged))
}
# ______________________________________________________________________________
#                             MDPDE for SCAD
# ______________________________________________________________________________
MDPDE_SCAD = function(X, Y, alpha=0.2, lambda=0, gamma=3.7, 
                      use_global=FALSE, tol = 1e-4, maxIter=100, 
                      verbose=TRUE){  
  # Penalized MDPDE for SCAD using iterative algorithm
  
  # X, Y: data
  # alpha: DPD parameter
  # lambda: penalty parameter
  # gamma: parameter for SCAD
  # use_global: (logical) if TURE, uses global X_global0, group_std_global, etc.
  # tol: tolerance limit for the iteration
  # maxIter: maximum number of iteration
  
  n = length(Y)
  p = ncol(X)
  
  if (!use_global) Create_global(X)
  
  #initial values
  if (alpha==0 | !exists("beta_global")){
    beta = XX_inv_std_global %*% Y
    sigma2 = mean((Y - cbind(1, std_X_global$X) %*% beta)^2)
    # beta = LAD(y=Y, X = std_X_global$X, intercept = TRUE)
    # sigma2 = sum((Y - cbind(1,std_X_global$X) %*% beta)^2) / (n - sum(beta != 0))
    # if (sigma2 < 0) sigma2 = sum((Y - cbind(1,std_X_global$X) %*% beta)^2) / n
  } else {
    beta = beta_global
    sigma2 = sigma2_global
  }
  
  if (exists("beta_global")){
    if (ncol(std_X_global$X) != length(beta_global[-1])){
      beta = XX_inv_std_global %*% Y
      sigma2 = mean((Y - cbind(1, std_X_global$X) %*% beta)^2)
    }
  }
  
  for (iter in 1:maxIter){
    
    beta_old = beta; sigma2_old = sigma2
    
    # density (beta in orthogonal scale)
    f_theta_Y = as.vector(dnorm(Y - cbind(1, std_X_global$X) %*% beta, 
                                0, sqrt(sigma2)))
    f_alpha = f_theta_Y^alpha
    D_alpha_half = diag(sqrt(f_alpha))
    
    XtD = D_alpha_half %*% std_X_global$X 
    ort_X = orthogonalize_X(X=XtD, group=group_std_global)
    XX = cbind(1, ort_X$X)
    group = ort_X$group
    beta_DPD = ginv(t(XX) %*% XX) %*% t(XX) %*% D_alpha_half %*% Y
    beta0 = beta_DPD[1]
    
    beta = SCAD(X=XX[,-1], Y=Y, beta_ols=beta_DPD, group=group, 
                lambda=lambda, gamma=gamma, tol=tol, maxIter=maxIter)
    
    beta = ort_X$A %*% beta
    beta = c(beta0, beta)
    
    sigma2_1st = sum((Y - cbind(1, std_X_global$X) %*% beta)^2 * f_alpha)
    sigma2_2nd = (n * alpha)/( (2*pi)^(alpha/2) * sigma2^(alpha/2 - 1) *
                                 (1 + alpha)^(3/2) )
    sigma2 = (sigma2_1st + sigma2_2nd)/(sum(f_alpha))
    
    if (sigma2 <= 0) stop("sigma must be positive.") 
    #if (max(c(abs(beta - beta_old), abs(sigma2 - sigma2_old))) < tol) break
    if (max( c(abs(beta - beta_old)/sum(abs(beta_old) ), 
               abs(sigma2 - sigma2_old)/sigma2 )) < tol) break
    
  }
  
  converged = ifelse(iter == maxIter, FALSE, TRUE)
  if (!converged & verbose) 
    cat(sprintf("MDPDE did not converge for alpha=%g.\n", alpha))
  if (converged){
    # global variables for initial parameters
    beta_global <<- beta
    sigma2_global <<- sigma2
  }
  
  # back to the original scale
  beta = beta[-1]/std_X_global$scale
  beta0 = beta_DPD[1] - sum(std_X_global$center * beta)
  beta1 = rep(0, p)
  beta1[std_X_global$non_const_X] = beta
  beta = c(beta0, beta)
  beta[is.na(beta)] = 0 #occurs when a column of X is 0
  
  #significant covariates
  X_col_selected = (abs(beta[-1]) > 0)
  
  return(list(X_col_selected=X_col_selected, sigma=sqrt(sigma2), df=sum(beta !=0),
              beta=beta, X_global=X_global0, converged=converged))
  
}
# ______________________________________________________________________________
#         Defining some global variables to avoid the same calculations
# ______________________________________________________________________________
Create_global = function(X){ 
  # Global variables to avoid recalculating 
  
  # X: data
  
  n = nrow(X)
  p = ncol(X)
  
  group = 1:p
  X_global = X
  
  X_global0 <<- X_global # renaming may be avoided
  
  std_X_global <<- standardize_X(X_global)
  XX = cbind(1, std_X_global$X)
  XX_inv_std_global <<- ginv(t(XX) %*% XX) %*% t(XX)
  non_const_X_global = std_X_global$non_const_X
  group_std_global <<- group[non_const_X_global]
  
}
# ______________________________________________________________________________
#                          Scale and center X  
# ______________________________________________________________________________
standardize_X = function(X) {
  # center and scale X
  
  # group: (numeric vector) group membership, must be ordered
  
  # returns: 
  #      X: scaled design matrix
  #      center: means of columns of the original X
  #      group: regrouping if columns of X_j are (near) constant
  
  n = nrow(X)
  p = ncol(X)
  mu = colMeans(X) #apply(X, 2, median) #
  X = X - matrix(mu, n, p, byrow = TRUE)
  sd_x = sqrt(colMeans(X^2)) #apply(X, 2, mad) #
  non_const_X = nz = which(sd_x > 1e-6)      # non-constant columns
  if (length(nz) != p) {
    X = X[, nz, drop=FALSE]
    mu = mu[nz]
    sd_x = sd_x[nz]
  }
  X = X %*% diag(1/sd_x)
  
  return(list(X=X, center=mu, scale=sd_x, non_const_X=non_const_X))
}
# ______________________________________________________________________________
#                      Orthogonalize X 
# ______________________________________________________________________________
orthogonalize_X = function(X, group) {
  # orthogonalize (group) X
  
  # group: (numeric vector) group membership, must be ordered
  
  # returns: 
  #      X: orthogonalized (group) matrix
  #      A: unorthogonalized (group) matrix
  #      group: regrouping if columns of X_j are dependent
  
  n = nrow(X)
  XX = matrix(0, n, ncol(X))
  group_number = unique(group) 
  n_group = length(group_number)
  A = vector("list", n_group)
  
  for (j in 1:n_group) {
    ind = (group==group_number[j])
    if (sum(ind)==1) {
      XX[,ind] = X[,ind]
      A[[j]] = 1
      next
    }
    SVD = svd(X[, ind, drop=FALSE], nu=0)
    r = which(SVD$d > 1e-10)
    A[[j]] = sweep(SVD$v[, r, drop=FALSE], 2, sqrt(n)/SVD$d[r], "*")
    # SVD$v[, r, drop=FALSE] %*% diag(sqrt(n)/SVD$d) if |r| > 1
    XX[, ind][, r] <- X[, ind] %*% A[[j]]
    
  }
  
  nz = !apply(XX==0, 2, all) # non-zero columns
  XX = XX[, nz, drop=FALSE]
  group = group[nz]
  A = as.matrix(bdiag(A))
  return(list(X=XX, A=A, group=group))
}
# ______________________________________________________________________________
#                             SCAD 
# ______________________________________________________________________________
SCAD = function(X, Y, beta_ols, group, lambda, gamma=3.7, 
                tol=1e-7, maxIter=100){
  # SCAD estimate for beta
  
  # X, Y: data
  # beta_ols: OLS estimator of beta
  # group: (numeric vector) group membership
  # lambda: penalty parameter
  # gamma: SCAD parameter
  # tol, maxIter: convergence criteria 
  
  group_number = unique(group)
  n = length(Y)
  
  # residual
  r = Y - beta_ols[1] - X %*% beta_ols[-1]
  beta_hat = beta_ols[-1]
  
  for (j in 1:maxIter){
    beta_hat_old = beta_hat
    
    beta_scad = numeric(0)
    for (i in group_number){
      
      index_i = (group == i)
      z = beta_hat[index_i] + t(r)  %*% X[,index_i, drop=FALSE]/n
      z_norm = sqrt(sum(z^2))
      lambda1 = sqrt(length(z)) * lambda
      
      if (z_norm <= lambda1) {
        beta_i = matrix(0, 1, length(z))
      } else if (z_norm <= 2*lambda1){
        beta_i = (z/z_norm) * (z_norm - lambda1) 
      } else if (z_norm <= gamma*lambda1){
        beta_i = (z/z_norm) * (z_norm-gamma*lambda1/(gamma-1))/(1-1/(gamma-1))
      } else 
        beta_i = z
      
      beta_scad = c(beta_scad, beta_i)
      r = r - X[,index_i, drop=FALSE] %*% t(beta_i - beta_hat[index_i])
    }
    
    beta_hat = beta_scad
    if (sum(abs(beta_hat - beta_hat_old)) < tol) break
  }
  
  #small component due to tol limit
  if (any(abs(beta_hat[beta_hat !=0])<1e-4))
    beta_hat[abs(beta_hat) < tol]=0
  
  return(beta_hat)
}
# ______________________________________________________________________________
#                    Setup lambda path for DPD SCAD
# ______________________________________________________________________________
SetupLambda = function(X, Y, gamma=3.7, 
                       alpha=0.2, beta_init=NULL, sigma_init=NULL,
                       n.lambda=100, lambda.min=NULL, lambda.max=NULL,
                       use_global=FALSE, tol=1e-4, maxIter=100,
                       verbose=TRUE){
  # Generates lambda path from SCAD algorithm
  
  # X, Y: data
  # gamma: SCAD parameter
  # use_global: (logical) if TURE, uses global X_global, group, etc.
  # n.lambda: returns n.lambda lambda values
  # if both lambda.min and lambda.max are provided, then it calculates
  # lambda path in log scale 
  
  if (!(is.null(lambda.min) | is.null(lambda.max))){
    lambda_vec = c(exp(seq(log(lambda.max), log(lambda.min), 
                           length.out = n.lambda-1)), 0)
    return(lambda_vec)
  }
  
  if (!use_global) Create_global(X)
  
  #MDPDE
  MDPDE = MDPDE_lm(X=std_X_global$X, Y=Y, alpha=alpha, beta_init=beta_init,
                   sigma_init=sigma_init, tol=tol, maxIter=maxIter, 
                   verbose=verbose)
  beta_hat = MDPDE$beta
  sigma_hat = MDPDE$sigma
  
  # density (beta in orthogonal scale)
  f_theta_Y = as.vector(dnorm(Y - cbind(1, std_X_global$X) %*% beta_hat, 
                              mean=0, sd=sigma_hat))
  f_alpha = f_theta_Y^alpha
  D_alpha_half = diag(sqrt(f_alpha))
  
  XtD = D_alpha_half %*% std_X_global$X 
  ort_X = orthogonalize_X(X=XtD, group=group_std_global)
  XX = ort_X$X
  group = ort_X$group
  
  group_number = unique(group)
  n = length(Y)
  r = Y - beta_hat[1] #mean(Y) 
  
  z_max = 0
  for (i in group_number){
    index_i = (group == i)
    z = t(r)  %*% XX[,index_i, drop=FALSE]/n
    z_norm = sqrt(sum(z^2))
    z_max_temp = z_norm / sqrt(length(z)) 
    if (z_max_temp > z_max) z_max = z_max_temp
  }
  
  lambda.max = z_max
  lambda.min = {if (nrow(XX) > ncol(XX)) 1e-4 else .05}
  lambda_vec = exp(seq(log(lambda.max), log(lambda.min*lambda.max), 
                       length=n.lambda))
  
  df = rep(NA, n.lambda)
  for (i in 1:n.lambda){
    df[i] = MDPDE_SCAD(X=X, Y=Y, alpha=alpha, lambda=lambda_vec[i], 
                       gamma=gamma, use_global=use_global, 
                       tol = tol, maxIter=maxIter, verbose=verbose)$df
  }
  
  min_index = (df == max(df))
  lambda.min = min(lambda_vec)
  if (sum(min_index) > 2) lambda.min =  lambda_vec[min_index][2]
  if (df[n.lambda] != max(df)) 
    lambda.min = {if (nrow(XX) > ncol(XX)) 1e-5 else .01}
  
  lambda.max = max(lambda_vec)
  max_index = (df == 1)
  if (sum(max_index) > 2) lambda.max =  lambda_vec[max_index][sum(max_index)]
  
  lambda_vec = c(exp(seq(log(lambda.max), log(lambda.min), 
                         length.out = n.lambda-1)), 0)
  
  return(lambda_vec)
}
# ______________________________________________________________________________
#                       Prediction for additive model
# ______________________________________________________________________________
predict_y = function(X, beta_hat_mat=NULL, Y_test, X_test, trim=0){
  # predict non-parametric additive model
  # X: covariates for the training data
  # beta_hat: matrix containing beta estimates (one column for one IC)
  # g_hat: estimate of g functions at X
  # X_test: covariates for the test data
  # Y_test: response for the test data
  # trim: trimming proportion
  
  # returns predicted values of y and trimmed root mean prediction errors
  
  n.test = nrow(X_test)
  p = ncol(X)
  n.IC = ncol(beta_hat_mat) 
  
  y_hat_mat = matrix(0, n.test, n.IC)
  RMPE_trim = rep(NA, n.IC)
  
  for (i in 1:n.IC){
    y_hat_mat[,i] = cbind(1, X_test) %*% beta_hat_mat[,i] 
    RMPE_trim[i] = sqrt(upper.trim.mean((y_hat_mat[,i] - Y_test)^2, trim=trim))
  }
  
  return(list(y_hat_mat=y_hat_mat, RMPE_trim=RMPE_trim))
}
# ______________________________________________________________________________
#                Upper trimmed mean
# ______________________________________________________________________________
upper.trim.mean = function(x,trim) {
  #trim: the fraction of observations to be trimmed from the top
  if (trim==0) return(mean(x))
  
  x <- sort(x) 
  mean(x[1:floor(length(x)*(1-trim))])
}
# ______________________________________________________________________________
# ______________________________________________________________________________