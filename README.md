<table>
<tbody>
<tr>
<td>rlm_SCAD_DPD</td>
</tr>
</tbody>
</table>

## Robust Linear Model Fitting Using SCAD Penalty and DPD Estimator

**Description**

rlm\_SCAD\_DPD fits a robust linear (regression) model (rlm)
using the SCAD-penalized density power divergent (DPD) estimator for
high-dimensional data.

**Usage**

    rlm_SCAD_DPD(X=NULL, Y=NULL, data=NULL, Y_index=NULL, X_test=NULL,
                 Y_test=NULL, newdata=NULL, validation_prop=NULL, 
                 trim=0.05, gamma=3.7, SIS_prop = 0.25, SIS_rho.min=NULL, 
                 alpha=seq(0, 0.2, length.out=21), alpha_robust=0.15, 
                 verbose = TRUE, tol=1e-4, maxIter=100, n.lambda=100,
                 seed=NULL)

**Arguments**

<table>
<colgroup>
<col style="width: 50%" />
<col style="width: 50%" />
</colgroup>
<tbody>
<tr>
<td>X</td>
<td>input matrix, of dimension n x p; each row is an observation
vector.</td>
</tr>
<tr>
<td>Y</td>
<td>response variable.</td>
</tr>
<tr>
<td>data</td>
<td>data frame containing the full data.</td>
</tr>
<tr>
<td>Y_index</td>
<td>(scalar) the column index or (character) column name for Y in data,
other columns in the data will be in X.<br><br>Either (Y, X) or (data,
Y_index) must be provided. (Y, X) will be in (vector, matrix) format
(data, Y_index) will be in (data frame, scalar, or character)
format.</td>
</tr>
<tr>
<td>X_test</td>
<td>test data for input matrix.</td>
</tr>
<tr>
<td>Y_test</td>
<td>test data for response variable.</td>
</tr>
<tr>
<td>newdata</td>
<td>(data frame) new test data. Either (Y_test, X_test) or newdata may
be provided (optional)</td>
</tr>
<tr>
<td>validation_prop</td>
<td>proportion of test data for the validation method. If Y_test,
X_test, and validation_prop are all NULL, the training mean absolute
deviation (MAD) and trimmed root mean prediction error (RMPE) are
calculated. If newdata is provided, validation_prop is ignored.</td>
</tr>
<tr>
<td>trim</td>
<td>trimming proportion for the trimmed RMPE.</td>
</tr>
<tr>
<td>gamma</td>
<td>SCAD parameter.</td>
</tr>
<tr>
<td>SIS_prop</td>
<td>proportion (q) of variables relative to the sample size (n) to be
used (i.e, n*q) after sure independence screening (SIS).</td>
</tr>
<tr>
<td>SIS_rho.min</td>
<td>threshold for correlation coefficient for SIS screening.</td>
</tr>
<tr>
<td>alpha</td>
<td>(vector) set of DPD parameters.</td>
</tr>
<tr>
<td>alpha_robust</td>
<td>robust alpha for estimating sigma_hat for the RCp (needs
tuning).</td>
</tr>
<tr>
<td>n.lambda</td>
<td>maximum number of lambda values (penalty parameter).</td>
</tr>
<tr>
<td>verbose</td>
<td>logical value for printing important information.</td>
</tr>
<tr>
<td>tol</td>
<td>tolerance limit for convergence of an MDPDE.</td>
</tr>
<tr>
<td>maxIter</td>
<td>maximum number of steps for convergence of an MDPDE.</td>
</tr>
<tr>
<td>seed</td>
<td>sets seed if CV is used</td>
</tr>
</tbody>
</table>

**Value**

A list containing summary and beta

<table>
<colgroup>
<col style="width: 50%" />
<col style="width: 50%" />
</colgroup>
<tbody>
<tr>
<td>summary</td>
<td>Summary table for the robust AIC and Cp. The optimum DPD parameter
is selected adaptively using the H-score.</td>
</tr>
<tr>
<td>beta</td>
<td>A (sparse) matrix containing estimated regression coefficient (beta)
based on the robust AIC and Cp.</td>
</tr>
</tbody>
</table>

**Warning**

It is an initial version of the function, updated on November 5, 2025.
Users are requested to take the final version of the code after it is
published.

**Note**

Required packages: MASS, expm, Matrix.

The full search may take a long time to run. The number of elements in
the ‘alpha’ input may be reduced, and the value of ‘tol’ may be
increased to speed up computation.

alpha\_robust = 0.15 or 0.2 is a good starting value, but it may be
tuned based on the optimum alpha from the output summary. SIS\_prop may
also be tuned.

**Author(s)**

Abhijit Mandal, University of Texas at El Paso, USA (Email:
<amandal@utep.edu>)

**References**

Mandal, A. and Ghosh, S. Robust Variable Selection Criteria for the
Penalized Regression. Journal of Multivariate Analysis (under review).

**See Also**

See Also: MASS::rlm, ncvreg::ncvreg.

**Examples**

Run examples

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
