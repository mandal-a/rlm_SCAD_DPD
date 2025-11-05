\name{rlm_SCAD_DPD}
\alias{rlm_SCAD_DPD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robust Linear Model Fitting Using SCAD Penalty and 
DPD Estimator
}
\description{
rlm_SCAD_DPD is used to fit a robust linear (regression) model (rlm) using the SCAD-penalized 
density power divergent (DPD) estimator for high-dimensional data.
}
\usage{
rlm_SCAD_DPD(X=NULL, Y=NULL, data=NULL, Y_index=NULL, X_test=NULL,
             Y_test=NULL, newdata=NULL, validation_prop=NULL, trim=0.05, 
             gamma=3.7, SIS_prop = 0.25, SIS_rho.min=NULL, 
             alpha=seq(0, 0.2, length.out=21), alpha_robust=0.15, 
             verbose = TRUE, tol=1e-4, maxIter=100, n.lambda=100,
             seed=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{x} here~~
input matrix, of dimension n x p; each row is an observation vector.
}

\item{Y}{
response variable.
}


\item{data}{data frame containing the full data.}


\item{Y_index}{(scalar) the column index or (character) column name for Y in data,
         other columns in the data will be in X.
         
         Either (Y, X) or (data, Y_index) must be provided.
       (Y, X) will be in (vector, matrix) format
       (data, Y_index) will be in (data frame, scalar, or character) format.}


\item{X_test}{test data for input matrix.}

\item{Y_test}{test data for response variable.}

\item{newdata}{(data frame) new test data.
         Either (Y_test, X_test) or newdata may be provided (optional)}

\item{validation_prop}{proportion of test data for the validation method.
        If Y_test, X_test, and validation_prop are all NULL, the training mean
        absolute deviation (MAD) and trimmed root mean prediction
        error (RMPE)  are calculated.
        If newdata is provided, validation_prop is ignored.}

\item{trim}{trimming proportion for the trimmed RMPE.}

\item{gamma}{SCAD parameter.}
\item{SIS_prop}{proportion (q) of variables relative to the sample size (n) to
            be used (i.e, n*q) after sure independence screening (SIS).}


\item{SIS_rho.min}{threshold for correlation coefficient for SIS screening.}


\item{alpha}{(vector) set of DPD parameters.}


\item{alpha_robust}{robust alpha for estimating sigma_hat for the RCp (needs tuning).}


\item{n.lambda}{maximum number of lambda values (penalty parameter).}

\item{verbose}{logical value for printing important information.}

\item{tol}{tolerance limit for convergence
              of an MDPDE.}

\item{maxIter}{maximum number of steps for convergence
              of an MDPDE.}

\item{seed}{sets seed if CV is used}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list containing summary and beta
\item{summary }{Summary table for the robust AIC and Cp. The optimum DPD
           parameter is selected adaptively using the H-score. }
\item{beta }{A (sparse) matrix containing estimated regression 
            coefficient (beta) based on the robust AIC and Cp.}
  
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Mandal, A. and Ghosh, S. Robust Variable Selection Criteria for the Penalized 
Regression. Journal of Multivariate Analysis (under review). 
}
\author{
%%  ~~who you are~~
Abhijit Mandal, University of Texas at El Paso, USA (Email: amandal@utep.edu)
}
\note{
Required packages: MASS, expm, Matrix. 

The full search may take a long time to run. 
The number of elements in the `alpha` input may be reduced, and the value of 
`tol` may be increased to speed up computation. 

alpha_robust = 0.15 or 0.2 is a good starting value, but it may be tuned based 
on the optimum alpha from the output summary. SIS_prop may also be tuned.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\section{Warning }{It is an initial version of the function, updated on 
November 5, 2025. Users are requested to take the final published version of 
the code.}

\seealso{
See Also \code{\link{MASS::rlm}}, \code{\link{ncvreg::ncvreg}}. 
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
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
Y = beta[1] + X \%*\% beta[-1] + epsilon

# Penalized DPD Fit
DPD_fit = rlm_SCAD_DPD(Y=Y, X=X)
DPD_fit$summary

non_zero_beta = cbind(beta, DPD_fit$beta)
colnames(non_zero_beta) = c("True beta", "Robust AIC", "Robust Cp")
(non_zero_beta = non_zero_beta[rowSums(abs(non_zero_beta)) != 0,])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
