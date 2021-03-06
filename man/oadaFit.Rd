% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oadaFit.R
\name{oadaFit}
\alias{oadaFit}
\title{Fit an order of acquisition diffusion analysis (OADA) model}
\usage{
oadaFit(nbdadata, type = "social", startValue = NULL, lower = NULL,
  upper = NULL, interval = c(0, 999), method = "nlminb",
  gradient = T, iterations = 150, standardErrors = "Numeric",
  formula = NULL, coxmeFit = NULL, SLdom = F)
}
\arguments{
\item{nbdadata}{an object of class nbdaData (\code{\link{nbdaData}}) to fit a model to a single diffusion or a list of
nbdaData objects to fit a model to multiple diffusions.}

\item{type}{a string specifying either "social" or "asocial" model. Usually asocial models have all s parameters
constrained =0 and all ILVs affecting only the rate of social learning are removed (i.e. those in the int_ilv slot(s)
of the nbdaData object(s)). However, if a non-zero offset is present on the social transmission component, e.g. when
constaining all s parameters to a specific value using \code{\link{constrainedNBDAdata}}, int_ilv variables are
retained. This situation occurs most commonly when the function is called internally by the \code{\link{profLikCI}}
function.}

\item{startValue}{optional numeric vector giving start values for the maximum likelihood optimization. Length to match
the number of parameters fitted in the model.}

\item{lower}{optional numeric vector giving lower values for the maximum likelihood optimization. Length to match
the number of parameters fitted in the model. By default taken to be 0 for all s parameters and -Inf for coefficients of
ILVs.}

\item{upper}{optional numeric vector giving upper values for the maximum likelihood optimization. Length to match
the number of parameters fitted in the model. By default taken to be Inf for all parameters.}

\item{interval}{currently non-functioning argument: can be ignored.}

\item{method}{character string determining which optimization algorithm is used, defaulting to "nlminb" using the
\code{\link[stats]{nlminb}} function. If set to "both" the optim method \code{\link[stats]{optim}} is also used and the
results returned for both optimization procedures.}

\item{gradient}{logical indicating whether the gradient function should be used during optimization.}

\item{iterations}{numerical determining the maximum iterations to be used during optimization. Increasing this may solve
convergence issues.}

\item{standardErrors}{a string indicating how standard errors should be calculated. Defaults to "Numeric" which uses the
\code{\link[numDeriv]{hessian}} function. In some cases an analytical solution can be provided using "Analytic": this can
increase the speed of the fit for large models. In cases where the analytic solution is not available in the package, the
standard errors default back to numeric. There is no advantage in accuracy to using the analytical solution. If any other
string is provided, no standard errors are returned.}

\item{formula}{a formula can be provided to customise the model if being fitted using the coxme function. At the current
time users are advised to leave this =NULL and the appropriate formula, including all effect specified in the nbdaData
objects, will be built internally.}

\item{coxmeFit}{logical indicating whether the \code{\link[coxme]{coxme}} function should be used to fit the model. This
is set to NULL by default, meaning the \code{\link[coxme]{coxme}} function will be used if random effects are included.
If random effects are included and coxmeFit=F the random effects will be ignored.}

\item{SLdom}{logical determining whether a "social learning dominant" model should be fitted. This is useful in cases
where the user suspects that s parameters are all =Inf- which can occur if individuals with zero connections to informed
individuals only learn when all other naive individuals also have zero connections to informed individuals. This means
that the strength of social learning relative to asocial learning will be estimated as Inf. In such cases
an SLdom model enables the user to judge the size of s parameters for different networks relative to one another.}
}
\value{
If coxmeFit=F or random effects are absent the function returns an object of class oadaFit. If coxmeFit=T or
random effects are included the function returns an object of class oadaFit_coxme.
}
\description{
\code{oadaFit} takes diffusion data in the form of an nbdaData object (\code{\link{nbdaData}}) or a list of nbdaData
objects (for multiple diffusions) and fits an order of acquisition diffusion analysis (OADA) model.
}
\details{
The model is fitted using maximum likelihood methods, for TADA models use \code{\link{tadaFit}}.The ILVs and random
effects included in the model are determined by those present in the nbdaData object(s). All nbdaData objects
must contain the same social networks (assMatrix must match in the third dimension), the same individual level
variables (ILVs) in each of the asoc_ilv, int_ilv and multi_ilv slots and the same random effects in the random_effects
slot. If random effects are included, the model is fitted by calls to the \code{\link[coxme]{coxme}} function in the
coxme package. Random effects are assumed to operate multiplicatively, i.e. affect asocial and social learning
differences among individuals by the same amount. If more complex effects are required then a Bayesian approach is
recommended (not implemented in the NBDA package). trueTies specified in nbdaData object(s) are accounted for unless
random effects are included. This is done by adding the likelihood across all orders of acquisition consistent with the
trueTies.l This is highly computationally intensive and inadvisable for any but a few true ties.
}
\section{oadaFit components}{

The following components of the oadaFit object are of key importance for
interpreting the output: \describe{
  \item{@outputPar}{The maximum likelihood estimates (MLEs) for the model parameters}
  \item{@varNames}{The name of the variable corresponding to each of the parameter estimates. These are numbered so
  the user can easily identify parameters when obtaining confidence intervals using \code{\link{profLikCI}}. The s
  parameters are labelled "Social transmission N" with N giving the number of the network. ILV effects on asocial
  learning are preceded with "Asocial:". ILV effects on social learning are preceded with "Social:". "Multiplicative"
  ILV effects constrained to be equal on asocial and social learning are preceded with "Social=Asocial".}
  \item{@se}{The standard error for each parameter. These can not always be derived so may be NaN. The user is advised
  to get confidence intervals for parameters using \code{\link{profLikCI}}.}
  \item{@aic}{The AIC for the model.}
  \item{@aicc}{The AICc for the model: AIC adjusted for sample size, with sample size taken to be the number of
  acquisition events.}
  \item{@loglik}{The -log-likelihood for the model. Can be used to conduct likelihood ratio tests to test hypotheses.}
}The oadaData object also contains the following components:
\describe{
  \item{@nbdadata}{The data the model is fitted to, as a list of nbdaData objects.}
  \item{@optimisation}{The output of the \code{\link[stats]{nlminb}} optimization alogorithm, useful for assessing
  convergence of the model.}
  \item{@optim}{The output of the \code{\link[stats]{optim}} optimization alogorithm, where used, useful for assessing
  convergence of the model.}
  \item{@hessian}{The hessian matrix- giving the value of the second partial derivatives of the -log-likelihood with
  respect to the model parameters at the maximum likelihood estimators. Used to dervive the standard errors.}
  \item{@type}{The model type: "asocial" or "social".}
  \item{@SLdom}{Logical showing whether a "social learning dominant" model was fitted (see above).}
  }
}

\section{Additional oadaFit_coxme components}{
 In addtion to the components described for oadaFit objects, the following
slots are present in an oadaFit_coxme object:
\describe{
  \item{@REvar}{The estimated variance of the random effects fitted by the \code{\link[coxme]{coxme}} function.}
  \item{@fixedFormula}{The formula for the fixed effects input to the \code{\link[coxme]{coxme}} function. This will
  include any multiplicative ILVs fitted.}
  \item{@randomFormula}{The formula for the random effects input to the \code{\link[coxme]{coxme}} function.}
  }
}

\seealso{
For TADA models use \code{\link{tadaFit}}. To obtain confidence intervals see \code{\link{profLikCI}}. For
further details about OADA see \url{https://www.sciencedirect.com/science/article/pii/S0022519310000081} and
\url{https://royalsocietypublishing.org/doi/full/10.1098/rstb.2016.0418}
}
