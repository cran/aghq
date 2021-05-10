## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(aghq)

## ----aghq1,warning = FALSE,message = FALSE------------------------------------
# Stable version:
# install.packages("aghq")
# Or devel version:
# devtools::install_github("awstringer1/aghq")

## ----aghqsimdata1-------------------------------------------------------------
set.seed(84343124)
y <- rpois(10,5) # True lambda = 5, n = 10

## ----aghq22-------------------------------------------------------------------
# Define the log posterior, log(pi(eta,y)) here
logpietay <- function(eta,y) {
  sum(y) * eta - (length(y) + 1) * exp(eta) - sum(lgamma(y+1)) + eta
}
# The objective function has to take a single argument, eta, so pass
# in the data:
objfunc <- function(x) logpietay(x,y)
# Simple numerical derivatives suffice here:
objfuncgrad <- function(x) numDeriv::grad(objfunc,x)
objfunchess <- function(x) numDeriv::hessian(objfunc,x)
# Now create the list to pass to aghq()
funlist <- list(
  fn = objfunc,
  gr = objfuncgrad,
  he = objfunchess
)

## ----doaghq1------------------------------------------------------------------
# AGHQ with k = 3
# Use eta = 0 as a starting value
thequadrature <- aghq::aghq(ff = funlist,k = 3,startingvalue = 0)

## ----doaghq2------------------------------------------------------------------
summary(thequadrature)
plot(thequadrature)

## ----aghqobjecttruth1---------------------------------------------------------
# The normalized posterior:
thequadrature$normalized_posterior$nodesandweights
# The log normalization constant:
options(digits = 6)
thequadrature$normalized_posterior$lognormconst
# Compare to the truth: 
lgamma(1 + sum(y)) - (1 + sum(y)) * log(length(y) + 1) - sum(lgamma(y+1))
# Quite accurate with only n = 10 and k = 3; this example is very simple.
# The mode found by the optimization:
thequadrature$optresults$mode
# The true mode:
log((sum(y) + 1)/(length(y) + 1))
options(digits = 3)

## ----postforlambda1,fig.cap="Approximate (---) and true (- - -) posterior for $\\lambda$"----
# Compute the pdf for eta
pdfwithlambda <- compute_pdf_and_cdf(
  thequadrature$marginals[[1]],
  transformation = list(
    totheta = function(x) log(x),
    fromtheta = function(x) exp(x)
  )
)
# Plot along with the true posterior
# lambdapostplot <- pdfwithlambda %>%
#   ggplot(aes(x = transparam,y = pdf_transparam)) +
#   theme_classic() +
#   geom_line() +
#   stat_function(fun = dgamma,
#                 args = list(shape = 1+sum(y),rate = 1 + length(y)),
#                 linetype = 'dashed') +
#   labs(x = expression(lambda),y = "Density")

with(pdfwithlambda,plot(transparam,pdf_transparam,type='l',xlab = expression(lambda),ylab = "Density"))
with(pdfwithlambda,lines(transparam,dgamma(transparam,shape = 1+sum(y),rate = 1+length(y)),lty='dashed'))


## ----postmean1----------------------------------------------------------------
options(digits = 6)
# Check if the posterior integrates to 1, by computing the "moment" of "1":
compute_moment(thequadrature$normalized_posterior,
               ff = function(x) 1)
# Posterior mean for eta:
compute_moment(thequadrature$normalized_posterior,
               ff = function(x) x)
# Posterior mean for lambda = exp(eta)
compute_moment(thequadrature$normalized_posterior,
               ff = function(x) exp(x))
# Compare to the truth:
(sum(y) + 1)/(length(y) + 1)
options(digits = 3)

## ----getquants1---------------------------------------------------------------
# Get the quantiles for eta:
etaquants <- compute_quantiles(
  thequadrature$marginals[[1]],
  q = c(.01,.25,.50,.75,.99)
)
# The quantiles for lambda = exp(eta) are obtained directly:
exp(etaquants)
# and compared with the truth:
qgamma(c(.01,.25,.50,.75,.99),shape = 1+sum(y),rate = 1+length(y))
# Can be used to get approximate posterior samples:
M <- 1e04 # 10,000 samples
samps <- compute_quantiles(thequadrature$marginals[[1]],runif(M))
hist(exp(samps),breaks = 50,freq=FALSE)
with(pdfwithlambda,lines(transparam,dgamma(transparam,shape = 1+sum(y),rate = 1+length(y)),lty='dashed'))

