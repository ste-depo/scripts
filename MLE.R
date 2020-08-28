# init = starting parameters
# model = function (takes "init" as first argument) and originates the modeled "mu"
# mu = the (mean) experimental observation
# var = variance of experimental observation
# range = multiplied by each parameter of the model returns the range where CI are computed
# ... = additional parameters given to "model"
MLE <- function(init, model, mu, var, range=c(0,10), ...) {
  LL <- function(x, mu, var) {
    sum(dnorm(x, mean=mu, sd=sqrt(var), log=T))
  }
  LLfun <- function(par, mu, var, ...) LL(x=model(par, ...), mu, var)
  suppressWarnings(opt_res <- optim(init, LLfun, mu=mu, var=var, ..., control=list(fnscale=-1)))
  ML <- opt_res$value
  par <- opt_res$par
  # per ogni parametro stima CI
  ci <- lapply( seq_along(par) , function(i) {
    CIFUN <- function(x, par, i, mu, var, ...) {
      par[i] <- x
      (LL(x=model(par, ...), mu, var) - ML + qchisq(.95, 1)/2)^2
    }
    suppressWarnings(left_ci <- optimize(CIFUN, c(range[1]*par[i], par[i]), par, i, mu, var, ...))
    suppressWarnings(right_ci <- optimize(CIFUN, c(par[i], range[2]*par[i]), par, i, mu, var, ...))
    return(cbind(left=unlist(left_ci), right=unlist(right_ci)))
  })
  # plot? va scelto lungo che parametero plottare nel caso ce ne siano più di 1
  # xx <- seq(range[1]*par[i], range[2]*par[i], length.out = 101)
  # sapply(xx, function(x) LLfun(x, ))
  return(list(opt=unlist(opt_res[c('par','value','convergence')]), ci=ci))
}
