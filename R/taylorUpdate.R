#' Taylor update of regularization parameters
#'
#' Compute new \eqn{\alpha_1} value in the latent observation model \deqn{Z_{kij} = \alpha_{0k} + \alpha_{1k} s_k(R_i(t_(kij)))+\varepsilon_{kij}}
#'
#' @details
#' The individual trajectory \eqn{(R_i)_{i\leq N}} is observed through \eqn{K > 0} latent processes, up to transformations \eqn{(s_k)_{k\leq K}} observed in timepoints  \eqn{(t_{kij})_{i\leq N,k\leq K,j\leq n_{kij}}} : \deqn{Z_{kij}=\alpha_{k0}+\alpha_{k1} s_k(R_i(t_{kij}))+\varepsilon_{kij}} where \eqn{\varepsilon_k\overset{iid}{\sim} \mathcal N(0,\sigma_k^2)}. The individual trajectories \eqn{(R_i)_{i\leq N}} depend on parameter \eqn{\theta}, in which we add the parameter \eqn{\alpha_0=(\alpha_{k0})_{k\leq K}} to differentiate the regularization parameters \eqn{\alpha_1=(\alpha_{k1})_{k\leq K}} from the others. For more details on the model, see the package description.
#'
#' For previous parameter \eqn{\alpha_1^{(l)} ,\theta^{(l)} }, we have the gradient and hessien of the log-likelihood in \eqn{(\theta^{(l)},\alpha_1^{(l)})} (contains in `ddLL` and `dLL`). Then, we update the \eqn{\alpha_1} parameter, for a fixed penalty parameter \eqn{\lambda} (`lambda`), by optimizing the lasso equation \deqn{\partial_{\alpha_{1}} LL(\theta^{(l)},\alpha_1)+\lambda |\alpha_1|=0}
#' Denoting \deqn{A=\partial_{\alpha_{1k}} LL(\theta^{(l)},\alpha_1)|_{\alpha_1=\alpha_1^{(l)}} - \sum_{l=1,\neq k}^K x_{lk}\alpha^{(l)}_{1l} + \sum_{c=1}^K x_{kc}\alpha_{1c}^{(l)}}
#' The update is given by the formula :
#' \deqn{\alpha_{1k}=\displaystyle\left\{\begin{matrix} \frac{A -\lambda}{x_{kk}} & if~~A<-\lambda, \\ \frac{A+\lambda}{x_{kk}} & if~~A > \lambda, \\  0 & otherwise.  \end{matrix}\right.}
#' with \eqn{-\partial_{\alpha_1}^2 LL(\theta^{(l)},\alpha_1)|_{\alpha_1^{(l)}}=-(x_{lc})_{l,c\leq K}}
#'
#' @param alpha previous estimation of the regulatization parameter \eqn{\alpha_1} ;
#' @param lambda penalty parameter ;
#' @param ddLL Hessienne matrix of the log-likelihood in previous estimated \eqn{\theta^{(l)},\alpha_1^{(l)}} ;
#' @param dLL Gradient of the log-likelihood in previous estimated \eqn{\theta^{(l)},\alpha_1^{(l)}}.
#'
#' @return New regularization parameter \eqn{\alpha_1}.
#' @export
#'
#' @examples
#' # TO DO
taylorUpdate <- function(alpha,lambda,ddLL,dLL){

  # check ; size of alpha and ddLL, dLL

  K = length(alpha)
  X = -ddLL

  value.test = dLL + diag(X)*alpha

  alpha.new = sapply(1:K,FUN = function(k){
    if( value.test[k] < - lambda){
      alpha = (value.test[k] - lambda)/X[k,k]
    }else if(value.test[k] >  lambda ){
      alpha = (value.test[k] + lambda)/X[k,k]
    }else{
      alpha = 0
    }
  })

  return(alpha.new)
}
