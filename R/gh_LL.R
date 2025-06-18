#' Adaptive Gauss-Hermite approximation of log-likelihood derivatives
#'
#' @description
#' Computes Adaptive Gauss-Hermite approximation of the log-likelihood and its derivatives in NLMEM with latent observation processes, see \code{\link{REMixed-package}} for details on the model.
#'
#' @details
#' Based on notation introduced \code{\link{REMixed-package}}. The log-likelihood of the model \eqn{LL(\theta,\alpha_1)} for a set of population parameters \eqn{\theta} and regulatization parameters \eqn{\alpha_1} is estimated using Adaptative Gausse-Hermite quadrature, using  conditional distribution estimation to locate the mass of the integrand. If the project has been initialized as a Monolix project, the user can use \code{\link{readMLX}} function to retrieve all the project information needed here.
#'
#' @param dynFUN function computing the dynamics of interest for a set of parameters. This function need to contain every sub-function that it may needs (as it is called in a foreach loop). The output of this function need to return a data.frame with \code{time} : as first columns and named dynamics in other columns. It must take in input : \itemize{\item\code{y} : a named vector with the initial condition. The names are the dynamics names.
#' \item\code{parms} : a named vector of parameter.
#' \item\code{time} : vector a timepoint.}
#'
#' See \code{\link{dynFUN_demo}}, \code{\link{model.clairon}}, \code{\link{model.pasin}} or \code{\link{model.pk}} for examples.
#' @param y initial condition of the mechanism model, conform to what is asked in \code{dynFUN}.
#' @param mu list of individuals random effects estimation (vector of r.e. need to be named by the parameter names), use to locate the density mass; (optional, see description).
#' @param Omega list of individuals estimated standard deviation diagonal matrix (matrix need to have rows and columns named by the parameter names), use to locate the density mass; (optional, see description).
#' @param theta list of model parameters containing (see details) \itemize{\item\code{phi_pop} : named vector with the population parameters with no r.e. \eqn{(\phi_{l\ pop})_{l\leq L}} (NULL if none) ;
#' \item\code{psi_pop} : named vector with the population parameters with r.e. \eqn{(\psi_{l\ pop})_{l\leq m}} ;
#' \item\code{gamma} : named list (for each parameters) of named vector (for each covariates) of covariate effects from parameters with no r.e. ;
#' \item\code{beta} : named list (for each parameters) of named vector (for each covariates) of covariate effects from parameters with r.e..
#' \item\code{alpha0} : named vector of \eqn{(\alpha_{0k})_{k\leq K}} parameters (names are identifier of the observation model, such as in a Monolix project);
#' \item\code{omega} : named vector of estimated r.e. standard deviation;}
#' (optional, see description).
#' @param alpha1 named vector of regulatization parameters \eqn{(\alpha_{1k})_{k\leq K}}, with identifier of observation model as names, (optional, see description).
#' @param covariates matrix of individual covariates (size N x n). Individuals must be sorted in the same order than in \code{mu} and \code{Omega}, (optional, see description).
#' @param ParModel.transfo named list of transformation functions \eqn{(h_l)_{l\leq m}} and \eqn{(s_k)_{k\leq K}} for the individual parameter model (names must be consistent with \code{phi_pop} and \code{psi_pop}, missing entries are set by default to the identity function ; optional, see description).
#' @param ParModel.transfo.inv Named list of inverse transformation functions for the individual parameter model (names must be consistent with \code{phi_pop} and \code{psi_pop} ; optional, see description).
#' @param Sobs list of individuals trajectories for the direct observation models \eqn{(Y_{pi})_{p \leq P,i\leq N}}. Each element \eqn{i\leq N} of the list, is a list of \eqn{p\leq P} data.frame with time  \eqn{(t_{pij})_{j\leq n_{ip}}} and observations  \eqn{(Y_{pij})_{j\leq n_{ip}}}. Each data.frame is named with the observation model identifiers.
#' @param Robs list of individuals trajectories for the latent observation models \eqn{(Z_{ki})_{k \leq K,i\leq N}}. Each element \eqn{i\leq N} of the list, is a list of \eqn{k\leq K} data.frame with time  \eqn{(t_{kij})_{j\leq n_{ik}}} and observations  \eqn{(Z_{kij})_{j\leq n_{ik}}}. Each data.frame is named with the observation model identifiers.
#' @param Serr named vector of the estimated error mocel constants \eqn{(\varsigma_p)_{p\leq P}} with observation model identifiers as names.
#' @param Rerr named vector of the estimated error mocel constants \eqn{(\sigma_k)_{k\leq K}} with observation model identifiers as names.
#' @param ObsModel.transfo list containing two lists of transformations and two vectors linking each transformations to their observation model name in the Monolix project. The list should include identity transformations and be named \code{S} and \code{R}. The two vectors should be named \code{linkS} and \code{linkR}.
#'
#' Both \code{S} (for the direct observation models) and \code{linkS}, as well as \code{R} (for latent process models) and \code{linkR}, must have the same length.
#'
#' \itemize{
#'   \item\code{S}: a list of transformations for the direct observation models. Each transformation corresponds to a variable \eqn{Y_p=h_p(S_p)}, where the name indicates which dynamic is observed (from \code{dynFUN});  \item\code{linkS} : a vector specifying the observation model names (that is used in the monolix project, \code{alpha1}, etc.) for each transformation, in the same order as in \code{S};
#'
#'   \item\code{R}: similarly, a list of transformations for the latent process models. Although currently there is only one latent dynamic, each \eqn{s_k, k\leq K} transformation corresponds to the same dynamic but may vary for each \eqn{Y_k} observed. The names should match the output from \code{dynFUN}; \item \code{linkR} : a vector specifying the observation model names for each transformation, in the same order as in \code{R}.
#' }
#' @param data output from \code{\link{readMLX}} containing parameters "\code{mu}", "\code{Omega}", "\code{theta}", "\code{alpha1}", "\code{covariates}", "\code{ParModel.transfo}", "\code{ParModel.transfo.inv}", "\code{Sobs}", "\code{Robs}", "\code{Serr}", "\code{Rerr}", "\code{ObsModel.transfo}" extract from a monolix project.
#' @param n number of points per dimension to use for the Gauss-Hermite quadrature rule.
#' @param prune integer between 0 and 1, percentage of pruning for the Gauss-Hermite quadrature rule (default NULL).
#' @param parallel logical, if computation should be done in parallel.
#' @param ncores number of cores to use for parallelization, default will detect the number of cores available.
#' @param onlyLL logical, if only the log-likelihood should be computed (and not \eqn{\partial_{\alpha_1} LL} or \eqn{\partial_{\alpha_1}^2 LL}).
#' @param verbose logical, if progress bar should be printed through the computation.
#'
#' @return A list with the approximation by Gauss-Hermite quadrature of the likelihood \code{L}, the log-likelihood \code{LL}, the gradient of the log-likelihood \code{dLL}, and the Hessian of the log-likelihood \code{ddLL} at the point \eqn{\theta, \alpha} provided.
#' @export
#'
#' @examples
#'\dontrun{
#' project <- getMLXdir()
#'
#'
#' ObsModel.transfo = list(S=list(AB=log10),
#'                         linkS="yAB",
#'                         R=rep(list(S=function(x){x}),5),
#'                         linkR = paste0("yG",1:5))
#'
#' alpha=list(alpha0=NULL,
#'            alpha1=setNames(paste0("alpha_1",1:5),paste0("yG",1:5)))
#'
#' data <- readMLX(project,ObsModel.transfo,alpha)
#'
#' LL <- gh.LL(dynFUN = dynFUN_demo,
#'             y = c(S=5,AB=1000),
#'             ObsModel.transfo=ObsModel.transfo,
#'             data = data)
#'
#' print(LL)
#' }
gh.LL <- function(
    dynFUN,
    y,
    mu=NULL,
    Omega=NULL,
    theta=NULL,
    alpha1=NULL,
    covariates=NULL,
    ParModel.transfo=NULL,
    ParModel.transfo.inv=NULL,
    Sobs=NULL,
    Robs=NULL,
    Serr=NULL,
    Rerr=NULL,
    ObsModel.transfo=NULL,
    data=NULL,
    n=NULL,
    prune=NULL,
    parallel=TRUE,
    ncores=NULL,
    onlyLL=FALSE,
    verbose=TRUE){

  if(is.null(data)){
    test <- sapply(c("mu","Omega","theta","alpha1","covariates","ParModel.transfo","ParModel.transfo.inv","Sobs","Robs","Serr","Rerr","ObsModel.transfo"),FUN=is.null)
    if(any(test))
      stop("Please provide all necessary arguments.")
  }else{
    if((is.null(mu) || is.null(Omega)) & !all(c("mu","Omega") %in% names(data))){
      stop("Please provide mu and Omega if these are missing from data.")
    }
    for(d in 1:length(data)){
      test <- eval(parse(text=(paste0("is.null(",names(data)[d],")"))))
      if(test){
        eval(parse(text=paste0(names(data[d]),"<- data[[d]]")))
      }
    }
  }
  if(is.null(n)){
    if(length(theta$psi_pop)==1){
      n <- 100
    }else if(length(theta$psi_pop)==2){
      n <- 10
    }else{
      n <- 7
    }
  }

  if(parallel){
    if(!is.null(ncores)){
      cluster <- snow::makeCluster(ncores)
    }else{
      cluster <- snow::makeCluster(parallel::detectCores())
    }
    doSNOW::registerDoSNOW(cluster)
  }

  N=length(mu)

  i = 1

  ntasks <- N
  if(verbose){
    pb <- txtProgressBar(max = ntasks, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  }else{
    opts <- NULL
  }

  res = foreach::foreach(i = 1:N,.packages = "REMixed",.options.snow=opts)%dopar%{
    if(0 %in% diag(Omega[[i]])){
      diag(Omega[[i]])[diag(Omega[[i]])==0] <- 10**(-5)
    }
    LLi = gh.LL.ind(mu_i = mu[[i]],
                    Omega_i = Omega[[i]],
                    theta = theta,
                    alpha1 = alpha1,
                    dynFUN = dynFUN,
                    y = y,
                    covariates_i = covariates[i,,drop=F],
                    ParModel.transfo = ParModel.transfo,
                    ParModel.transfo.inv = ParModel.transfo.inv,
                    Sobs_i = lapply(Sobs,FUN=function(S){S[S$id==i,]}),
                    Robs_i = lapply(Robs,FUN=function(R){R[R$id==i,]}),
                    Serr = Serr,
                    Rerr = Rerr,
                    ObsModel.transfo = ObsModel.transfo,
                    ind = i,
                    n = n,
                    prune = prune,
                    onlyLL = onlyLL)


    return(LLi)
  }

  L = prod(sapply(res,FUN=function(ri){ri$Li}))

  LL = sum(sapply(res,FUN=function(ri){ri$LLi}))


  if(verbose)
    close(pb)
  if(parallel){
    snow::stopCluster(cluster)
  }

  if(!onlyLL){
    if(length(alpha1)!=1){
      dLL = rowSums(sapply(res,FUN=function(ri){ri$dLLi}))
      ddLL = matrix(rowSums(sapply(res,FUN=function(ri){ri$ddLLi})),ncol=length(dLL))
    }else{
      dLL = sum(sapply(res,FUN=function(ri){ri$dLLi}))
      ddLL = matrix(sum(sapply(res,FUN=function(ri){ri$ddLLi})),ncol=length(dLL))
    }
    return(list(L=L,LL=LL,dLL=dLL,ddLL=ddLL))
  }else{
    return(LL)
  }
}

gh.LL.ind <- function(
    dynFUN,
    y,
    mu_i=NULL,
    Omega_i=NULL,
    theta=NULL,
    alpha1=NULL,
    covariates_i=NULL,
    ParModel.transfo=NULL,
    ParModel.transfo.inv=NULL,
    Sobs_i=NULL,
    Robs_i=NULL,
    Serr=NULL,
    Rerr=NULL,
    ObsModel.transfo=NULL,
    data=NULL,
    ind = NULL,
    n = NULL,
    prune=NULL,
    onlyLL=FALSE){

  mu <- Omega <- Sobs <- Robs <- covariates <- NULL

  if(is.null(data)){
    test <- sapply(c("mu_i","Omega_i","theta","alpha1","covariates_i","ParModel.transfo","ParModel.transfo.inv","Sobs_i","Robs_i","Serr","Rerr","ObsModel.transfo"),FUN=is.null)
    if(any(test))
      stop("Please provide all necessary arguments.")
  }else{
    if((is.null(mu_i) || is.null(Omega_i)) & !all(c("mu","Omega") %in% names(data))){
      stop("Please provide mu and Omega if these are missing from data.")
    }
    argsNONind = c("theta","alpha1","ParModel.transfo","ParModel.transfo.inv","Serr","Rerr","ObsModel.transfo")
    for(d in 1:length(argsNONind)){
      test <- eval(parse(text=paste0("is.null(",argsNONind[d],")")))
      if(test){
        eval(parse(text=paste0(argsNONind[d],"<- data[[argsNONind[d]]]")))
      }
    }
    argsInd = c("mu","Omega","Robs","Sobs","covariates")
    for(d in 1:length(argsInd)){
      test <- eval(parse(text=paste0("is.null(",paste0(argsInd[d],"_i"),")")))
      if(test){
        eval(parse(text=paste0(argsInd[d],"<- data[[argsInd[d]]]")))

      }
    }
    if((is.null(mu_i) || is.null(Omega_i) || is.null(Robs_i) || is.null(Sobs_i) || is.null(covariates_i)) & is.null(ind)){
      stop("Please provide individual information (mu_i,Omega_i,Sobs_i,Robs_i,covariates_i), or the individual id ind.")
    }
    if(is.null(mu_i)){
      mu_i = mu[[ind]]
    }
    if(is.null(Omega_i)){
      Omega_i = Omega[[ind]]
    }
    if(is.null(Sobs_i)){
      Sobs_i = lapply(Sobs,FUN=function(S){S[S$id==ind,]})
    }
    if(is.null(Robs_i)){
      Robs_i = lapply(Robs,FUN=function(R){R[R$id==ind,]})
    }
    if(is.null(covariates_i)){
      covariates_i = covariates[ind,,drop=F]
    }
  }

  if(is.null(n)){
    if(length(theta$psi_pop)==1){
      n <- 100
    }else if(length(theta$psi_pop)==2){
      n <- 10
    }else{
      n <- 7
    }
  }

  dm = length(mu_i)

  if(dm!=1){
    mh.parm <- amgauss.hermite(n,mu=mu_i,Omega=Omega_i,prune=prune)
    detsq.omega = prod(theta$omega)
    root.omega = diag(1/theta$omega**2)
  }else{
    mh.parm <- agauss.hermite(n,mu = mu_i,sd=sqrt(as.numeric(Omega_i)))
    detsq.omega = as.numeric(theta$omega)
    root.omega = 1/as.numeric(theta$omega)**2
  }

  nd = length(mh.parm$Weights)
  R.sz = length(Rerr)
  S.sz = length(Serr)
  all.tobs = sort(union(unlist(lapply(Robs_i,FUN=function(x){x$time})),unlist(lapply(Sobs_i,FUN=function(x){x$time}))))



  # Need to compute, for each eta_i x individual i, the dynamics of the model
  # eta_i = split(mh.parm$Points,1:nd)[[1]]
  dyn <- setNames(lapply(split(mh.parm$Points,1:nd),FUN=function(eta_i){
    PSI_i  = indParm(theta[c("phi_pop","psi_pop","gamma","beta")],covariates_i,setNames(eta_i,colnames(Omega_i)),ParModel.transfo,ParModel.transfo.inv)
    dyn_eta_i <- dynFUN(all.tobs,y,unlist(unname(PSI_i)))

    return(dyn_eta_i)
  }),paste0("eta_",1:nd))

  # Compute marginal latent density for each individual and value of RE
  R.margDensity <- setNames(lapply(1:nd,FUN=function(ei){
    dyn.ei = dyn[[paste0("eta_",ei)]]

    # .Machine$double.xmin

    decomp_marg = lapply(1:R.sz,FUN=function(k){
      yGk = ObsModel.transfo$linkR[k]
      sig = Rerr[[yGk]]
      tki = Robs_i[[yGk]]$time
      Zki = Robs_i[[yGk]][,yGk]

      a0 = theta$alpha0[[yGk]]
      a1 = alpha1[[yGk]]
      trs = ObsModel.transfo$R[which(ObsModel.transfo$linkR==yGk)]
      Rki = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tki,names(trs)])

      aux <- 1/(sig*sqrt(2*pi))*exp(-1/2*((Zki-a0-a1*Rki)/sig)**2)
      if(any(aux)==0){
        pw <- Rmpfr::mpfr(-1/2*((Zki-a0-a1*Rki)/sig)**2,precBits = 32)

        aux = 1/(sig*sqrt(2*pi))*exp(pw)
      }

      decomp_aux = setNames(sapply(decomp(prod(aux)),as.numeric),c("exponent","mantissa"))

      return(decomp_aux)
    })

    if(any(sapply(decomp_marg,function(x){x["mantissa"]})==0)){
      stop("[Error] in log-likelihood computation, latent part of the marginal likelihood is equal to zero.")
    }

    # decomp_marg = lapply(little_marg,decomp)

    decomp_intermediaire = decomp(prod(sapply(decomp_marg,FUN=function(decomp){decomp["mantissa"]})))

    mantissa_res = decomp_intermediaire["mantissa"]
    exponent_res = sum(sapply(decomp_marg,FUN=function(decomp){decomp["exponent"]})) + decomp_intermediaire["exponent"]

    return(c(exponent_res,mantissa_res))
  }),paste0("eta_",1:nd))

 # Compute marginal density for each individual and value of RE
  if(S.sz!=0){
    S.margDensity <- setNames(sapply(1:nd,FUN=function(ei){
      eta_i = split(mh.parm$Points,1:nd)[[ei]]
      dyn.ei = dyn[[paste0("eta_",ei)]]

      res = prod(sapply(1:S.sz,FUN=function(p){
        Yp = ObsModel.transfo$linkS[p]
        sig = Serr[[Yp]]
        tpi = Sobs_i[[Yp]]$time
        Ypi = Sobs_i[[Yp]][,Yp]

        trs = ObsModel.transfo$S[which(ObsModel.transfo$linkS==Yp)]
        Spi = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tpi,names(trs)])

        prod(1/(sig*sqrt(2*pi))*exp(-1/2*((Ypi-Spi)/sig)**2))
      }))*(1/((2*pi)**(dm/2)*detsq.omega)*exp(-1/2*eta_i%*%root.omega%*%eta_i))
    }),paste0("eta_",1:nd))
  }else{
    S.margDensity = setNames(sapply(split(mh.parm$Points,1:nd),FUN=function(eta_i){
      (1/((2*pi)**(dm/2)*detsq.omega)*exp(-1/2*eta_i%*%root.omega%*%eta_i))
    }),paste0("eta_",1:nd))
  }



  # Compute individual log-Likelihood
  Li.aux = lapply(1:nd,function(ei){
    decomp_intermediaire = decomp(mh.parm$Weights[[ei]]*R.margDensity[[ei]][["mantissa"]]*S.margDensity[[ei]])

    mantissa_res = decomp_intermediaire["mantissa"]
    exponent_res = R.margDensity[[ei]][["exponent"]] + decomp_intermediaire["exponent"]

    return(c(exponent_res,mantissa_res))
  })
  max_exponent <- max(sapply(Li.aux, function(x) x[["exponent"]]))
  decomp.aux <- decomp(sum(sapply(Li.aux,function(li){
    return(li[["mantissa"]]*10**(li[["exponent"]] - max_exponent))
  })))

  Li = c(exponent=(max_exponent+decomp.aux[["exponent"]]),mantissa=decomp.aux[["mantissa"]])
    # sum(mh.parm$Weights*R.margDensity*S.margDensity)
  invLi.aux = decomp(1/Li[["mantissa"]])
  invLi = c(exponent=-Li[["exponent"]]+invLi.aux[["exponent"]],mantissa=invLi.aux[["mantissa"]])
  # ln(p*10**o)=ln(p)+ln(10**o)=ln(p)+o*ln(10)
  LLi = log(Li[["mantissa"]])+Li[["exponent"]]*log(10)

  if(!onlyLL){
    # Compute  gradient of individual log-Likelihood
    # dim.alpha = length(alpha1) # = K = length(Robs_i )
    dfact = lapply(1:R.sz,FUN=function(k){
      f = setNames(sapply(1:nd,FUN=function(ei){
        dyn.ei = dyn[[paste0("eta_",ei)]]
        yGk = ObsModel.transfo$linkR[k]
        sig = Rerr[[yGk]]
        tki = Robs_i[[yGk]]$time
        Zki = Robs_i[[yGk]][,yGk]

        a0 = theta$alpha0[[yGk]]
        a1 = alpha1[[yGk]]
        trs = ObsModel.transfo$R[which(ObsModel.transfo$linkR==yGk)]
        Rki = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tki,names(trs)])

        return(1/(sig**2)*(sum(Rki*(Zki-a0-a1*Rki))))
      }),paste0("eta_",1:nd))
    })

    dLLi = sapply(dfact,FUN=function(f){
      res.aux <- lapply(1:nd,function(ei){
        decomp_intermediaire = decomp(mh.parm$Weights[[ei]]*f[[ei]]*R.margDensity[[ei]][["mantissa"]]*S.margDensity[[ei]])

        mantissa_res = decomp_intermediaire["mantissa"]
        exponent_res = R.margDensity[[ei]][["exponent"]] + decomp_intermediaire["exponent"]

        return(c(exponent_res,mantissa_res))
      })
      max_exponent <- max(sapply(res.aux, function(x) x[["exponent"]]))
      decomp.aux <- decomp(sum(sapply(res.aux,function(ri){
        return(ri[["mantissa"]]*10**(ri[["exponent"]] - max_exponent))
      })))

      aux =c(exponent=(max_exponent+decomp.aux[["exponent"]]),mantissa=decomp.aux[["mantissa"]])

     return(invLi[["mantissa"]]*aux[["mantissa"]]*10**(invLi[["exponent"]]+aux[["exponent"]]))
       # 1/Li*sum(mh.parm$Weights*f*R.margDensity*S.margDensity))
    })

    # Compute hessienne of individual log-Likelihood
    ddLLi = matrix(0,ncol=R.sz,nrow=R.sz)
    # diagonal term
    ddfact = lapply(1:R.sz,FUN=function(k){
      setNames(sapply(1:nd,FUN=function(ei){
        dyn.ei = dyn[[paste0("eta_",ei)]]
        yGk = ObsModel.transfo$linkR[k]
        sig = Rerr[[yGk]]
        tki = Robs_i[[yGk]]$time
        Zki = Robs_i[[yGk]][,yGk]

        a0 = theta$alpha0[[yGk]]
        a1 = alpha1[[yGk]]
        trs = ObsModel.transfo$R[which(ObsModel.transfo$linkR==yGk)]
        Rki = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tki,names(trs)])

        return(
          (1/(sig**2)*(sum(Rki*(Zki-a0-a1*Rki))))**2 - 1/(sig**2)*(sum(Rki**2))
        )
      }),paste0("eta_",1:nd))
    })

    diag(ddLLi) = sapply(ddfact,FUN=function(f){
      res.aux <- lapply(1:nd,function(ei){
        decomp_intermediaire = decomp(mh.parm$Weights[[ei]]*f[[ei]]*R.margDensity[[ei]][["mantissa"]]*S.margDensity[[ei]])

        mantissa_res = decomp_intermediaire["mantissa"]
        exponent_res = R.margDensity[[ei]][["exponent"]] + decomp_intermediaire["exponent"]

        return(c(exponent_res,mantissa_res))
      })
      max_exponent <- max(sapply(res.aux, function(x) x[["exponent"]]))
      decomp.aux <- decomp(sum(sapply(res.aux,function(ri){
        return(ri[["mantissa"]]*10**(ri[["exponent"]] - max_exponent))
      })))

      aux =c(exponent=(max_exponent+decomp.aux[["exponent"]]),mantissa=decomp.aux[["mantissa"]])

      return(invLi[["mantissa"]]*aux[["mantissa"]]*10**(invLi[["exponent"]]+aux[["exponent"]]))
      # return(1/Li*sum(mh.parm$Weights*f*R.margDensity*S.margDensity))
    })-(dLLi)**2


    # non diagonal term
    for(k1 in 1:R.sz){
      for(k2 in 1:R.sz){
        if(k1<k2){
          ddfact = setNames(sapply(1:nd,FUN=function(ei){
            dyn.ei = dyn[[paste0("eta_",ei)]]
            prod(sapply(c(k1,k2),FUN=function(k){
              yGk = ObsModel.transfo$linkR[k]

              sig = Rerr[[yGk]]
              tki = Robs_i[[yGk]]$time
              Zki = Robs_i[[yGk]][,yGk]

              a0 = theta$alpha0[[yGk]]
              a1 = alpha1[[yGk]]
              trs = ObsModel.transfo$R[which(ObsModel.transfo$linkR==yGk)]
              Rki = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tki,names(trs)])

              return(
                1/sig**2*(sum(Rki*(Zki-a0-a1*Rki)))
              )
            }))
          }),paste0("eta_",1:nd))

          # ddLLi[k1,k2] <- ddLLi[k2,k1] <- 1/Li*sum(mh.parm$Weights*ddfact*R.margDensity*S.margDensity)-prod(dLLi[c(k1,k2)])
          res.aux <- lapply(1:nd,function(ei){
            decomp_intermediaire = decomp(mh.parm$Weights[[ei]]*ddfact[[ei]]*R.margDensity[[ei]][["mantissa"]]*S.margDensity[[ei]])

            mantissa_res = decomp_intermediaire["mantissa"]
            exponent_res = R.margDensity[[ei]][["exponent"]] + decomp_intermediaire["exponent"]

            return(c(exponent_res,mantissa_res))
          })
          max_exponent <- max(sapply(res.aux, function(x) x[["exponent"]]))
          decomp.aux <- decomp(sum(sapply(res.aux,function(ri){
            return(ri[["mantissa"]]*10**(ri[["exponent"]] - max_exponent))
          })))

          aux =c(exponent=(max_exponent+decomp.aux[["exponent"]]),mantissa=decomp.aux[["mantissa"]])

          ddLLi[k1,k2] <- ddLLi[k2,k1] <- (invLi[["mantissa"]]*aux[["mantissa"]]*10**(invLi[["exponent"]]+aux[["exponent"]])) - prod(dLLi[c(k1,k2)])
        }
      }
    }
  }else{
    dLLi <- ddLLi <- NULL
  }
  return(list(Li=Li[["mantissa"]]*10**(Li[["exponent"]]),LLi=LLi,dLLi=dLLi,ddLLi=ddLLi))
}


# One dimensional Adaptative Hermite Gauss  ------------------------------------------
agauss.hermite <- function(n,mu=0,sd=1){
  gh <- fastGHQuad::gaussHermiteData(n)
  names(gh) <- c("Points","Weights")

  ww <- gh$Weights * sqrt(2) * sd * exp(gh$Points**2)
  xx <- mu + sd * sqrt(2) * gh$Points

  return(data.frame(Points=xx,Weights=ww))
}

# Multi-dimensional Adaptative Hermite Gauss -----------------------------------------
amgauss.hermite <- function(n,mu=rep(0,ncol(Omega)),Omega=diag(rep(1,length(mu))),prune=NULL){
  dm = length(mu)

  gh  <- as.data.frame(fastGHQuad::gaussHermiteData(n))

  idx <- as.matrix(expand.grid(rep(list(1:n),dm)))
  pts <- matrix(gh[idx,1],nrow(idx),dm)
  wts <- apply(matrix(gh[idx,2],nrow(idx),dm), 1, prod)
  # plot(pts, cex=-5/log(wts), pch=19)


  if(!is.null(prune)) {
    qwt <- quantile(wts, probs=prune)
    pts <- pts[wts > qwt,]
    wts <- wts[wts > qwt]
  }

  eig <- eigen(Omega)
  S = eig$vectors
  L = diag(sqrt(eig$values))
  rot <- S %*% L

  ppts  <-t( sqrt(2) * rot %*% t(pts) + mu)
  wwts  <- wts * sqrt(2)**dm * abs(det(rot)) * exp(apply(pts,1,FUN=function(x){t(x)%*%x}))


  return(list(Points=ppts, Weights=wwts))
}

# individual contribution to log-lik gradient in alpha1=0 ----------------------------
lambda.max.ind <- function(
    dynFUN,
    y,
    mu_i=NULL,
    Omega_i=NULL,
    theta=NULL,
    alpha1=NULL,
    covariates_i=NULL,
    ParModel.transfo=NULL,
    ParModel.transfo.inv=NULL,
    Sobs_i=NULL,
    Robs_i=NULL,
    Serr=NULL,
    Rerr=NULL,
    ObsModel.transfo=NULL,
    data=NULL,
    ind = NULL,
    n = NULL,
    prune=NULL){



  mu <- Omega <- Sobs <- Robs <- covariates <- NULL

  if(is.null(data)){
    test <- sapply(c("mu_i","Omega_i","theta","alpha1","covariates_i","ParModel.transfo","ParModel.transfo.inv","Sobs_i","Robs_i","Serr","Rerr","ObsModel.transfo"),FUN=is.null)
    if(any(test))
      stop("Please provide all necessary arguments.")
  }else{
    if((is.null(mu_i) || is.null(Omega_i)) & !all(c("mu","Omega") %in% names(data))){
      stop("Please provide mu and Omega if these are missing from data.")
    }
    argsNONind = c("theta","alpha1","ParModel.transfo","ParModel.transfo.inv","Serr","Rerr","ObsModel.transfo")
    for(d in 1:length(argsNONind)){
      test <- eval(parse(text=paste0("is.null(",argsNONind[d],")")))
      if(test){
        eval(parse(text=paste0(argsNONind[d],"<- data[[argsNONind[d]]]")))
      }
    }
    argsInd = c("mu","Omega","Robs","Sobs","covariates")
    for(d in 1:length(argsInd)){
      test <- eval(parse(text=paste0("is.null(",paste0(argsInd[d],"_i"),")")))
      if(test){
        eval(parse(text=paste0(argsInd[d],"<- data[[argsInd[d]]]")))

      }
    }
    if((is.null(mu_i) || is.null(Omega_i) || is.null(Robs_i) || is.null(Sobs_i) || is.null(covariates_i)) & is.null(ind)){
      stop("Please provide individual information (mu_i,Omega_i,Sobs_i,Robs_i,covariates_i), or the individual id ind.")
    }
    if(is.null(mu_i)){
      mu_i = mu[[ind]]
    }
    if(is.null(Omega_i)){
      Omega_i = Omega[[ind]]
    }
    if(is.null(Sobs_i)){
      Sobs_i = lapply(Sobs,FUN=function(S){S[S$id==ind,]})
    }
    if(is.null(Robs_i)){
      Robs_i = lapply(Robs,FUN=function(R){R[R$id==ind,]})
    }
    if(is.null(covariates_i)){
      covariates_i = covariates[ind,,drop=F]
    }
  }

  if(is.null(n)){
    n <- floor(100**(1/length(theta$psi_pop)))
  }

  dm = length(mu_i)

  if(dm!=1){
    mh.parm <- amgauss.hermite(n,mu=mu_i,Omega=Omega_i,prune=prune)
    detsq.omega = prod(theta$omega)
    root.omega = diag(1/theta$omega**2)
  }else{
    mh.parm <- agauss.hermite(n,mu = mu_i,sd=sqrt(as.numeric(Omega_i)))
    detsq.omega = as.numeric(theta$omega)
    root.omega = 1/as.numeric(theta$omega)**2
  }

  nd = length(mh.parm$Weights)
  R.sz = length(Rerr)
  S.sz = length(Serr)
  all.tobs = sort(union(unlist(lapply(Robs_i,FUN=function(x){x$time})),unlist(lapply(Sobs_i,FUN=function(x){x$time}))))

  dyn <- setNames(lapply(split(mh.parm$Points,1:nd),FUN=function(eta_i){
    PSI_i  = indParm(theta[c("phi_pop","psi_pop","gamma","beta")],covariates_i,setNames(eta_i,colnames(Omega_i)),ParModel.transfo,ParModel.transfo.inv)
    dyn_eta_i <- dynFUN(all.tobs,y,unlist(unname(PSI_i)))

    return(dyn_eta_i)
  }),paste0("eta_",1:nd))

  # Compute marginal latent density for each individual and value of RE
  R.margDensity <- setNames(sapply(1:nd,FUN=function(ei){
    dyn.ei = dyn[[paste0("eta_",ei)]]


    res = prod(sapply(1:R.sz,FUN=function(k){
      yGk = ObsModel.transfo$linkR[k]
      sig = Rerr[[yGk]]
      tki = Robs_i[[yGk]]$time
      Zki = Robs_i[[yGk]][,yGk]

      a0 = theta$alpha0[[yGk]]
      a1 = alpha1[[yGk]]
      trs = ObsModel.transfo$R[which(ObsModel.transfo$linkR==yGk)]
      Rki = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tki,names(trs)])

      prod(1/(sig*sqrt(2*pi))*exp(-1/2*((Zki-a0-a1*Rki)/sig)**2))
    }))
  }),paste0("eta_",1:nd))

  # A(Y|theta0)
  if(S.sz!=0){
    S.margDensity <- setNames(sapply(1:nd,FUN=function(ei){
      eta_i = split(mh.parm$Points,1:nd)[[ei]]
      dyn.ei = dyn[[paste0("eta_",ei)]]

      res = prod(sapply(1:S.sz,FUN=function(p){
        Yp = ObsModel.transfo$linkS[p]
        sig = Serr[[Yp]]
        tpi = Sobs_i[[Yp]]$time
        Ypi = Sobs_i[[Yp]][,Yp]

        trs = ObsModel.transfo$S[which(ObsModel.transfo$linkS==Yp)]
        Spi = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tpi,names(trs)])

        prod(1/(sig*sqrt(2*pi))*exp(-1/2*((Ypi-Spi)/sig)**2))
      }))*(1/((2*pi)**(dm/2)*detsq.omega)*exp(-1/2*eta_i%*%root.omega%*%eta_i))
    }),paste0("eta_",1:nd))
  }else{
    S.margDensity = setNames(sapply(split(mh.parm$Points,1:nd),FUN=function(eta_i){
      (1/((2*pi)**(dm/2)*detsq.omega)*exp(-1/2*eta_i%*%root.omega%*%eta_i))
    }),paste0("eta_",1:nd))
  }

  # s'k(Ri(krij))A(Y|theta0)
  S2.margDensity <- lapply(1:R.sz,FUN=function(k){
    setNames(sapply(1:nd,FUN=function(ei){
      dyn.ei = dyn[[paste0("eta_",ei)]]

      yGk = ObsModel.transfo$linkR[k]
      sig = Rerr[[yGk]]
      tki = Robs_i[[yGk]]$time
      Zki = Robs_i[[yGk]][,yGk]

      a0 = theta$alpha0[[yGk]]
      a1 = alpha1[[yGk]]
      trs = ObsModel.transfo$R[which(ObsModel.transfo$linkR==yGk)]
      Rki = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tki,names(trs)])

      return(1/sig**2*sum((Zki-a0)*Rki)*S.margDensity[[ei]])
    }),paste0("eta_",1:nd))})


  indi.contrib <- sapply(1:R.sz,FUN=function(k){
    sum(mh.parm$Weights*S2.margDensity[[k]])
  })/sum(mh.parm$Weights*S.margDensity)

  return(indi.contrib)
}

# log-lik gradient in alpha1=0 for maximum penalization parameter to reach ----------
lambda.max  <- function(
    dynFUN,
    y,
    mu=NULL,
    Omega=NULL,
    theta=NULL,
    alpha1=NULL,
    covariates=NULL,
    ParModel.transfo=NULL,
    ParModel.transfo.inv=NULL,
    Sobs=NULL,
    Robs=NULL,
    Serr=NULL,
    Rerr=NULL,
    ObsModel.transfo=NULL,
    data=NULL,
    n = NULL,
    prune=NULL,
    parallel = TRUE,
    ncores=NULL,
    onlyLL=FALSE,
    verbose=TRUE){

  if(is.null(data)){
    test <- sapply(c("mu","Omega","theta","alpha1","covariates","ParModel.transfo","ParModel.transfo.inv","Sobs","Robs","Serr","Rerr","ObsModel.transfo"),FUN=is.null)
    if(any(test))
      stop("Please provide all necessary arguments.")
  }else{
    if((is.null(mu) || is.null(Omega)) & !all(c("mu","Omega") %in% names(data))){
      stop("Please provide mu and Omega if these are missing from data.")
    }
    for(d in 1:length(data)){
      test <- eval(parse(text=(paste0("is.null(",names(data)[d],")"))))
      if(test){
        eval(parse(text=paste0(names(data[d]),"<- data[[d]]")))
      }
    }
  }
  if(is.null(n)){
    n <- floor(100**(1/length(theta$psi_pop)))
  }

  if(parallel){
    if(!is.null(ncores)){
      cluster <- parallel::makeCluster(ncores)
    }else{
      cluster <- parallel::makeCluster(parallel::detectCores())
    }
    doSNOW::registerDoSNOW(cluster)
    }

  N=length(mu)



  i = 1
  if(verbose){
    ntasks <- N
    pb <- txtProgressBar(max = ntasks, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  }else{
    opts <- NULL
  }

  res = foreach::foreach(i = 1:N,.packages = "REMixed",.export = c("lambda.max.ind","amgauss.hermite"),.options.snow=opts)%dopar%{
    if(0 %in% diag(Omega[[i]])){
      diag(Omega[[i]])[diag(Omega[[i]])==0] <- 10**(-5)
    }
    lambda.max.ind(mu_i = mu[[i]],
                   Omega_i = Omega[[i]],
                   theta = theta,
                   alpha1 = alpha1,
                   dynFUN = dynFUN,
                   y = y,
                   covariates_i = covariates[i,,drop=F],
                   ParModel.transfo = ParModel.transfo,
                   ParModel.transfo.inv = ParModel.transfo.inv,
                   Sobs_i = lapply(Sobs,FUN=function(S){S[S$id==i,]}),
                   Robs_i = lapply(Robs,FUN=function(R){R[R$id==i,]}),
                   Serr = Serr,
                   Rerr = Rerr,
                   ObsModel.transfo = ObsModel.transfo,
                   n = n,
                   prune = prune)
    }

  if(verbose)
    close(pb)
  if(parallel)
    snow::stopCluster(cluster)
  return(max(Reduce("+",res)))
}

decomp <- function(x){
  if(x==0){
    return(c(exponent=0,mantissa=0))
  }
  exponent_x <- floor(log10(abs(x))) # Exponent
  exponent_lim <- abs(floor(log10(.Machine$double.xmin)))

  if(abs(exponent_x)>exponent_lim){
    mantissa_x <- x*10**(-sign(exponent_x)*(exponent_lim-1))
    mantissa_x <- mantissa_x*10**(-(exponent_x+(exponent_lim-1)))
  }else{
    mantissa_x <- x*10**(-exponent_x)
  }

  return(c(exponent=exponent_x,mantissa=mantissa_x))
}
#
# gh.LL.ind_OLD <- function(
#     dynFUN,
#     y,
#     mu_i=NULL,
#     Omega_i=NULL,
#     theta=NULL,
#     alpha1=NULL,
#     covariates_i=NULL,
#     ParModel.transfo=NULL,
#     ParModel.transfo.inv=NULL,
#     Sobs_i=NULL,
#     Robs_i=NULL,
#     Serr=NULL,
#     Rerr=NULL,
#     ObsModel.transfo=NULL,
#     data=NULL,
#     ind = NULL,
#     n = NULL,
#     prune=NULL,
#     onlyLL=FALSE){
#
#   mu <- Omega <- Sobs <- Robs <- covariates <- NULL
#
#   if(is.null(data)){
#     test <- sapply(c("mu_i","Omega_i","theta","alpha1","covariates_i","ParModel.transfo","ParModel.transfo.inv","Sobs_i","Robs_i","Serr","Rerr","ObsModel.transfo"),FUN=is.null)
#     if(any(test))
#       stop("Please provide all necessary arguments.")
#   }else{
#     if((is.null(mu_i) || is.null(Omega_i)) & !all(c("mu","Omega") %in% names(data))){
#       stop("Please provide mu and Omega if these are missing from data.")
#     }
#     argsNONind = c("theta","alpha1","ParModel.transfo","ParModel.transfo.inv","Serr","Rerr","ObsModel.transfo")
#     for(d in 1:length(argsNONind)){
#       test <- eval(parse(text=paste0("is.null(",argsNONind[d],")")))
#       if(test){
#         eval(parse(text=paste0(argsNONind[d],"<- data[[argsNONind[d]]]")))
#       }
#     }
#     argsInd = c("mu","Omega","Robs","Sobs","covariates")
#     for(d in 1:length(argsInd)){
#       test <- eval(parse(text=paste0("is.null(",paste0(argsInd[d],"_i"),")")))
#       if(test){
#         eval(parse(text=paste0(argsInd[d],"<- data[[argsInd[d]]]")))
#
#       }
#     }
#     if((is.null(mu_i) || is.null(Omega_i) || is.null(Robs_i) || is.null(Sobs_i) || is.null(covariates_i)) & is.null(ind)){
#       stop("Please provide individual information (mu_i,Omega_i,Sobs_i,Robs_i,covariates_i), or the individual id ind.")
#     }
#     if(is.null(mu_i)){
#       mu_i = mu[[ind]]
#     }
#     if(is.null(Omega_i)){
#       Omega_i = Omega[[ind]]
#     }
#     if(is.null(Sobs_i)){
#       Sobs_i = lapply(Sobs,FUN=function(S){S[S$id==ind,]})
#     }
#     if(is.null(Robs_i)){
#       Robs_i = lapply(Robs,FUN=function(R){R[R$id==ind,]})
#     }
#     if(is.null(covariates_i)){
#       covariates_i = covariates[ind,,drop=F]
#     }
#   }
#
#   if(is.null(n)){
#     n <- floor(100**(1/length(theta$psi_pop)))
#   }
#
#   dm = length(mu_i)
#
#   if(dm!=1){
#     mh.parm <- amgauss.hermite(n,mu=mu_i,Omega=Omega_i,prune=prune)
#     detsq.omega = prod(theta$omega)
#     root.omega = diag(1/theta$omega**2)
#   }else{
#     mh.parm <- agauss.hermite(n,mu = mu_i,sd=sqrt(as.numeric(Omega_i)))
#     detsq.omega = as.numeric(theta$omega)
#     root.omega = 1/as.numeric(theta$omega)**2
#   }
#
#   nd = length(mh.parm$Weights)
#   R.sz = length(Rerr)
#   S.sz = length(Serr)
#   all.tobs = sort(union(unlist(lapply(Robs_i,FUN=function(x){x$time})),unlist(lapply(Sobs_i,FUN=function(x){x$time}))))
#
#   # Need to compute, for each eta_i x individual i, the dynamics of the model
#   # eta_i = split(mh.parm$Points,1:nd)[[1]]
#   dyn <- setNames(lapply(split(mh.parm$Points,1:nd),FUN=function(eta_i){
#     PSI_i  = indParm(theta[c("phi_pop","psi_pop","gamma","beta")],covariates_i,setNames(eta_i,colnames(Omega_i)),ParModel.transfo,ParModel.transfo.inv)
#     dyn_eta_i <- dynFUN(all.tobs,y,unlist(unname(PSI_i)))
#
#     return(dyn_eta_i)
#   }),paste0("eta_",1:nd))
#
#   # Compute marginal latent density for each individual and value of RE
#   R.margDensity <- setNames(sapply(1:nd,FUN=function(ei){
#     dyn.ei = dyn[[paste0("eta_",ei)]]
#
#
#     res = prod(sapply(1:R.sz,FUN=function(k){
#       yGk = ObsModel.transfo$linkR[k]
#       sig = Rerr[[yGk]]
#       tki = Robs_i[[yGk]]$time
#       Zki = Robs_i[[yGk]][,yGk]
#
#       a0 = theta$alpha0[[yGk]]
#       a1 = alpha1[[yGk]]
#       trs = ObsModel.transfo$R[which(ObsModel.transfo$linkR==yGk)]
#       Rki = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tki,names(trs)])
#
#       prod(1/(sig*sqrt(2*pi))*exp(-1/2*((Zki-a0-a1*Rki)/sig)**2))
#     }))
#   }),paste0("eta_",1:nd))
#
#   # Compute marginal density for each individual and value of RE
#   if(S.sz!=0){
#     S.margDensity <- setNames(sapply(1:nd,FUN=function(ei){
#       eta_i = split(mh.parm$Points,1:nd)[[ei]]
#       dyn.ei = dyn[[paste0("eta_",ei)]]
#
#       res = prod(sapply(1:S.sz,FUN=function(p){
#         Yp = ObsModel.transfo$linkS[p]
#         sig = Serr[[Yp]]
#         tpi = Sobs_i[[Yp]]$time
#         Ypi = Sobs_i[[Yp]][,Yp]
#
#         trs = ObsModel.transfo$S[which(ObsModel.transfo$linkS==Yp)]
#         Spi = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tpi,names(trs)])
#
#         prod(1/(sig*sqrt(2*pi))*exp(-1/2*((Ypi-Spi)/sig)**2))
#       }))*(1/((2*pi)**(dm/2)*detsq.omega)*exp(-1/2*eta_i%*%root.omega%*%eta_i))
#     }),paste0("eta_",1:nd))
#   }else{
#     S.margDensity = setNames(sapply(split(mh.parm$Points,1:nd),FUN=function(eta_i){
#       (1/((2*pi)**(dm/2)*detsq.omega)*exp(-1/2*eta_i%*%root.omega%*%eta_i))
#     }),paste0("eta_",1:nd))
#   }
#
#   # Compute individual log-Likelihood
#   Li = sum(mh.parm$Weights*R.margDensity*S.margDensity)
#   LLi = log(Li)
#
#   if(!onlyLL){
#     # Compute  gradient of individual log-Likelihood
#     # dim.alpha = length(alpha1) # = K = length(Robs_i )
#     dfact = lapply(1:R.sz,FUN=function(k){
#       setNames(sapply(1:nd,FUN=function(ei){
#         dyn.ei = dyn[[paste0("eta_",ei)]]
#         yGk = ObsModel.transfo$linkR[k]
#         sig = Rerr[[yGk]]
#         tki = Robs_i[[yGk]]$time
#         Zki = Robs_i[[yGk]][,yGk]
#
#         a0 = theta$alpha0[[yGk]]
#         a1 = alpha1[[yGk]]
#         trs = ObsModel.transfo$R[which(ObsModel.transfo$linkR==yGk)]
#         Rki = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tki,names(trs)])
#
#         return(1/(sig**2)*(sum(Rki*(Zki-a0-a1*Rki))))
#       }),paste0("eta_",1:nd))
#     })
#
#     dLLi = sapply(dfact,FUN=function(f){
#       return(1/Li*sum(mh.parm$Weights*f*R.margDensity*S.margDensity))
#     })
#
#     # Compute hessienne of individual log-Likelihood
#     ddLLi = matrix(0,ncol=R.sz,nrow=R.sz)
#     # diagonal term
#     ddfact = lapply(1:R.sz,FUN=function(k){
#       setNames(sapply(1:nd,FUN=function(ei){
#         dyn.ei = dyn[[paste0("eta_",ei)]]
#         yGk = ObsModel.transfo$linkR[k]
#         sig = Rerr[[yGk]]
#         tki = Robs_i[[yGk]]$time
#         Zki = Robs_i[[yGk]][,yGk]
#
#         a0 = theta$alpha0[[yGk]]
#         a1 = alpha1[[yGk]]
#         trs = ObsModel.transfo$R[which(ObsModel.transfo$linkR==yGk)]
#         Rki = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tki,names(trs)])
#
#         return(
#           (1/(sig**2)*(sum(Rki*(Zki-a0-a1*Rki))))**2 - 1/(sig**2)*(sum(Rki**2))
#         )
#       }),paste0("eta_",1:nd))
#     })
#
#     diag(ddLLi) = sapply(ddfact,FUN=function(f){
#       return(1/Li*sum(mh.parm$Weights*f*R.margDensity*S.margDensity))
#     })-(dLLi)**2
#
#     # non diagonal term
#     for(k1 in 1:R.sz){
#       for(k2 in 1:R.sz){
#         if(k1<k2){
#           ddfact = setNames(sapply(1:nd,FUN=function(ei){
#             dyn.ei = dyn[[paste0("eta_",ei)]]
#             prod(sapply(c(k1,k2),FUN=function(k){
#               yGk = ObsModel.transfo$linkR[k]
#
#               sig = Rerr[[yGk]]
#               tki = Robs_i[[yGk]]$time
#               Zki = Robs_i[[yGk]][,yGk]
#
#               a0 = theta$alpha0[[yGk]]
#               a1 = alpha1[[yGk]]
#               trs = ObsModel.transfo$R[which(ObsModel.transfo$linkR==yGk)]
#               Rki = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tki,names(trs)])
#
#               return(
#                 1/sig**2*(sum(Rki*(Zki-a0-a1*Rki)))
#               )
#             }))
#           }),paste0("eta_",1:nd))
#
#           ddLLi[k1,k2] <- ddLLi[k2,k1] <- 1/Li*sum(mh.parm$Weights*ddfact*R.margDensity*S.margDensity)-prod(dLLi[c(k1,k2)])
#         }
#       }
#     }
#   }else{
#     dLLi <- ddLLi <- NULL
#   }
#   return(list(Li=Li,LLi=LLi,dLLi=dLLi,ddLLi=ddLLi))
# }
