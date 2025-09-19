#' Line search
#' @noRd
searpas <- function(vw,step,b,delta,funcpa,res.out.error,stored=NULL,print=FALSE,...){

  valfpa <- function(vw,b,delta,funcpa,stored,...){
    if(is.na(vw)) return(list(fim=-Inf,stored=stored))
    bk <- b + (exp(vw)*delta)

    res = funcpa(bk,stored,...)
    funcpa.out = res$pen
    if(is.null(stored)){
      stored = res$stored
    }
    return(list(fim=-funcpa.out,stored=stored))
  }

  goto50  <- function(step,vlw2,fi1,fi2,fi3,b,delta,funcpa,stored,...){
    vm <- vlw2-(step*(fi1-fi3))/(2*(fi1-2*fi2+fi3))
    res <- valfpa(vm,b,delta,funcpa,stored,...)
    if(is.null(stored)){
      return(list(vm=vm,fim=res$fim,stored=res$stored))
    }else{
      return(list(vm=vm,fim=res$fim,stored=stored))
    }
  }
  vlw1 <- log(vw)
  vlw2 <- vlw1+step
  fi1 <- valfpa(vlw1,b,delta,funcpa,stored,...)
  stored <- fi1$stored
  fi1 <- fi1$fim
  fi2 <- valfpa(vlw2,b,delta,funcpa,stored,...)$fim

  if((sum(!is.finite(fi1)) > 0) || (sum(!is.finite(fi2)) > 0)){
    if(print){
      cat("Probably too much accuracy requested...\n")
      cat("Last step values :\n")
      cat("      b :",res.out.error$old.b,"\n")
      cat("      function value :",res.out.error$old.rl,"\n")
      cat("      Convergence criteria: parameters stability=", res.out.error$old.ca, "\n")
      cat("                          : function stability=", res.out.error$old.cb, "\n")
      cat("                          : best relative distance to maximum obtained (RDM)=", res.out.error$old.dd, "\n")
      stop("")
    }
  }

  if((fi2 >= fi1)){
    vlw3 <- vlw2
    vlw2 <- vlw1
    fi3 <- fi2
    fi2 <- fi1
    step <- -step
    vlw1 <- vlw2+step
    fi1 <- valfpa(vlw1,b,delta,funcpa,stored,...)$fim
    gt50 <- goto50(step,vlw2,fi1,fi2,fi3,b,delta,funcpa,stored,...)
    vm <- gt50$vm
    fim <- gt50$fim
    if(is.na(fim)) fim <- Inf
    if(fim <= fi2){
      vw <- exp(vm)
    }else{
      vm <- vlw2
      fim <- fi2
      vw <- exp(vm)
    }

  }else{
    vlw <- vlw1
    vlw1 <- vlw2
    vlw2 <- vlw
    fim <- fi1
    fi1 <- fi2
    fi2 <- fim

    for(i in 1:40){
      vlw3 <- vlw2
      vlw2 <- vlw1
      fi3 <- fi2
      fi2 <- fi1
      vlw1=vlw2+step
      fi1 <- valfpa(vlw1,b,delta,funcpa,stored,...)$fim
      if(fi1 > fi2){
        gt50 <- goto50(step,vlw2,fi1,fi2,fi3,b,delta,funcpa,stored,...)
        out <- 1
        break
      }
      if(fi1 == fi2){
        fim <- fi2
        vm <- vlw2
        vw <- exp(vm)
        out <- 1
        break
      }
    }
  }
  return(list(vw=vw,fim=fim))
}

funcpa <- function(b,stored,dynFUN, y, data, n, prune,to.recalibrate, ncores,parallel,lambda){

  data$alpha1[to.recalibrate] <- b
  res = fim.searpas(dynFUN = dynFUN, y = y, data = data, n = n, prune = prune, ncores = ncores,parallel=parallel,stored=stored, to.recalibrate=to.recalibrate)
  if(is.null(stored)){
    stored = res$stored
  }
  LL = res$LL
  return(list(pen=LL-lambda*sum(abs(data$alpha1)),stored=stored))
}

# gh.LL.fim ---------------------------------------------------------------
fim.searpas <- function(
    dynFUN,
    y,
    stored,
    to.recalibrate,
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
    ncores=NULL){

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
  pb <- utils::txtProgressBar(max = ntasks, style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  res = foreach::foreach(i = 1:N,.packages = "REMixed",.export = "fim.searpas.ind")%dopar%{
    if(0 %in% diag(Omega[[i]])){
      diag(Omega[[i]])[diag(Omega[[i]])==0] <- 10**(-5)
    }
    fim.searpas.ind(mu_i = mu[[i]],
              Omega_i = Omega[[i]],
              theta = theta,
              alpha1 = alpha1,
              dynFUN = dynFUN,
              y = y,
              covariates_i = covariates[i,,drop=FALSE],
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
              to.recalibrate = to.recalibrate,
              stored=stored[[i]])
  }

  close(pb)
  if(parallel)
    snow::stopCluster(cluster)

  LL = sum(sapply(res,FUN=function(ri){ri$LL}))

  if(is.null(stored)){
    stored = lapply(res,FUN=function(ri){ri$stored})
  }

  return(list(LL=LL,stored=stored))
}

fim.searpas.ind <- function(
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
    to.recalibrate,
    stored=NULL,
    precBits=10){
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
      covariates_i = covariates[ind,,drop=FALSE]
    }
  }

  if(is.null(n)){
    n <- floor(100**(1/length(theta$psi_pop)))
  }

  dm = length(mu_i)

  # STORED ----------------------------------------------------------------
  if(!is.null(stored)){
    mh.parm = stored$agh$mh.parm
    detsq.omega = stored$agh$detsq.omega
    root.omega = stored$agh$root.omega

    nd = length(mh.parm$Weights)
    R.sz = length(Rerr)
    not.computed.again = setdiff(1:R.sz,to.recalibrate)
    S.sz = length(Serr)
    all.tobs = sort(union(unlist(lapply(Robs_i,FUN=function(x){x$time})),unlist(lapply(Sobs_i,FUN=function(x){x$time}))))

    dyn <- stored$dyn


    R.margDensity <-  setNames(lapply(1:nd,FUN=function(ei){
      dyn.ei = dyn[[paste0("eta_",ei)]]

      decomp_marg = lapply(to.recalibrate,FUN=function(k){
        yGk = ObsModel.transfo$linkR[k]
        sig = Rerr[[yGk]]
        tki = Robs_i[[yGk]]$time
        Zki = Robs_i[[yGk]][,yGk]

        a0 = theta$alpha0[[yGk]]
        a1 = alpha1[[yGk]]
        trs = ObsModel.transfo$R[which(ObsModel.transfo$linkR==yGk)]
        Rki = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tki,names(trs)])

        aux <- 1/(sig*sqrt(2*pi))*exp(-1/2*((Zki-a0-a1*Rki)/sig)**2)
        if(any(aux==0)){
          pw <- Rmpfr::mpfr(-1/2*((Zki-a0-a1*Rki)/sig)**2,precBits = precBits)

          aux = 1/(sig*sqrt(2*pi))*exp(pw)

          decomp_aux = setNames(sapply(decomp(prod(aux)),as.numeric),c("exponent","mantissa"))
        }else{

          decomp_aux <- lapply(aux,decomp)

          decomp_intermediaire = decomp(prod(sapply(decomp_aux,FUN=function(decomp){decomp["mantissa"]})))

          mantissa_aux = decomp_intermediaire["mantissa"]
          exponent_aux = sum(sapply(decomp_aux,FUN=function(decomp){decomp["exponent"]})) + decomp_intermediaire["exponent"]

          decomp_aux <- c(exponent_aux,mantissa_aux)
        }
        return(decomp_aux)
      })

      if(any(sapply(decomp_marg,function(x){x["mantissa"]})==0)){
        stop("[Error] in log-likelihood computation, latent part of the marginal likelihood is equal to zero.")
      }

      # decomp_marg = lapply(little_marg,decomp)

      decomp_intermediaire = decomp(prod(sapply(decomp_marg,FUN=function(decomp){decomp["mantissa"]})))

      mantissa_res = decomp_intermediaire["mantissa"]
      exponent_res = sum(sapply(decomp_marg,FUN=function(decomp){decomp["exponent"]})) + decomp_intermediaire["exponent"]

      decomp_intermediaire2 = decomp(unname(stored$R.marg[[ei]]["mantissa"]*mantissa_res))

      mantissa_res = decomp_intermediaire2["mantissa"]
      exponent_res = exponent_res + stored$R.marg[[ei]]["exponent"] + decomp_intermediaire2["exponent"]

      return(c(exponent_res,mantissa_res))
    }),paste0("eta_",1:nd))

    S.margDensity <- stored$S.marg




  # Not STORED --------------------------------------------------------------
  }else{
    if(dm!=1){
      mh.parm <- amgauss.hermite(n,mu=mu_i,Omega=Omega_i,prune=prune)
      detsq.omega = prod(theta$omega)
      root.omega = diag(1/theta$omega**2)
    }else{
      mh.parm <- agauss.hermite(n,mu = mu_i,sd=sqrt(as.numeric(Omega_i)))
      detsq.omega = as.numeric(theta$omega)
      root.omega = 1/as.numeric(theta$omega)**2
    }

    stored <- list()
    stored$agh <-list(mh.parm = mh.parm,
                      detsq.omega = detsq.omega,
                      root.omega = root.omega)

    nd = length(mh.parm$Weights)
    R.sz = length(Rerr)
    not.computed.again = setdiff(1:R.sz,to.recalibrate)
    S.sz = length(Serr)
    all.tobs = sort(union(unlist(lapply(Robs_i,FUN=function(x){x$time})),unlist(lapply(Sobs_i,FUN=function(x){x$time}))))

    dyn <- setNames(lapply(split(mh.parm$Points,1:nd),FUN=function(eta_i){
      PSI_i  = indParm(theta[c("phi_pop","psi_pop","gamma","beta")],covariates_i,setNames(eta_i,colnames(Omega_i)),ParModel.transfo,ParModel.transfo.inv)

      yi <- sapply(y,function(yk){
        if(length(yk)==1){
          return(yk)
        }else{
          return(yk[[ind]])
        }
      })

      dyn_eta_i <- dynFUN(all.tobs,yi,unlist(unname(PSI_i)))

      return(dyn_eta_i)
    }),paste0("eta_",1:nd))

    stored$dyn <- dyn

    stored$R.marg <- setNames(lapply(1:nd,FUN=function(ei){
      dyn.ei = dyn[[paste0("eta_",ei)]]

      # .Machine$double.xmin
      decomp_marg = lapply(not.computed.again,FUN=function(k){
        yGk = ObsModel.transfo$linkR[k]
        sig = Rerr[[yGk]]
        tki = Robs_i[[yGk]]$time
        Zki = Robs_i[[yGk]][,yGk]

        a0 = theta$alpha0[[yGk]]
        a1 = alpha1[[yGk]]
        trs = ObsModel.transfo$R[which(ObsModel.transfo$linkR==yGk)]
        Rki = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tki,names(trs)])

        aux <- 1/(sig*sqrt(2*pi))*exp(-1/2*((Zki-a0-a1*Rki)/sig)**2)
        if(any(aux==0)){
          pw <- Rmpfr::mpfr(-1/2*((Zki-a0-a1*Rki)/sig)**2,precBits = precBits)

          aux = 1/(sig*sqrt(2*pi))*exp(pw)

          decomp_aux = setNames(sapply(decomp(prod(aux)),as.numeric),c("exponent","mantissa"))
        }else{

          decomp_aux <- lapply(aux,decomp)

          decomp_intermediaire = decomp(prod(sapply(decomp_aux,FUN=function(decomp){decomp["mantissa"]})))

          mantissa_aux = decomp_intermediaire["mantissa"]
          exponent_aux = sum(sapply(decomp_aux,FUN=function(decomp){decomp["exponent"]})) + decomp_intermediaire["exponent"]

          decomp_aux <- c(exponent_aux,mantissa_aux)
        }
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

    R.margDensity <- setNames(lapply(1:nd,FUN=function(ei){
      dyn.ei = dyn[[paste0("eta_",ei)]]

      decomp_marg = lapply(to.recalibrate,FUN=function(k){
        yGk = ObsModel.transfo$linkR[k]
        sig = Rerr[[yGk]]
        tki = Robs_i[[yGk]]$time
        Zki = Robs_i[[yGk]][,yGk]

        a0 = theta$alpha0[[yGk]]
        a1 = alpha1[[yGk]]
        trs = ObsModel.transfo$R[which(ObsModel.transfo$linkR==yGk)]
        Rki = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tki,names(trs)])

        aux <- 1/(sig*sqrt(2*pi))*exp(-1/2*((Zki-a0-a1*Rki)/sig)**2)
        if(any(aux==0)){
          pw <- Rmpfr::mpfr(-1/2*((Zki-a0-a1*Rki)/sig)**2,precBits = precBits)

          aux = 1/(sig*sqrt(2*pi))*exp(pw)

          decomp_aux = setNames(sapply(decomp(prod(aux)),as.numeric),c("exponent","mantissa"))
        }else{

          decomp_aux <- lapply(aux,decomp)

          decomp_intermediaire = decomp(prod(sapply(decomp_aux,FUN=function(decomp){decomp["mantissa"]})))

          mantissa_aux = decomp_intermediaire["mantissa"]
          exponent_aux = sum(sapply(decomp_aux,FUN=function(decomp){decomp["exponent"]})) + decomp_intermediaire["exponent"]

          decomp_aux <- c(exponent_aux,mantissa_aux)
        }
        return(decomp_aux)
      })

      if(any(sapply(decomp_marg,function(x){x["mantissa"]})==0)){
        stop("[Error] in log-likelihood computation, latent part of the marginal likelihood is equal to zero.")
      }

      # decomp_marg = lapply(little_marg,decomp)

      decomp_intermediaire = decomp(prod(sapply(decomp_marg,FUN=function(decomp){decomp["mantissa"]})))

      mantissa_res = decomp_intermediaire["mantissa"]
      exponent_res = sum(sapply(decomp_marg,FUN=function(decomp){decomp["exponent"]})) + decomp_intermediaire["exponent"]

      decomp_intermediaire2 = decomp(unname(stored$R.marg[[ei]]["mantissa"]*mantissa_res))

      mantissa_res = decomp_intermediaire2["mantissa"]
      exponent_res = exponent_res + stored$R.marg[[ei]]["exponent"] + decomp_intermediaire2["exponent"]

      return(c(exponent_res,mantissa_res))
    }),paste0("eta_",1:nd))

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

    stored$S.marg <- S.margDensity
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

  return(list(LL=LLi,stored=stored))
}


