#' Adaptive Gauss-Hermite approximation of log-likelihood derivatives
#'
#' @description
#' Computes Adaptive Gauss-Hermite approximation of the log-likelihood and its derivatives in NLMEM with latent observation processes, see \code{\link{REMix-package}} for details on the model.
#'
#' @details
#' Based on notation introduced \code{\link{REMix-package}}. The log-likelihood of the model \eqn{LL(\theta,\alpha_1)} for a set of population parameters \eqn{\theta} and regulatization parameters \eqn{\alpha_1} is estimated using Adaptative Gausse-Hermite quadrature, using  conditional distribution estimation to locate the mass of the integrand. If the project has been initialized as a Monolix project, the user can use \code{\link{readMLX}} function to retrieve all the project information needed here.
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
#' @param n number of points to use for the Gauss-Hermite quadrature rule (see \code{\link{gaussHermiteData}}).
#' @param prune integer between 0 and 1, percentage of pruning for the Gauss-Hermite quadrature rule (default NULL).
#' @param parallel logical, if computation should be done in parallel.
#' @param ncores number of cores to use for parallelization, default will detect the number of cores available (see \code{\link{detectCores}}).
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
#'
# ---------- CALCUL DES L, LL, DLL, DDLL  ----------

# ---------- 1. Calcul des Li, LLi, dLLi, ddLLi ----------

gh.LL.ind_multi <- function(
    dynFUN,
    y,
    mu_i=NULL,
    Omega_i=NULL,
    theta=NULL,
    alphas=NULL,
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
  #data <- sortie de readMLX qui contient tout ce qu'il faut
  if(is.null(data)){ #si data est nul on verifie qu'on a quand meme tous les arguments
    test <- sapply(c("mu_i","Omega_i","theta","alphas","covariates_i","ParModel.transfo","ParModel.transfo.inv","Sobs_i","Robs_i","Serr","Rerr","ObsModel.transfo"),FUN=is.null)
    if(any(test)) #si certains sont nuls
      stop("Please provide all necessary arguments to gh.LL.ind.") #no demande tous les arguments
  }
  else{ #si data est non nul
    #si mu_i ou Omega_i sont nuls et si mu et omega ne sont pas presents dans data
        #(si mu et Omega ne sont pas dans data, all(c("mu","Omega") %in% names(data)) --> FALSE
        #donc, !all(c("mu","Omega") %in% names(data)) --> TRUE)
    if((is.null(mu_i) || is.null(Omega_i)) & !all(c("mu","Omega") %in% names(data))){
      stop("Please provide mu and Omega if these are missing from data.")
    }
    #arguments non individuels
    argsNONind = c("theta","alphas","ParModel.transfo","ParModel.transfo.inv","Serr","Rerr","ObsModel.transfo")
    #pour chacun des arguments,
    for(d in 1:length(argsNONind)){
      #on evalue si l'argument est nul ou non
        #eval(...) execute l'expression
        #parse(text = ...) transforme le texte ... en expression R
      test <- eval(parse(text=paste0("is.null(",argsNONind[d],")")))
      #si certains arguments sont nuls
      if(test){
        #on recupere l'info dans data, qui n'est autre que l'execution de readMLX
        eval(parse(text=paste0(argsNONind[d],"<- data[[argsNONind[d]]]")))
      }
    }
    #arguments individuels
    argsInd = c("mu","Omega","Robs","Sobs","covariates")
    #pour chacun d'eux
    for(d in 1:length(argsInd)){
      #on regarde s'ils sont nuls
      test <- eval(parse(text=paste0("is.null(",paste0(argsInd[d],"_i"),")")))
      #si oui, on recupere leurs valeurs dans data
      if(test){
        eval(parse(text=paste0(argsInd[d],"<- data[[argsInd[d]]]")))

      }
    }
    #si pour un individu, un des argsInd est nul et si l'entree ind est nulle aussi
    if((is.null(mu_i) || is.null(Omega_i) || is.null(Robs_i) || is.null(Sobs_i) || is.null(covariates_i)) & is.null(ind)){
      stop("Please provide individual information (mu_i,Omega_i,Sobs_i,Robs_i,covariates_i), or the individual id ind.")
    }
    #si mu_i est nul mais ind ne l'est pas :
    if(is.null(mu_i)){
      mu_i = mu[[ind]] #on recupere la bonne valeur
    }
    if(is.null(Omega_i)){
      Omega_i = Omega[[ind]]
    }
    if(is.null(Sobs_i)){
      Sobs_i = lapply(Sobs,FUN=function(S){S[S$id==ind,]}) #on recupere toutes les observations de l'individu ind
    }
    if(is.null(Robs_i)){
      Robs_i = lapply(Robs,FUN=function(R){R[R$id==ind,]})
    }
    if(is.null(covariates_i)){
      covariates_i = covariates[ind,,drop=F]
    }
  }
  #n : nb de point pour la quadrature de GH
  #s'il n'est pas donne, on choisit arbitrairement (modifie dans la derniere version d'Auriane)
  if(is.null(n)){
    if(length(theta$psi_pop)==1){
      n <- 100
    }else if(length(theta$psi_pop)==2){
      n <- 10
    }else{
      n <- 7
    }
  }
  #dimension : taille du vecteur moyenne, autant qu'il y a de parametres avec variabilite
  dm = length(mu_i)

  #Omega : sd individuelle (deviation de la valeur individuelle par rapport a la vraie valeur)
  #omega : sd globale (deviation globale par rapport a la vraie valeur)

  #si on est en multidimensionnel
  if(dm!=1){
    mh.parm <- amgauss.hermite(n,mu=mu_i,Omega=Omega_i,prune=prune) #on recupere les points et les poids
    detsq.omega = prod(theta$omega) #produit de tous els elements du vecteur
    root.omega = diag(1/theta$omega**2) #omega^-1
  }
  #si on a un seul parametre avec variabilite
  else{
    mh.parm <- agauss.hermite(n,mu = mu_i,sd=sqrt(as.numeric(Omega_i)))
    detsq.omega = as.numeric(theta$omega)
    root.omega = 1/as.numeric(theta$omega)**2
  }

  #nombre de poids
  nd = length(mh.parm$Weights)

  #Rerr : vecteur des sd des observation
      #ex : data <- readMLX
      #data$Rerr :
      # yyG1      yyG2      yyG3      yyG4      yyG5
    # 0.2753227 0.3351119 0.3877975 0.3567795 0.3192995
  R.sz = length(Rerr)
  S.sz = length(Serr)
  #on extrait les temps d'observation de Robs et Sobs
  #on prend tous les temps recuperes et on les reorganise
  all.tobs = sort(union(unlist(lapply(Robs_i,FUN=function(x){x$time})),unlist(lapply(Sobs_i,FUN=function(x){x$time}))))

  # Need to compute, for each eta_i x individual i, the dynamics of the model
  #tout est deterministe a eta_i donne

  #indParm : fonction qui genere les param indiv ayant les covariable 'covariates" et les effets aleatoires 'eta_i'
  dyn <- setNames(lapply(split(mh.parm$Points,1:nd),FUN=function(eta_i){
    #on calcule les parametres indiv de maniere deterministe
    PSI_i  = indParm(theta[c("phi_pop","psi_pop","gamma","beta")],covariates_i,setNames(eta_i,colnames(Omega_i)),ParModel.transfo,ParModel.transfo.inv)
   #on calcule la dynamique deterministe avec tous les nouveaux parametres
    #all.tobs : temps d'observation
    #y : points de depart
    #unlist(unname(PSI_i)) : valeurs des parametres
    dyn_eta_i <- dynFUN(all.tobs,y,unlist(unname(PSI_i)))
    return(dyn_eta_i)
  }),paste0("eta_",1:nd))
  #dyn est une liste nommee (eta_1, ..., eta_nd) de dynamiques

  # Compute marginal latent density for each individual and value of Random Effect
  R.margDensity <- setNames(lapply(1:nd,FUN=function(ei){
    #on recupere les points qui correspondent au bon eta_i
    dyn.ei = dyn[[paste0("eta_",ei)]]

    #pour chaque biomarqueur (R.sz = length(Rerr)) :
    decomp_marg = sapply(1:R.sz,FUN=function(k){
      yGk = ObsModel.transfo$linkR[k] #on recupere le nom du biomarqueur
      sig = Rerr[[yGk]] #on recupere son ecart type
      tki = Robs_i[[yGk]]$time #on recupere les temps d'observations pour ce biomarqueur
      Zki = Robs_i[[yGk]][,yGk] #on recupere les valeurs observees pour ce biomarqueur

      a0 = theta$alpha0[[yGk]]
      #on recupere la valeur alpha0k

      #on recupere une liste NOMMEE (alpha1, alpha2...) avec tous les alpha_l pour le k donne
      atot = get_all_alpha_l_for_one_k(alphas, yGk, L)


      #on recupere la liste NOMMEE (R1, R2...) des transformations pour un biomarqueur donne
      trs = get_all_transfo_l_for_one_k(ObsModel.transfo, yGk, L)

      #names(trs) = "S" "L" ...(remplis dans ObsModel.transfo)
      #dyn.ei[,"time"] %in% tki : on filtre les lignes pour lesquels les temps correspondent aux temps d'observation des biomarqueurs
      #dyn.ei[..., names(trs)] :
      #trs[[1]](...) : on applique la premiere transformation

      #on obtient une liste numerotee (pour le compartiment 1, [[1]]...)
      Rki = lapply(names(trs), function(nom) {
        #on applique la transformation du compartiment trs[[nom]] (S, L, ...)
        #a la colonne de donnees pour le dit compartiment
        trs[[nom]](dyn.ei[dyn.ei[,"time"] %in% tki, nom])
      })
      sum<-0
      for (l in 1:length(atot)){
      alphal<-paste0("alpha",l)
      sum<-sum + atot[[alphal]]*Rki[[l]]
      }
      #pour chaque effet aleatoire ei,
        #pour chaque biomarqueur k,
          #on a une valeur de Rki
      aux <-  prod(1/(sig*sqrt(2*pi))*exp(-1/2*((Zki-a0-sum)/sig)**2))
      if(any(aux==0)){
        pw <- Rmpfr::mpfr(-1/2*((Zki-a0-sum)/sig)**2,precBits = 32)

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

    decomp_intermediaire = decomp(prod(sapply(decomp_marg,FUN=function(decomp){decomp["mantissa"]})))

    mantissa_res = decomp_intermediaire["mantissa"]
    exponent_res = sum(sapply(decomp_marg,FUN=function(decomp){decomp["exponent"]})) + decomp_intermediaire["exponent"]

    return(c(exponent_res,mantissa_res))
  }),paste0("eta_",1:nd))
  #a ce stade on a un exposant et une mantisse pour les densites marginales de chaque ei

  # Compute marginal density for each individual and value of Random Effects
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

        #on suppose qu'on a qu'une seule transformation dans S
        Spi = trs[[1]](dyn.ei[dyn.ei[,"time"] %in% tpi,names(trs)])
        #si c'est NULL, regarder les noms de dynFUN et de ObsModel.transfo
        #fY x phi
        prod(1/(sig*sqrt(2*pi))*exp(-1/2*((Ypi-Spi)/sig)**2))
      })*(1/((2*pi)**(dm/2)*detsq.omega)*exp(-1/2*eta_i%*%root.omega%*%eta_i)))
    }),paste0("eta_",1:nd))
    # print(S.margDensity)
    # eta_1        eta_2        eta_3        eta_4        eta_5        eta_6        eta_7        eta_8        eta_9       eta_10
    # 2.859058e+03 3.700500e+03 3.713187e-06 4.808526e-06 1.702137e-03 2.202174e-03 1.316331e+03 1.703761e+03 2.369124e-03 3.084847e-03
    #on a une valeur par effet aleatoire
    }else{
      #juste phi
    S.margDensity = setNames(sapply(split(mh.parm$Points,1:nd),FUN=function(eta_i){
      (1/((2*pi)**(dm/2)*detsq.omega)*exp(-1/2*eta_i%*%root.omega%*%eta_i))
    }),paste0("eta_",1:nd))
  }

  # Compute individual log-Likelihood
  Li.aux = lapply(1:nd,function(ei){
    #pour chaque effet aleatoire ei
    #on recupere une decomposition expo-man
    decomp_intermediaire = decomp(mh.parm$Weights[[ei]]*R.margDensity[[ei]][["mantissa"]]*S.margDensity[[ei]])
    #on recupere les mantisses
    mantissa_res = decomp_intermediaire["mantissa"]
    #on recupere les exposants (somme car c'etait un produit)
    exponent_res = R.margDensity[[ei]][["exponent"]] + decomp_intermediaire["exponent"]
    #on a la decomposition du produit entre les poids fG*fY*phi
    return(c(exponent_res,mantissa_res))
  })

  #on trouve le plus grand exposant dans Li.aux
  max_exponent <- max(sapply(Li.aux, function(x) x[["exponent"]]))
  #on divise par 10**le plus grand exposant
  #on somme tous les termes de Li.aux 'renormalises"
  #--> on obtient la vraisemblance individuelle a l'echelle max_exponent pres
  #on refait la decomposition mantisse-exposant
  decomp.aux <- decomp(sum(sapply(Li.aux,function(li){
    return(li[["mantissa"]]*10**(li[["exponent"]] - max_exponent))

  })))
  #on recupere la vraisemblance en exponent-mantisse
  Li = c(exponent=(max_exponent+decomp.aux[["exponent"]]),mantissa=decomp.aux[["mantissa"]])
  #inverse de la vraisemblance :
  ##mantisse
  invLi.aux = decomp(1/Li[["mantissa"]])
  ##total
  invLi = c(exponent=-Li[["exponent"]]+invLi.aux[["exponent"]],mantissa=invLi.aux[["mantissa"]])
  #log-vraisemblance :
  # ln(p*10**o)=ln(p)+ln(10**o)=ln(p)+o*ln(10)
  LLi = log(Li[["mantissa"]])+Li[["exponent"]]*log(10)
  #print('LLi oke')

  #si on veut autre chose que juste la log-vraisemblance :
  if(!onlyLL){
    # Compute  gradient of individual log-Likelihood
    # dim.alpha = length(alpha1) # = K = length(Robs_i )
    dfact = lapply(1:R.sz,FUN=function(k){
      lapply(1:L,FUN=function(comp){
        f = setNames(sapply(1:nd,FUN=function(ei){
          dyn.ei = dyn[[paste0("eta_",ei)]]
          yGk = ObsModel.transfo$linkR[k]
          sig = Rerr[[yGk]]
          tki = Robs_i[[yGk]]$time
          Zki = Robs_i[[yGk]][,yGk]

          a0 = theta$alpha0[[yGk]]

          #rappel : alphas = (alpha11, alpha12, alpha13, alpha21, alpha22, alpha23, ...)
          #pour k=1, on veut alpha11 et alpha 21
          #on recupere une liste NOMMEE (alpha1, alpha2...) avec tous les alpha_l pour le k donne
          atot = get_all_alpha_l_for_one_k(alphas, yGk, L)
          #on recupere la liste NOMMEE (R1, R2...) des transformations pour un biomarqueur donne
          trs = get_all_transfo_l_for_one_k(ObsModel.transfo, yGk, L)
          Rki = lapply(names(trs), function(nom) {
            #on applique la transformation du compartiment trs[[nom]] (S, L, ...)
            #a la colonne de donnees pour le dit compartiment
          trs[[nom]](dyn.ei[dyn.ei[,"time"] %in% tki, nom])
         })
          somme<-0
          for (l in 1:length(atot)){
            alphal<-paste0("alpha",l)
            somme<-somme + atot[[alphal]]*Rki[[l]]
         }

        return(1/(sig**2)*(sum(Rki[[comp]]*(Zki-a0-somme))))
      }),paste0("eta_",1:nd))
    })})
    #on recupere les derivees premires selon alpha_l1_k1
    dLLi = lapply(1:R.sz, function(k){
              lapply(1:L, function(l){
                 res.aux <- lapply(1:nd,function(ei){
                  #pour chaque biomarqueur, pour chaque compartiment, on calcule le terme de l'integrale
                 decomp_intermediaire = decomp(mh.parm$Weights[[ei]]*dfact[[k]][[l]][[ei]]*R.margDensity[[ei]][["mantissa"]]*S.margDensity[[ei]])

                 mantissa_res = decomp_intermediaire["mantissa"]
                 exponent_res = R.margDensity[[ei]][["exponent"]] + decomp_intermediaire["exponent"]

        return(c(exponent_res,mantissa_res))
      })
      #on refait l'astuce de renormaliser les exposants
      #on calcule l'integrale (GH)
      max_exponent <- max(sapply(res.aux, function(x) x[["exponent"]]))
      decomp.aux <- decomp(sum(sapply(res.aux,function(ri){
        return(ri[["mantissa"]]*10**(ri[["exponent"]] - max_exponent))
      })))
      #on recupere la version decomposee de l'integrale
      aux =c(exponent=(max_exponent+decomp.aux[["exponent"]]),mantissa=decomp.aux[["mantissa"]])
      #on calcule la derivee :
      return(invLi[["mantissa"]]*aux[["mantissa"]]*10**(invLi[["exponent"]]+aux[["exponent"]]))
      # 1/Li*sum(mh.parm$Weights*f*R.margDensity*S.margDensity))
    })})
    #print('dLLi oke')

    # Compute hessienne of individual log-Likelihood
    #taille de la matrice : L x K
    ddLLi = matrix(0,ncol=R.sz*L,nrow=R.sz*L)
    # diagonal term
    ddfact = lapply(1:R.sz,FUN=function(k){
      lapply(1:L,FUN=function(comp){
      setNames(sapply(1:nd,FUN=function(ei){
        dyn.ei = dyn[[paste0("eta_",ei)]]
        yGk = ObsModel.transfo$linkR[k]
        sig = Rerr[[yGk]]
        tki = Robs_i[[yGk]]$time
        Zki = Robs_i[[yGk]][,yGk]

        a0 = theta$alpha0[[yGk]]
        atot = get_all_alpha_l_for_one_k(alphas, yGk, L)
        #on recupere la liste NOMMEE (R1, R2...) des transformations pour un biomarqueur donne
        trs = get_all_transfo_l_for_one_k(ObsModel.transfo, yGk, L)
        Rki = lapply(names(trs), function(nom) {
          #on applique la transformation du compartiment trs[[nom]] (S, L, ...)
          #a la colonne de donnees pour le dit compartiment
          trs[[nom]](dyn.ei[dyn.ei[,"time"] %in% tki, nom])
        })
        somme<-0
        for (l in 1:length(atot)){
          alphal<-paste0("alpha",l)
          somme<-somme + atot[[alphal]]*Rki[[l]]
        }
        return(
          (1/(sig**2)*(sum(Rki[[comp]]*(Zki-a0-somme))))**2 - 1/(sig**2)*(sum(Rki[[comp]]**2))
        )
      }),paste0("eta_",1:nd))
    })})


    termes_diag<-lapply(1:R.sz, function(k){
      lapply(1:L, function(l){
        res.aux <- lapply(1:nd,function(ei){
          #pour chaque biomarqueur, pour chaque compartiment, on calcule le terme de l'integrale
          decomp_intermediaire = decomp(mh.parm$Weights[[ei]]*ddfact[[k]][[l]][[ei]]*R.margDensity[[ei]][["mantissa"]]*S.margDensity[[ei]])

        mantissa_res = decomp_intermediaire["mantissa"]
        exponent_res = R.margDensity[[ei]][["exponent"]] + decomp_intermediaire["exponent"]

        return(c(exponent_res,mantissa_res))
      })

      max_exponent <- max(sapply(res.aux, function(x) x[["exponent"]]))
      decomp.aux <- decomp(sum(sapply(res.aux,function(ri){
        return(ri[["mantissa"]]*10**(ri[["exponent"]] - max_exponent))
      })))

      aux =c(exponent=(max_exponent+decomp.aux[["exponent"]]),mantissa=decomp.aux[["mantissa"]])

      return(invLi[["mantissa"]]*aux[["mantissa"]]*10**(invLi[["exponent"]]+aux[["exponent"]])-dLLi[[k]][[l]]**2)
      # return(1/Li*sum(mh.parm$Weights*f*R.margDensity*S.margDensity))
    })})
    #ici on a les termes diagonaux mais sous la forme suivante :
    # [[1]]
    # [[1]][[1]]
    # [1] -287.0224
    #
    # [[1]][[2]]
    # [1] -0.0728066
    #
    #
    # [[2]]
    # [[2]][[1]]
    # [1] -688.7226
    #
    # [[2]][[2]]
    # [1] -0.179974
    for (k in 1:R.sz) {
      for (l in 1:L) {
        # Position dans la Hessienne
        indice <- (l - 1) * R.sz + k
        ddLLi[indice, indice] <- termes_diag[[k]][[l]]
      }
    }
  #print('ddLLi diag oke')
    #sous-diagonale : meme k, different l
   # if (L>1){
    sddfact <- lapply(1:R.sz, function(k) {
      lapply(1:L, function(l1) {
        lapply(1:L, function(l2) {
          if (l1 == l2) return(NULL)

          setNames(sapply(1:nd, function(ei) {
            dyn.ei <- dyn[[paste0("eta_", ei)]]
            yGk <- ObsModel.transfo$linkR[k]
            sig <- Rerr[[yGk]]
            tki <- Robs_i[[yGk]]$time
            Zki <- Robs_i[[yGk]][, yGk]

            a0 <- theta$alpha0[[yGk]]
            atot <- get_all_alpha_l_for_one_k(alphas, yGk, L)
            trs <- get_all_transfo_l_for_one_k(ObsModel.transfo, yGk, L)

            Rki <- lapply(names(trs), function(nom) {
              trs[[nom]](dyn.ei[dyn.ei[, "time"] %in% tki, nom])
            })

            somme <- 0
            for (l in 1:L) {
              alphal <- paste0("alpha", l)
              somme <- somme + atot[[alphal]] * Rki[[l]]
            }

            s1 <- sum(Rki[[l1]] * (Zki - a0 - somme))
            s2 <- sum(Rki[[l2]] * (Zki - a0 - somme))

            return((1 / sig^4) * s1 * s2 - sum(Rki[[l1]] * Rki[[l2]]) / sig^2)
          }), paste0("eta_", 1:nd))
        })
      })
    })

    termes_sous_diag <- lapply(1:R.sz, function(k) {
      lapply(1:L, function(l1) {
        lapply(1:L, function(l2) {
          if (l1 == l2) return(NULL)

          res.aux <- lapply(1:nd, function(ei) {
            decomp_intermediaire <- decomp(mh.parm$Weights[[ei]] * sddfact[[k]][[l1]][[l2]][[ei]] *R.margDensity[[ei]][["mantissa"]] *S.margDensity[[ei]])
            mantissa_res <- decomp_intermediaire["mantissa"]
            exponent_res <- R.margDensity[[ei]][["exponent"]] + decomp_intermediaire["exponent"]
            return(c(exponent_res, mantissa_res))
          })

          max_exp <- max(sapply(res.aux, function(x) x[["exponent"]]))
          decomp.aux <- decomp(sum(sapply(res.aux, function(ri) {
            ri[["mantissa"]] * 10^(ri[["exponent"]] - max_exp)
          })))

          aux <- c(
            exponent = max_exp + decomp.aux[["exponent"]],
            mantissa = decomp.aux[["mantissa"]]
          )

          return(
            invLi[["mantissa"]] * aux[["mantissa"]] * 10^(invLi[["exponent"]] + aux[["exponent"]])
            - dLLi[[k]][[l1]] * dLLi[[k]][[l2]]
          )
        })
      })
    })
    for (k in 1:R.sz) {
      for (l1 in 1:L) {
        for (l2 in 1:L) {
          if (l1 == l2) next
          i <- (l1 - 1) * R.sz + k
          j <- (l2 - 1) * R.sz + k
          ddLLi[i, j] <- termes_sous_diag[[k]][[l1]][[l2]]
        }
      }
    }
    print('ddLLi sous diag oke')
    #le reste
    reste_ddfact <- lapply(1:R.sz, function(k1) {
      lapply(1:R.sz, function(k2) {
        lapply(1:L, function(l1) {
          lapply(1:L, function(l2) {

            #on oublie les cas deja remplis ie
            if ((k1 == k2 && l1 == l2) ||               #la diagonale
                (k1 == k2 && l1 == l2 - 1) ||           #la sous-diagonale
                (k1 > k2 || (k1 == k2 && l1 > l2))) {   #la partie triangulaire superieure
              return(NULL)
            }

            setNames(lapply(1:nd, function(ei) {
              dyn.ei <- dyn[[paste0("eta_", ei)]]
              yGk1 <- ObsModel.transfo$linkR[k1]
              sig1 <- Rerr[[yGk1]]
              tki1 <- Robs_i[[yGk1]]$time
              Zki1 <- Robs_i[[yGk1]][, yGk1]

              a01 <- theta$alpha0[[yGk1]]
              atot1 <- get_all_alpha_l_for_one_k(alphas, yGk1, L)
              trs1 <- get_all_transfo_l_for_one_k(ObsModel.transfo, yGk1, L)

              Rki1 <- lapply(names(trs1), function(nom) {
                trs1[[nom]](dyn.ei[dyn.ei[, "time"] %in% tki1, nom])
              })

              somme1 <- 0
              for (l in 1:L) {
                alphal <- paste0("alpha", l)
                somme1 <- somme1 + atot1[[alphal]] * Rki1[[l]]
              }

              yGk2 <- ObsModel.transfo$linkR[k2]
              sig2 <- Rerr[[yGk2]]
              tki2 <- Robs_i[[yGk2]]$time
              Zki2 <- Robs_i[[yGk2]][, yGk2]

              a02 <- theta$alpha0[[yGk2]]
              atot2 <- get_all_alpha_l_for_one_k(alphas, yGk2, L)
              trs2 <- get_all_transfo_l_for_one_k(ObsModel.transfo, yGk2, L)

              Rki2 <- lapply(names(trs2), function(nom) {
                trs2[[nom]](dyn.ei[dyn.ei[, "time"] %in% tki2, nom])
              })

              somme2 <- 0
              for (l in 1:L) {
                alphal <- paste0("alpha", l)
                somme2 <- somme2 + atot2[[alphal]] * Rki2[[l]]
              }

              s1 <- (1 / sig1^2) * sum(Rki1[[l1]] * (Zki1 - a01 - somme1))
              s2 <- (1 / sig2^2) * sum(Rki2[[l2]] * (Zki2 - a02 - somme2))

              return(s1 * s2)
            }), paste0("eta_", 1:nd))

          })
        })
      })
    })


    termes_restants <- lapply(1:R.sz, function(k1) {
      lapply(1:R.sz, function(k2) {
        lapply(1:L, function(l1) {
          lapply(1:L, function(l2) {

            bloc <- reste_ddfact[[k1]][[k2]][[l1]][[l2]]
            if (is.null(bloc)) return(NULL)

            res.aux <- lapply(1:nd, function(ei) {
              valeur <- bloc[[ei]]
              decomp_intermediaire <- decomp(mh.parm$Weights[[ei]] * valeur * R.margDensity[[ei]][["mantissa"]] * S.margDensity[[ei]])
              mantissa_res <- decomp_intermediaire["mantissa"]
              exponent_res <- R.margDensity[[ei]][["exponent"]] + decomp_intermediaire["exponent"]
              return(c(exponent_res, mantissa_res))
            })

            max_exp <- max(sapply(res.aux, function(x) x[["exponent"]]))
            decomp.aux <- decomp(sum(sapply(res.aux, function(ri) {
              ri[["mantissa"]] * 10^(ri[["exponent"]] - max_exp)
            })))

            aux <- c(
              exponent = max_exp + decomp.aux[["exponent"]],
              mantissa = decomp.aux[["mantissa"]]
            )

            return(
              invLi[["mantissa"]] * aux[["mantissa"]] * 10^(invLi[["exponent"]] + aux[["exponent"]])
              - dLLi[[k1]][[l1]] * dLLi[[k2]][[l2]]
            )
          })
        })
      })
    })


    for (k1 in 1:R.sz) {
      for (k2 in 1:R.sz) {
      for (l1 in 1:L) {
        for (l2 in 1:L) {
          if ((k1 == k2 && l1 == l2) ||               #la diagonale
              (k1 == k2 && l1 == l2 - 1) ||           #la sous-diagonale
              (k1 > k2 || (k1 == k2 && l1 > l2)))    #la partie triangulaire superieure
           next
          i <- (l1 - 1) * R.sz + k1
          j <- (l2 - 1) * R.sz + k2
          ddLLi[i, j] <- termes_restants[[k1]][[k2]][[l1]][[l2]]
          #symetrie
          ddLLi[j, i] <- termes_restants[[k1]][[k2]][[l1]][[l2]]
        }
      }
      }
    }
    print('ddLLi oke')
  }
    #si on a pas onlyLL
    else{
    dLLi <- ddLLi <- NULL
  }
  return(list(Li=Li[["mantissa"]]*10**(Li[["exponent"]]),LLi=LLi,dLLi=unlist(lapply(dLLi, function(dk) unlist(dk))),ddLLi=ddLLi))
}

# ---------- 2. Calcul global ----------

gh.LL_multi <- function(
    dynFUN,
    y,
    mu=NULL,
    Omega=NULL,
    theta=NULL,
    alphas=NULL,
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

  #si data est nul
  if(is.null(data)){
    #on cherche si un de ces arguments est nul
    test <- sapply(c("mu","Omega","theta","alphas","covariates","ParModel.transfo","ParModel.transfo.inv","Sobs","Robs","Serr","Rerr","ObsModel.transfo"),FUN=is.null)
    if(any(test))
      stop("Please provide all necessary arguments.")
    #si data n'est pas nul
  }else{
    #si mu, Omega sont nuls ou absents dans data
    if((is.null(mu) || is.null(Omega)) & !all(c("mu","Omega") %in% names(data))){
      stop("Please provide mu and Omega if these are missing from data.")
    }
    #si certaines variables sont definies mais nulles
    #on leur donne la bonne valeur de data
    for(d in 1:length(data)){
      test <- eval(parse(text=(paste0("is.null(",names(data)[d],")"))))
      if(test){
        eval(parse(text=paste0(names(data[d]),"<- data[[d]]")))
      }
    }
  }
  #si n, le nb de points de la quadrature, est nul
  if(is.null(n)){
    if(length(theta$psi_pop)==1){
      n <- 100
    }else if(length(theta$psi_pop)==2){
      n <- 10
    }else{
      n <- 7
    }
  }

  #si on veut faire les calcules en paralelle
  if(parallel){
    #si lutilisateur a precise le nb de coeurs a utiliser
    if(!is.null(ncores)){
      cluster <- snow::makeCluster(ncores)
      #sinon on detecte automatiquement le nb total de coeurs dispo et on les attirbue tous
    }else{
      cluster <- snow::makeCluster(parallel::detectCores())

    }
    doSNOW::registerDoSNOW(cluster)
  }

  #dimension de la quadrature (en focntion du nb de param avec des effets aleatoires)
  N=length(mu)

  i = 1

  #affichage d'une barre de progression
  ntasks <- N #nb de taches
  if(verbose){
    pb <- txtProgressBar(max = ntasks, style = 3) #on cree une barre de progression de style [======= ]
    #on def une fonction de mise a jour de la barre, fonction du numero d'iterztion n
    progress <- function(n) setTxtProgressBar(pb, n)
    #on def l'appel
    opts <- list(progress = progress)
  }else{
    opts <- NULL
  }

  #on lance en parallele, on charge le package, on lance la barre de progression
  #pour chaque individu i
  res = foreach::foreach(i = 1:N,.packages = "REMixed",.options.snow=opts)%dopar%{
    #si il y a des 0 sur la diagonale de Omega
    if(0 %in% diag(Omega[[i]])){
      #on remplace par quelque chose de tres petit
      diag(Omega[[i]])[diag(Omega[[i]])==0] <- 10**(-5)
    }
    #sinon, on recupere la LLi
    LLi = gh.LL.ind_multi(mu_i = mu[[i]],
                    Omega_i = Omega[[i]],
                    theta = theta,
                    alphas = alphas,
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
  #dans res on a tous les Li, on en fait le produit
  L = prod(sapply(res,FUN=function(ri){ri$Li}))

  #pour la LL on somme les LLi
  LL = sum(sapply(res,FUN=function(ri){ri$LLi}))

  #on ferme la barre de progression ou le parallele
  if(verbose)
    close(pb)
  if(parallel){
    snow::stopCluster(cluster)
  }

  #si on ne veut pas seuleemnt la LL
  if(!onlyLL){
    if(length(alpha$alpha1)!=1){
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

#--------------sans parallelisation-------------
gh.LL1 <- function(
    dynFUN,
    y,
    mu = NULL,
    Omega = NULL,
    theta = NULL,
    alphas = NULL,
    covariates = NULL,
    ParModel.transfo = NULL,
    ParModel.transfo.inv = NULL,
    Sobs = NULL,
    Robs = NULL,
    Serr = NULL,
    Rerr = NULL,
    ObsModel.transfo = NULL,
    data = NULL,
    n = NULL,
    prune = NULL,
    parallel = FALSE,
    ncores = NULL,
    onlyLL = FALSE,
    verbose = TRUE
) {

  if (is.null(data)) {
    test <- sapply(c("mu", "Omega", "theta", "alphas", "covariates",
                     "ParModel.transfo", "ParModel.transfo.inv", "Sobs",
                     "Robs", "Serr", "Rerr", "ObsModel.transfo"), FUN = function(x) is.null(get(x)))
    if (any(test))
      stop("Please provide all necessary arguments.")
  } else {
    if ((is.null(mu) || is.null(Omega)) & !all(c("mu", "Omega") %in% names(data))) {
      stop("Please provide mu and Omega if these are missing from data.")
    }
    for (d in names(data)) {
      if (is.null(get(d))) {
        assign(d, data[[d]])
      }
    }
  }
  # print(mu)
  # print(Omega)
  # print(theta)
  # print(alphas)
  # print(covariates)
  # print(ParModel.transfo)
  # print(ParModel.transfo.inv)
  # print(Sobs)
  # print(Robs)
  # print(Serr)
  # print(Rerr)


  if (is.null(n)) {
    if (length(theta$psi_pop) == 1) {
      n <- 100
    } else if (length(theta$psi_pop) == 2) {
      n <- 10
    } else {
      n <- 7
    }
  }

  N <- length(mu)

  res <- vector("list", N)

  if (verbose) {
    pb <- txtProgressBar(max = N, style = 3)
  }

  for (i in 1:N) {
    if (0 %in% diag(Omega[[i]])) {

      diag(Omega[[i]])[diag(Omega[[i]]) == 0] <- 1e-5
    }


   res[[i]] <- gh.LL.ind(
      mu_i = mu[[i]],
      Omega_i = Omega[[i]],
      theta = theta,
      alphas = alphas,
      dynFUN = dynFUN,
      y = y,
      covariates_i = covariates[i, , drop = FALSE],
      ParModel.transfo = ParModel.transfo,
      ParModel.transfo.inv = ParModel.transfo.inv,
      Sobs_i = lapply(Sobs, function(S) S[S$id == i, ]),
      Robs_i = lapply(Robs, function(R) R[R$id == i, ]),
      Serr = Serr,
      Rerr = Rerr,
      ObsModel.transfo = ObsModel.transfo,
      ind = i,
      n = n,
      prune = prune,
      onlyLL = FALSE
    )

    print('i=')
    print(i)
    print(res[[i]])

    if (verbose) setTxtProgressBar(pb, i)
  }

  if (verbose) close(pb)

  L <- prod(sapply(res, function(ri) ri$Li))
  LL <- sum(sapply(res, function(ri) ri$LLi))

  if (!onlyLL) {
    if (length(alphas$alpha1) != 1) {
      dLL <- rowSums(sapply(res, function(ri) ri$dLLi))
      ddLL <- matrix(rowSums(sapply(res, function(ri) ri$ddLLi)), ncol = length(dLL))
    } else {
      dLL <- sum(sapply(res, function(ri) ri$dLLi))
      ddLL <- matrix(sum(sapply(res, function(ri) ri$ddLLi)), ncol = length(dLL))
    }
    return(list(L = L, LL = LL, dLL = dLL, ddLL = ddLL))
  } else {
    return(LL)
  }
}


# ---------- ONE-DIMENSIONAL ADAPTATIVE GAUSS HERMITE  ----------
agauss.hermite <- function(n,mu=0,sd=1){
  #on recupere les poids et les racines des polynomes d'hermite
  gh <- fastGHQuad::gaussHermiteData(n)
  #on les appelle comme suit :
  names(gh) <- c("Points","Weights")

  #on calcule les bons poids et points
  ww <- gh$Weights * sqrt(2) * sd * exp(gh$Points**2)
  xx <- mu + sd * sqrt(2) * gh$Points

  #on les met dans un tableau
  return(data.frame(Points=xx,Weights=ww))
}

# ---------- MULTIDIMENSIONAL ADAPTATIVE GAUSS HERMITE  ----------
amgauss.hermite <- function(n,mu=rep(0,ncol(Omega)),Omega=diag(rep(1,length(mu))),prune=NULL){
  #dimension
  dm = length(mu)
  #on recupere les points et les poids pour un polynome de degre n en dimension 1
  gh  <- as.data.frame(fastGHQuad::gaussHermiteData(n))
  #on a :
  #         x       w
  # 1       x1      w1
  # ...
  # n       xn      wn


  idx <- as.matrix(expand.grid(rep(list(1:n),dm)))
  #rep(list(1:n),dm) : liste de dm elements avec chacun n composantes
  #ce qui colle bien avec l'idee de faire dm GH de degre n
  #expand.grid(...) : cree toutes les combinaisons de tous les vecteurs de la liste (dans un data.frame)
  #exemple
  #expand.grid(list(1:2, 1:2))
  #    Var1 Var2
  # 1    1    1
  # 2    2    1
  # 3    1    2
  # 4    2    2
  #as.matrix(...) : convertit de dataframe en matrice)
  #donc on a une matrice avec autant de lignes que de combinaisons possibles de dm listes de n elements

  pts <- matrix(gh[idx,1],nrow(idx),dm)
  #gh[idx,1] : acceder aux elements de la premiere colonne de gh (les x)
  #en utilisant les indices de idx

  #on fait une matrice
  #avec dm colonnes (taille de la dimension), chaque ligne est un point
  #avec autant de ligne qu'on peut faire de combinaison de dm listes de n elements (ie nb de lignes de idx)

  #pour la ligne L :
  #on regarde la Le combinaison des dm elements : 1 2 5 ... dm
  #on obtient la Le ligne de la matrice pts :  x1 x2 x5 ... xdm

  wts <- apply(matrix(gh[idx,2],nrow(idx),dm), 1, prod)
  #matrix(gh[idx,2],nrow(idx),dm) : matrice des poids de GH selon les combinaisons de idx
  #apply(...,1,prod) : applique prod() a chauqe ligne de la matrice
  #ie : on a un vecteur colonne ou chaque composante est le produit (de taille dm) des poids de la combinaison d'idx

  #si on a un coefficient de pruning non nul
  if(!is.null(prune)) {
    qwt <- quantile(wts, probs=prune) #quantile d'ordre prune ie valeur du poids tq prune% des poids soient inferieurs a cette valeur

    pts <- pts[wts > qwt,]
    #wts > qwt renvoie TRUE si le poids est superieur au seuil qwt, FALSE sinon
    #pts[wts > qwt] selectionne les lignes de pts pour lesquelles le poids correspondant est > qwt.
    # , : indique que toutes les colonnes sont conservees
    wts <- wts[wts > qwt]
    #on ne garde que les poids superieurs au quantile
  }

  eig <- eigen(Omega)
  #eigen() : fonction qui calcule les valeurs propres (values) et les vecteurs propres (vectors) d'une matrice carree
  S = eig$vectors #matrice des vecteurs propres
  L = diag(sqrt(eig$values)) #matrice diagonale racine des valeurs propres
  #Omega = SL^2S^T
  rot <- S %*% L #multiplication matricielle entre S et L

  ppts  <-t( sqrt(2) * rot %*% t(pts) + mu)
  wwts  <- wts * sqrt(2)**dm * abs(det(rot)) * exp(apply(pts,1,FUN=function(x){t(x)%*%x}))


  return(list(Points=ppts, Weights=wwts))
}

# ---------- INDIVIDUAL CONTRIBUTION TO LOG-LIKELIHOOD GRADIENT IN ALPHAl=0 ----------

# ---------- 1. Pour chaque individu  ----------
lambda.max.ind_multi <- function(
    dynFUN,
    y,
    mu_i=NULL,
    Omega_i=NULL,
    theta=NULL,
    alphas=NULL,
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


#on reinitialise toutes els variables
  mu <- Omega <- Sobs <- Robs <- covariates <- NULL

#si on n'a pas fourni de 'data'
#on verifie qu'on a fourni les autres arguments qu'il faut
  if(is.null(data)){
    test <- sapply(c("mu_i","Omega_i","theta","alphas","covariates_i","ParModel.transfo","ParModel.transfo.inv","Sobs_i","Robs_i","Serr","Rerr","ObsModel.transfo"),FUN=is.null)
    if(any(test))
      stop("Please provide all necessary arguments.")
  }else{
    #si on a fournit 'data'
    #on verifie qu'on a bien fait tourne EBE dans monolix et qu'on a tout ce qu'il faut
    if((is.null(mu_i) || is.null(Omega_i)) & !all(c("mu","Omega") %in% names(data))){
      stop("Please provide mu and Omega if these are missing from data.")
    }
  #arguments non indivudels
    argsNONind = c("theta","alphas","ParModel.transfo","ParModel.transfo.inv","Serr","Rerr","ObsModel.transfo")
    #on attribue la bonne valeur a la bonne variable
    for(d in 1:length(argsNONind)){
      test <- eval(parse(text=paste0("is.null(",argsNONind[d],")")))
      if(test){
        eval(parse(text=paste0(argsNONind[d],"<- data[[argsNONind[d]]]")))
      }
    }
  #arguments indiv
    argsInd = c("mu","Omega","Robs","Sobs","covariates")
    for(d in 1:length(argsInd)){
      test <- eval(parse(text=paste0("is.null(",paste0(argsInd[d],"_i"),")")))
      if(test){
        eval(parse(text=paste0(argsInd[d],"<- data[[argsInd[d]]]")))

      }
    }
  #s'il manque un argument individuel ET le numero de l'individu -> impossible
    if((is.null(mu_i) || is.null(Omega_i) || is.null(Robs_i) || is.null(Sobs_i) || is.null(covariates_i)) & is.null(ind)){
      stop("Please provide individual information (mu_i,Omega_i,Sobs_i,Robs_i,covariates_i), or the individual id ind.")
    }
  #s'il manque des arguments individuels mais qu'on a le numero de l'indiv
  #on les recupere depuis monolix
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

  #dimension de la quadrature = nb de parametres avec des effets aleatoires
  dm = length(mu_i)

  #on recupere les points et les poids de la quadrature
  #si on est en dimension >1
  if(dm!=1){
    mh.parm <- amgauss.hermite(n,mu=mu_i,Omega=Omega_i,prune=prune)
    detsq.omega = prod(theta$omega)
    root.omega = diag(1/theta$omega**2)
  #et si on est en dimension < 1
  }else{
    mh.parm <- agauss.hermite(n,mu = mu_i,sd=sqrt(as.numeric(Omega_i)))
    detsq.omega = as.numeric(theta$omega)
    root.omega = 1/as.numeric(theta$omega)**2
  }

  #nombre de points de la quadrature
  nd = length(mh.parm$Weights)

  #nombre de biomarqueurs
  R.sz = length(Rerr)

  #nombre de comapritments observes
  S.sz = length(Serr)

  #tous les temps d'observation
  all.tobs = sort(union(unlist(lapply(Robs_i,FUN=function(x){x$time})),unlist(lapply(Sobs_i,FUN=function(x){x$time}))))

  #on simule la dynamique deterministe selon les EDO a chaque instant pour chaque effet aleatoire
  dyn <- setNames(lapply(split(mh.parm$Points,1:nd),FUN=function(eta_i){
    PSI_i  = indParm(theta[c("phi_pop","psi_pop","gamma","beta")],covariates_i,setNames(eta_i,colnames(Omega_i)),ParModel.transfo,ParModel.transfo.inv)
    dyn_eta_i <- dynFUN(all.tobs,y,unlist(unname(PSI_i)))

    return(dyn_eta_i)
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
  lapply(1:L,FUN=function(comp){
    setNames(sapply(1:nd,FUN=function(ei){
      dyn.ei = dyn[[paste0("eta_",ei)]]
      yGk = ObsModel.transfo$linkR[k]
      sig = Rerr[[yGk]]
      tki = Robs_i[[yGk]]$time
      Zki = Robs_i[[yGk]][,yGk]

      a0 = theta$alpha0[[yGk]]

      atot = get_all_alpha_l_for_one_k(alphas, yGk, L)

      trs = get_all_transfo_l_for_one_k(ObsModel.transfo, yGk, L)
      Rki = lapply(names(trs), function(nom) {

        trs[[nom]](dyn.ei[dyn.ei[,"time"] %in% tki, nom])
      })

      return(1/sig**2*sum(Rki[[comp]]*(Zki-a0))*S.margDensity[[ei]])

    }),paste0("eta_",1:nd))
  })})

  indi.contrib <- sapply(1:R.sz,FUN=function(k){
    sapply(1:L,FUN=function(comp){
    sum(mh.parm$Weights*S2.margDensity[[k]][[comp]])
  })})/sum(mh.parm$Weights*S.margDensity)

  return(indi.contrib)
}

# ---------- 2. Calcul du lambda max pour la grille de penalisation  ----------
lambda.max_multi  <- function(
    dynFUN,
    y,
    mu=NULL,
    Omega=NULL,
    theta=NULL,
    alphas=NULL,
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
    test <- sapply(c("mu","Omega","theta","alphas","covariates","ParModel.transfo","ParModel.transfo.inv","Sobs","Robs","Serr","Rerr","ObsModel.transfo"),FUN=is.null)
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

  res = foreach::foreach(i = 1:N,.packages = "REMixed",.export = c("lambda.max.ind_multi","amgauss.hermite"),.options.snow=opts)%dopar%{
    if(0 %in% diag(Omega[[i]])){
      diag(Omega[[i]])[diag(Omega[[i]])==0] <- 10**(-5)
    }
    lambda.max.ind_multi(mu_i = mu[[i]],
                   Omega_i = Omega[[i]],
                   theta = theta,
                   alphas = alphas,
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

#-------------3. sans parallelisation -------------
lambda.max1  <- function(
    dynFUN,
    y,
    mu=NULL,
    Omega=NULL,
    theta=NULL,
    alphas=NULL,
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
    test <- sapply(c("mu","Omega","theta","alphas","covariates","ParModel.transfo","ParModel.transfo.inv","Sobs","Robs","Serr","Rerr","ObsModel.transfo"),FUN=is.null)
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
      cluster <- parallel::makeCluster(ncores)
    }else{
      cluster <- parallel::makeCluster(parallel::detectCores())
    }
    doSNOW::registerDoSNOW(cluster)
  }

  N=length(mu)

  res <- vector("list", N)

  if (verbose) {
    pb <- txtProgressBar(max = N, style = 3)
  }

  for (i in 1:N) {
    if (0 %in% diag(Omega[[i]])) {

      diag(Omega[[i]])[diag(Omega[[i]]) == 0] <- 1e-5
    }


    res[[i]] <- lambda.max.ind_multi(mu_i = mu[[i]],
                               Omega_i = Omega[[i]],
                               theta = theta,
                               alphas = alphas,
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


    if (verbose) setTxtProgressBar(pb, i)
  }


  if(verbose)
    close(pb)
  if(parallel)
    snow::stopCluster(cluster)
  return(max(Reduce("+",res)))
}


# ---------- FUNCTIONS TO DECOMPOSE IN EXP & MANT ----------
decomp <- function(x){
  #si on veut decomposer 0
  if (x == 0) {
    return(c(exponent = 0, mantissa = 0))
  }
  #on calcule l'exposant en base 10 et on arrondi en dessous
  exponent_x <- floor(log10(abs(x)))
  #on prend le plus petit exposant que la machine tolere
  exponent_lim <- abs(floor(log10(.Machine$double.xmin)))
  #si l'exposant depasse la limite de la machine
  if(abs(exponent_x)>exponent_lim){
    #on prned la decomposition qu'on aurait si exponent_x etait la limite de la machine
    mantissa_x <- x*10**(-sign(exponent_x)*(exponent_lim-1))
    #puis pn fait la difference entre la vraie valeur de exponent_x et celle qu'on a deja utilisee pour la premiere decomposition
    mantissa_x <- mantissa_x*10**(-(exponent_x+(exponent_lim-1)))
  }
  #sinon
  else{
    #on recupere l'ecriture scientifique de x
    mantissa_x <- x*10**(-exponent_x)
  }

  return(c(exponent=exponent_x,mantissa=mantissa_x))
}



# ---------- FUNCTIONS TO GET ALL THE _l FOR ONE FIXED k  ----------
# ---------- 1. to get all the alpha_l for one k  ----------
#on doit recuperer pour un k donne les alpha_lk
get_all_alpha_l_for_one_k <- function(alphas, k, L) {
  if(L != length(alphas)){
    stop(paste0("The number of unobserved compartments L = ",L,"  must be the same as the length of the list alphas = ",length(alphas)," (that contains the L lists of the K values alpha, i.e., for each compartment l, one value per biomarker k)."))
  }
    a <- list()
    for (l in 1:L) {
      alphal <- paste0("alpha", l)
      #on recupere la valeur k de alpha l
      a[[alphal]] <- alphas[[alphal]][[k]]
    }
    return(a)
    #ça renvoie une liste contenant tous les alpha_l pour un k donne
}

# ---------- 2. to get all the transformations for one k ----------

get_all_transfo_l_for_one_k <- function(ObsModel.transfo, k, L) {
  #on recupere le nombre de compartiments latent Rl dans ObsModel.transfo
  Rl_names <- grep("^R[0-9]+$", names(ObsModel.transfo), value = TRUE)
  #ce qu'on filtre avec grep :
    # ^	: debut de la chaine
    # R	: lettre R
    # [0-9]+	: 1 ou plusieurs chiffres (+ signifie au moins 1)
    # $	: fin de la chaîne
  #value = TRUE car de base grap renvoie les indices et nous on veut vraiment les noms

  #si on n'a pas autant de compartiments Rl que L :
  if (length(Rl_names) != L) {
    stop(paste0("The number of unobserved compartments  L = ", L, " must be the same as the length of the number of comprtments Rl = ", length(Rl_names), " in ObsModel.transfo."))
  }
  trs <- list()
  for (l in 1:L) {
    #pour chaque compartiment
    Rl <- paste0("R", l)
    #on recupere la transformation qui correspond au compartiment et au biomarqueur
    transfo <- ObsModel.transfo[[Rl]][[which(ObsModel.transfo$linkR == k)]]
    #on nomme les composantes de trs comme il faut
    #cad avec le nom du compartiment auquel la transfo correspond
    for (fun_name in names(transfo)) {
      trs[[fun_name]] <- transfo[[fun_name]]
    }
  }
  return(trs)
}
