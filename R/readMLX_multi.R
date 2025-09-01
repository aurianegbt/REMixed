#' Extract Data for REMix Algorithm from a Monolix Project
#'
#' @description
#' This function retrieves all necessary information from a Monolix project file to format the input for the REMix package. It gathers all relevant data required for the REMix algorithm.
#'
#' @details
#' To simplify its use, functions \code{\link{remix}}, \code{\link{cv.remix}}, \code{\link{gh.LL}} can be used with arguments \code{data} rather than all necessary informations "\code{theta}", "\code{alpha1}", "\code{covariates}", "\code{ParModel.transfo}", "\code{ParModel.transfo.inv}", "\code{Sobs}", "\code{Robs}", "\code{Serr}", "\code{Rerr}", "\code{ObsModel.transfo}" that could be extract from a monolix project. If the SAEM task of the project hasn't been launched, it's the initial condition and not the estimated parameters that are returned. If the conditional distribution estimation task has been launched, parameters "\code{mu}" and "\code{Omega}" are returned too.
#'
#' @param project directory of the Monolix project (in .mlxtran). If NULL, the current loaded project is used (default is NULL).
#' @param ObsModel.transfo list containing two lists of transformations and two vectors linking each transformations to their observation model name in the Monolix project. The list should include identity transformations and be named \code{S} and \code{R}. The two vectors should be named \code{linkS} and \code{linkR}.
#'
#' Both \code{S} (for the direct observation models) and \code{linkS}, as well as \code{R} (for latent process models) and \code{linkR}, must have the same length.
#'
#' \itemize{
#'   \item\code{S}: a list of transformations for the direct observation models. Each transformation corresponds to a variable \eqn{Y_p=h_p(S_p)}, where the name indicates which dynamic is observed (from \code{dynFUN});  \item\code{linkS} : a vector specifying the observation model names (that is used in the monolix project, \code{alpha1}, etc.) for each transformation, in the same order as in \code{S};
#'
#'   \item\code{R}: similarly, a list of transformations for the latent process models. Although currently there is only one latent dynamic, each \eqn{s_k, k\leq K} transformation corresponds to the same dynamic but may vary for each \eqn{Y_k} observed. The names should match the output from \code{dynFUN}; \item \code{linkR} : a vector specifying the observation model names for each transformation, in the same order as in \code{R}.
#' }
#' @param alpha named list of named vector "\code{alpha0}", "\code{alpha1}" (all \code{alpha1} are mandatory). The names of \code{alpha$alpha0} and \code{alpha$alpha1} are the observation model names from the monolix project to which they are linked (if the observations models are defined whithout intercept, alpha$alpha0 need to be set to the vector NULL).
#'
#' @return A list containing parameters, transformations, and observations from the Monolix project in the format needed for the REMix algorithm :
#' \itemize{
#'   \item\code{mu} list of individuals random effects estimation (vector of r.e. need to be named by the parameter names), use to locate the density mass (if conditional distribution estimation through Monolix has been launched);
#'   \item\code{Omega} list of individuals estimated standard deviation diagonal matrix (matrix need to have rows and columns named by the parameter names), use to locate the density mass (if conditional distribution estimation through Monolix has been launched);
#'   \item\code{theta} list of model parameters containing i\itemize{\item\code{phi_pop} : named vector with the population parameters with no r.e. \eqn{(\phi_{l\ pop})_{l\leq L}} (NULL if none) ;
#' \item\code{psi_pop} : named vector with the population parameters with r.e. \eqn{(\psi_{l\ pop})_{l\leq m}} ;
#' \item\code{gamma} : named list (for each parameters) of named vector (for each covariates) of covariate effects from parameters with no r.e. ;
#' \item\code{beta} : named list (for each parameters) of named vector (for each covariates) of covariate effects from parameters with r.e..
#' \item\code{alpha0} : named vector of \eqn{(\alpha_{0k})_{k\leq K}} parameters (names are identifier of the observation model, such as in a Monolix project);
#' \item\code{omega} : named vector of estimated r.e. standard deviation;}
#'   \item\code{alpha1}  named vector of regulatization parameters \eqn{(\alpha_{1k})_{k\leq K}}, with identifier of observation model as names;
#'   \item\code{covariates} matrix of individual covariates (size N x n). Individuals must be sorted in the same order than in \code{mu} and \code{Omega};
#'   \item\code{ParModel.transfo} named list of transformation functions \eqn{(h_l)_{l\leq m}} and \eqn{(s_k)_{k\leq K}} for the individual parameter model (names must be consistent with \code{phi_pop} and \code{psi_pop}, missing entries are set by default to the identity function ;
#'   \item \code{ParModel.transfo.inv} named list of inverse transformation functions for the individual parameter model (names must be consistent with \code{phi_pop} and \code{psi_pop} ;
#'   \item\code{Sobs} ist of individuals trajectories for the direct observation models \eqn{(Y_{pi})_{p \leq P,i\leq N}}. Each element \eqn{i\leq N} of the list, is a list of \eqn{p\leq P} data.frame with time  \eqn{(t_{pij})_{j\leq n_{ip}}} and observations  \eqn{(Y_{pij})_{j\leq n_{ip}}}. Each data.frame is named with the observation model identifiers ;
#'   \item\code{Robs} list of individuals trajectories for the latent observation models \eqn{(Z_{ki})_{k \leq K,i\leq N}}. Each element \eqn{i\leq N} of the list, is a list of \eqn{k\leq K} data.frame with time  \eqn{(t_{kij})_{j\leq n_{ik}}} and observations  \eqn{(Z_{kij})_{j\leq n_{ik}}}. Each data.frame is named with the observation model identifiers ;
#'   \item\code{Serr}  named vector of the estimated error mocel constants \eqn{(\varsigma_p)_{p\leq P}} with observation model identifiers as names ;
#'   \item\code{Rerr} named vector of the estimated error mocel constants \eqn{(\sigma_k)_{k\leq K}} with observation model identifiers as names ;
#'   \item\code{ObsModel.transfo} same as input\code{ObsModel.transfo} list.
#' }
#' @export
#'
#' @seealso \code{\link{remix}}, \code{\link{cv.remix}}, \code{\link{getContinuousObservationModel}}, \code{\link{getObservationInformation}}
#'
#' @examples
#' project <- getMLXdir()
#'
#' ObsModel.transfo = list(S=list(AB=log10),
#'                         linkS="yAB",
#'                         R=rep(list(S=function(x){x}),5),
#'                         linkR = paste0("yG",1:5))
#'
#' alpha=list(alpha0=NULL,
#'            alpha1=setNames(paste0("alpha_1",1:5),paste0("yG",1:5)))
#'
#' res <- readMLX(project,ObsModel.transfo,alpha)
#'

# ---------- readMLX ----------
readMLX_multi <- function(project=NULL,
                    ObsModel.transfo,
                    alpha,L
                    ){

###VERIFICATIONS
  #on verifie bien la structure de chaque vecteur
  check.readMLX(ObsModel.transfo,alpha,L)
  #on charge le projet monolix
  #supressMessages : empeche l'affichage de messages du type "bonne importation du projet"
  suppressMessages({
    if (!is.null(project)){
      #si l'argument 'project' est fourni par l'utilisateur
      project <- prcheck(project)$project
      #sinon
    }else{
      project <- lixoftConnectors::getProjectSettings()$project #on recupere le chemin du projet mxtran
    }
  #on charge le projet
    lixoftConnectors::loadProject(project)
  })

  #Reduce(union) permet de faire d'un coup la meme chose que :
  #union(union(union(ObsModel.transfo$linkS, ObsModel.transfo$linkR), names(alpha$alpha0)), names(alpha$alpha1))

  #on verifie que tous les alphal sont bien definis et que leurs noms sont identiques a ceux dans monolix
  for(l in 1:L){
    a<- paste0("alpha",l)
    if(!all(Reduce(union,list(ObsModel.transfo$linkS,ObsModel.transfo$linkR,names(alpha$alpha0),names(alpha[[a]]))) %in% names(lixoftConnectors::getObservationInformation()))){
    stop("Inconsistent observation model between the one provided in the monolix project and the one in function inputs.")
    }}

###PARAMETRES POP
#getLaunchedTasks :
  #tasks = getLaunchedTasks()
  #tasks
  #-> $populationParameterEstimation = TRUE
  #$conditionalModeEstimation = TRUE
  #$standardErrorEstimation = "linearization"

#si l'estimation des parametres a ete faite
  if(lixoftConnectors::getLaunchedTasks()$populationParameterEstimation){
    #on stocke les estimations dans est
    est = lixoftConnectors::getEstimatedPopulationParameters()
    #exemple : getEstimatedPopulationParameters() -> [V_pop = 0.5, Cl_pop = 0.25, ka_pop = 0.05]
    #on fait un tableau avec les valeurs des parametres
    value.params = data.frame(name=names(est),initialValue=unname(est))
#si l'estimation n'a pas ete faite
  }else{
    #on les recupere mais on considere seulement les 2 premieres colonnes : name et initialValue
    #on prend toutes les lignes et toutes les colonnes sauf la 3e (la colonne method)
    value.params = lixoftConnectors::getPopulationParameterInformation()[,-3]
  }

###COVARIABLES
  #getCovariateInformation() renvoie une liste avec name, type (continous...), modalityNumber, covariate (un tableau avec les valeurs)
  #on prend toutes les colonnes sauf la premiere : id
  covariates = lixoftConnectors::getCovariateInformation()$covariate[,-1]
  #on regarde si on a ou pas des covariables
  cov.exists = !is.null(covariates) #renvoie TRUE si ce n'est pas null

  #PARAMETRES INDIVIDUELS
  IndividualParameterModel <- lixoftConnectors::getIndividualParameterModel()
  #detail de la fonction 'getIndividualParameterModel()' :
    #formula: (character) the formula used for each parameter
    #covariateModel: (logical) a list giving, for each individual parameter, a vector with TRUE for each covariate that is included in the model for that parameter or FALSE if not.
    #If there are no covariates in the model, this is an empty list.
  #si on a bien des covariables et donc si cov.exists=TRUE
  if(cov.exists){
    #on recupere la formule et covariateModel
    formula = IndividualParameterModel$formula
    CovariateModel <- IndividualParameterModel$covariateModel
    #si les parametres n'ont pas les bons noms (test avec la fonction ok.beta)
    if(!ok.beta(formula,CovariateModel)){
      stop("Please take care of having standard covariates coefficients names in your monolix project, of form 'beta_param_cov'.")
    }
    #si les omega ne sont pas present pour chaque covariable (test avec ok.omega)
    if(!ok.omega(IndividualParameterModel$variability$id,value.params)){
      stop("Please take care of having standard covariates coefficients names in your monolix project, of form 'omega_param'.")
    }
  }

  #on recupere le nom des parametres de pop (exemple :  delta_S   phi_S delta_AB alpha_11 alpha_12 alpha_13 alpha_14 alpha_15)
  var.param <- IndividualParameterModel$variability$id #la on a les noms et les valeurs
  param <- names(var.param) #on garde les noms

  #!var.param : inverse le vecteur booleen de si oui ou non les variables ont de la variabilite
  #on prend les noms parametres qui n'ont pas de variabilite (!var.param=TRUE) et ceux qui ne sont pas dans alpha
  phi.names = param[!var.param & !(param %in% unlist(alpha))]

  #on prend ceux qui ont de la variabilité et qui ne sont pas dans alpha
  psi.names = param[var.param & !(param %in% unlist(alpha))]

  phi_pop = sapply(phi.names,FUN=function(p){ #pour chaque nom de parametre p dans phi.names
    value.params[value.params$name==paste0(p,"_pop"),"initialValue"]}) #on nomme la valeur 'nom_pop' et on prend la colonne "initialValue"
  #pareil pour les param avec variabilite
  psi_pop =  sapply(psi.names,FUN=function(p){
    value.params[value.params$name==paste0(p,"_pop"),"initialValue"]
  })

  #si cov.exists=TRUE,
  if(cov.exists){
    #la liste gamma est nommee
    gamma = setNames(lapply(phi.names,FUN=function(p){ #on construit une liste gamma et pour chaque parametre p de phi.names
      g <-  sapply(names(CovariateModel[[p]]),FUN=function(c){
        #si un parametre a une covariable (ex : l'age), alors CovariateModel renvoie pour cette variable list(age=TRUE) contre list() pour les autres
        if(CovariateModel[[p]][[c]]){ #donc si un des parametres a des covariables
          value.params[value.params$name==paste0("beta_",p,"_",c), #on cree le nom du coef et on lui attribue la valeur init
                       "initialValue"]
        }
        #sinon aucun parametre n'a de covariable
        else{
          0
        }},USE.NAMES = F) #USE.NAMES=FALSE => les elements de g ne sont pas nommes
      #si tous les beta valent 0 on met g a NULL (au lieu d'a voir un vecteur de zeros)
      if(identical(g,rep(0,ncol(covariates)))){
        g <- NULL
      }
      return(g)
    }),phi.names)#on nomme la liste gamma comme phi.names

    #exactement pareil pour psi :
    beta = setNames(lapply(psi.names,FUN=function(p){
      b <-  sapply(names(CovariateModel[[p]]),FUN=function(c){
        if(CovariateModel[[p]][[c]]){
          value.params[value.params$name==paste0("beta_",p,"_",c),
                       "initialValue"]
        }else{
          0
        }},USE.NAMES = F)
      if(identical(b,rep(0,ncol(covariates)))){
        b <- NULL
      }
      return(b)
    }),psi.names)
  }else{
    gamma = NULL
    beta = NULL
  }

  #vecteur de standard error pour les parametres qui ont de la variabilite
  omega = if(length(psi.names)!=0){ #si il y a des param qui ont de la variabilite
    setNames(sapply(psi.names,FUN=function(p){ #pour chaque param p de la liste
      value.params[value.params$name==paste0("omega_",p),"initialValue"] #on recupere les valeur initiales des omega
    }),psi.names) #chaque sd est nomme par le nom de son param
  }else{NULL} #si il n'y a pas de param avec variabilite, le vecteur des omega est nul

  #on nomme les termes de alpha0 avec les noms des biomarqueurs
  alpha0 <- setNames(rep(0,length(ObsModel.transfo$linkR)),ObsModel.transfo$linkR)
  #si le vecteur n'est pas NULL
  if(!is.null(alpha$alpha0)){
    #on prend les alpha0 non nuls
    #on prend ceux qui sont dans le projet monolix et aussi dans le vecteur alpha0
    #on prend leur valeur init
    #on les nomme 'yG...'
    alpha0[!is.null(alpha$alpha0)] <- setNames(value.params[value.params$name %in% paste0(alpha$alpha0,"_pop"),"initialValue"],names(alpha$alpha0))
  }

  alphas <- list()
  for (l in 1:L) {
    alphal <- paste0("alpha", l)
    alphas[[alphal]] <- setNames(
      sapply(alpha[[alphal]], function(p) {
        value.params[value.params$name == paste0(p, "_pop"), "initialValue"]
      }, USE.NAMES = FALSE),
      names(alpha[[alphal]])
    )
  }

  #a ce stade, le vecteur alpha compose de alpha0, alpha1 et alpha2 est toujours comme suit :
  #$alpha0
  #yyG1       yyG2       yyG3
  #"alpha_01" "alpha_02" "alpha_03"

  #$alpha1
  #yyG1       yyG2       yyG3
  #"alpha_11" "alpha_12" "alpha_13"

  #$alpha2
  #yyG1       yyG2       yyG3
  #"alpha_21" "alpha_22" "alpha_23"

  #par contre, alphas est une liste avec les alphal
  #ou alphal est un vecteur nomme avec 'yGl' en nom et la valeur initiale de Monolix en valeur

  #on cree notre vecteur theta
  theta = setNames(list(phi_pop,psi_pop,gamma,beta,alpha0,omega),c("phi_pop","psi_pop","gamma","beta","alpha0","omega"))

###Serr, Rerr
  #on recupere les noms des sd de monolyx (ex : ayAB, ayG1...)
  ErrorModel <- lixoftConnectors::getContinuousObservationModel()$parameters

  # ObsModel.transfo need to be name after ??? model name in monolix ??? or dyn ??? # need to look at that closer

  #OBSERVATIONS DIRECTES (dans l'exemple, seulement AB)
  #on attribue a chaque valeur de sd le nom qui correspond au biomarqueur
  Serr <- sapply(ObsModel.transfo$linkS,FUN=function(x){
    #on recupere la valeur de chaque sd (ex : value.param$ayAB)
    #on lui attrive la valeur initiale de Monolix
    value.params[value.params$name==ErrorModel[[x]],"initialValue"]
  },USE.NAMES = T)

  #OBSERVATIONS BIOMARQUEURS (dans l'exemple G1 à G5)
  #pareil
  Rerr <- sapply(ObsModel.transfo$linkR,FUN=function(x){
    value.params[value.params$name==ErrorModel[[x]],"initialValue"]
  },USE.NAMES = T)

###ParModel.transfo and inv
  #on recupere les distribution (normal, lognormal...) des parametres avec et sans variabilite
  ParModel <- IndividualParameterModel$distribution[union(phi.names,psi.names)]

  #on initialise deux listes
  ParModel.transfo = list()
  ParModel.transfo.inv =list()
  #pour chqaue parametre considere
  for(i in 1:length(ParModel)){
    #si la distribution du parametre est 'normal'
    if(ParModel[[i]]=="normal"){
      ParModel.transfo[[i]] <-  function(x){x}
      ParModel.transfo.inv[[i]] <- function(x){x}
    }#on garde l'identite comme transformation

    #si la distribution du parametre est 'lognormal'
    else if(ParModel[[i]]=="logNormal" || ParModel[[i]]=="lognormal"){
      ParModel.transfo[[i]] <- log #la transformation par rapport à la normale est le log
      ParModel.transfo.inv[[i]] <- exp #son inverse : exp
    }
    #si la distribution du parametre est 'logitnormal'
    else if(ParModel[[i]]=="logitNormal"){
      ParModel.transfo[[i]] <- function(x){log(x/(1-x))} #on garde la transformation logit
      ParModel.transfo.inv[[i]] <- function(x){exp(x)/(1+exp(x))} #on prend son inverse
    }
  }
  #chaque composante des deux liste est nommee comme le parametre considere
  names(ParModel.transfo) <-  names(ParModel.transfo.inv) <- names(ParModel)

###Sobs Robs
  Observation <- lixoftConnectors::getObservationInformation()

  #on nomme la colonne de valeurs osbervees (du fichier data) comme dans ObsModel.transfo
  #pour les donnees du compartiment observe
  Sobs <- setNames(lapply(ObsModel.transfo$linkS,FUN=function(x){Observation[[x]]}),ObsModel.transfo$linkS)
  #et des biomarqueurs
  Robs <- setNames(lapply(ObsModel.transfo$linkR,FUN=function(x){Observation[[x]]}),ObsModel.transfo$linkR)


#si on n'a pas lance l'estimation de la distrib conditionnelle
  #ie, si on l'a lance : getLaunchedTasks()$conditionalDistributionSampling = TRUE
  #donc !getLaunchedTasks()$conditionalDistributionSampling = FALSE
  #donc on ne rentre pas dans la boucle
  if(!lixoftConnectors::getLaunchedTasks()$conditionalDistributionSampling){
    #si on n'a pas estime mu et Omega, on renvoie tout sauf mu et Omega
    return(list(
      theta = theta,
      alphas = alphas,
      covariates = covariates,
      ParModel.transfo = ParModel.transfo,
      ParModel.transfo.inv = ParModel.transfo.inv,
      Sobs = Sobs,
      Robs = Robs,
      Serr= Serr,
      Rerr=Rerr,
      ObsModel.transfo = ObsModel.transfo
      ))
  }
  #si on a lance l'estimation des distrib conditionnelles
  else{
    random.effect = lixoftConnectors::getEstimatedRandomEffects(method="conditionalMean")
    #on recupere deux tableaux : conditionalMean et conditionalSD
    #pour chaque individu, on a le conditionalMean et SD de chaque parametre de pop
    N = length(random.effect$conditionalMean$id)

    #rappel : psi.names = vecteur de noms des params ayant de la variabilite (ex : "delta_S"  "phi_S"    "delta_AB")
    #pour chaque individu i
    Omega = lapply(1:N,FUN=function(i){
      #pour chaque parametre ayant de la variabilite (ie etant dans psi.names)
      eta = sapply(psi.names,FUN=function(p){

        random.effect$conditionalSD[i,paste0("eta_",p)]
      })
      #eta :
      #pour chaque individu, on a un vecteur de taille 3, nomme avec la valeur de la sd pour chaque parametre de psi.names
      #exemple pour les deux premiers indiv :
      #delta_S      phi_S   delta_AB
      #0.11438682 0.18654756 0.04162139
      #delta_S      phi_S   delta_AB
      #0.11794636 0.21585289 0.03949503


      #pour l'individu i,
      #si il n'y a qu'un seul parametre avec de la variabilite
      if(length(eta)==1){
        Omega_i = matrix(eta**2,nrow=1,ncol=1)
      }
      #Omega est un reel, c'est la variance du parametre (ecart type au carre)

      #s'il y en aplus d'un, on doit construire la matrice de variance-covariance
      else{
        #matrice diagonale avec la variance sur la diagonale
        #eta_i : vecteur a autant de commposantes que de parametres dans psi.names
        #eta_i**2 : pareil
        #Omega_i : matrice dont la diagonale est le vecteur eta_i**2
        Omega_i = diag(eta**2)
      }
      #on appelle les lignes et les colonnes avec les noms de chaque parametre
      rownames(Omega_i) <- colnames(Omega_i) <- psi.names
      #exemple de ce qu'on peut obtenir pour un individu :
      #.        delta_S      phi_S    delta_AB
      #delta_S  0.01550855 0.00000000 0.000000000
      #phi_S    0.00000000 0.03172105 0.000000000
      #delta_AB 0.00000000 0.00000000 0.002191794
      return(Omega_i)
    })
#Omega : resence les Omega_i

    #meme demarche pour mu mais cette fois c'est juste la moyenne
    mu = lapply(1:N,FUN=function(i){
      mu_i = sapply(psi.names,FUN=function(p){
        random.effect$conditionalMean[i,paste0("eta_",p)]
      })
      return(mu_i)
    })

    #on renvoie aussi ces estimations en sorties de readMLX
    return(list(
      mu=mu,
      Omega=Omega,
      theta = theta,
      alphas = alphas,
      covariates = covariates,
      ParModel.transfo = ParModel.transfo,
      ParModel.transfo.inv = ParModel.transfo.inv,
      Sobs = Sobs,
      Robs = Robs,
      Serr= Serr,
      Rerr=Rerr,
      ObsModel.transfo = ObsModel.transfo
    ))
  }
}

# ---------- OK FUNCTIONS ----------

#EXEMPLE
#s_kl : transformation du l-ieme compartiment du k-ieme biomarqueur

#Gij1 = alpha_01 + alpha_11 x s_11(R1) + alpha_21 x s_12(R2)
#Gij2 = alpha_02 + alpha_12 x s_21(R1) + alpha_22 x s_22(R2)`
#Gij3 = alpha_03 + alpha_13 x s_31(R1) + alpha_23 x s_32(R2)`

#ici, on voudrait ObsModel.transfo : R1=(s_11, s_12, s_13); R2=(s21, s_22, s_32);,linkR=(yG1,yG2,yG3)
#alpha = ( alpha0=(alpha_01, alpha_02, alpha_03),
#          alpha1=(alpha_11, alpha_12, alpha_13),
#          alpha2=(alpha_21, alpha_22, alpha_23) )

# ---------- 1. check.readMLX ----------
#fonction qui verifie que les ObsModel.transfo et alpha sont bien definis par l'utilisateur
#L : nb de compartiments latents
check.readMLX <- function(ObsModel.transfo,alpha,L){
  #meme taille pour le nombre de transformations du compartiment observe et le nombre de compartiments observes
  if(length(ObsModel.transfo$S)!=length(ObsModel.transfo$linkS)){
    stop(paste0("Error in length of ObsModel.transfo arguments. S (length=",length(ObsModel.transfo$S),") and linkS (length=",length(ObsModel.transfo$linkS),") should have the same length."))
  }

  #on verifie qu'il y a le bon nombre d'alpha_0 (autant que de biomarqueurs observes)
  else if(!is.null(alpha$alpha0) && length(alpha$alpha0)!=length(ObsModel.transfo$linkR)){
    stop(paste0("Error in length of alpha0 parameters. alpha$alpha0 (length=",length(alpha$alpha0),") and ObsModel.transfo$R (length=",length(ObsModel.transfo$R),") should have the same size. "))
  }

  #on verifie que les alpha0 et les compartiments latents ont bien le meme nom
  else if(!is.null(alpha$alpha0) && !all(names(alpha$alpha0)==ObsModel.transfo$linkR)){
    stop("Error in names of ObsModel.transfo$linkR and alpha0 parameters. They must be consistent (and in the same order).")
  }

  for(l in 1:L){
    #pour chaque compartiment l
    comp<-paste0("R",l) #comp devient R1, R2, ..., RL
    a<-paste0("alpha",l) #a deveint alpha1, ..., alphaL
    #on verifife que chaque Rl a bien la meme taille que le nb de biomarqueurs
    if(length(ObsModel.transfo[[comp]])!=length(ObsModel.transfo$linkR)){
      stop(paste0("Error in length of ObsModel.transfo arguments. R",l," (length=",length(ObsModel.transfo[[comp]]),") and linkR (length=",length(ObsModel.transfo$linkR),") should have the same length."))
      }
    #on verifie que chaque vecteur alphal a bien la meme taille que les Rl
    else if(length(alpha[[a]])!=length(ObsModel.transfo[[comp]])){
      stop(paste0("Error in length of alpha",l," parameters. alpha$alpha",l," (length=",length(alpha[[a]]),") and ObsModel.transfo$R",l," (length=",length(ObsModel.transfo[[comp]]),") should have the same size. "))
      }
    #on verifie que chaque alphal a bien le meme nom que les biomarqueurs observes
    else if(!all(names(alpha[[a]])==ObsModel.transfo$linkR)){
      stop(paste0("Error in names of ObsModel.transfo$linkR and alpha",l," parameters. They must be consistent (and in the same order)."))
      }
    }

  return(invisible(TRUE))
}

# ---------- 2. ok.beta ----------
#formula : getIndividualParameterModel()$formula,
#dans le getIndividualParameterModel()$formula, les formules sont separees d'un \n (cf fichier exemples)

#CovariateModel : getIndividualParameterModel()$covariateModel,

#pour verifier que les noms des param de la formule collent bien avec les covariables
ok.beta <- function(formula,CovariateModel){
  form <- strsplit(formula, "\n")[[1]] #avec cette fonction on decoupe la chaine formula en plusieurs morceaux, en utilisant \n comme separateur
  #pour chaque formule dans 'formula' ie pour chaque param indiv
  for(p in 1:length(form)){
    f <- form[p] #on stocke l'expression de la covariable p dans f
    cov <- CovariateModel[[p]] #(souvent named list())
    param <- names(CovariateModel)[p] #nom de la covariable p dans param

    if(!identical(cov,list())){
      for(c in names(cov)[cov]){
        before_c <- sub(paste0("\\*",c,".*"), "", f) #on enleve la covariable et tout ce qui suit
        beta <- sub(".*\\s", "", before_c) #on prend juste le nom du param considere
        if(beta != paste0("beta_",param,"_",c)){ #on regarde si le param s'appelle correctement
          return(FALSE)
        }
      }
    }
  }
  return(TRUE)
}

# ---------- 3. ok.omega ----------
#verififer pour chaque parametre avec une variabilite qu'il y a bien un omega_param
ok.omega <- function(variability,value.params){
  for(p in names(variability)){
    #si le parametre p a une variabilite ie variability[[p]]=TRUE
    #et si il n'y a pas d'omega_p correspondant
    if(variability[[p]] && !(paste0("omega_",p) %in% value.params$name)){
      #on a un porbleme
      return(FALSE)
    }
  }
  return(TRUE)
}

# ---------- 4. prcheck ----------

prcheck <- function (project, f = NULL, settings = NULL, model = NULL,
                     paramToUse = NULL, parameters = NULL, level = NULL, tests = NULL,
                     nboot = NULL, method = NULL)
{
  #on met tout le monde en NULL
  RsmlxDemo1.project <- RsmlxDemo2.project <- warfarin.data <- resMonolix <- NULL
  #si le nom du projet commence par RsmlxDemo
  if (identical(substr(project, 1, 9), "RsmlxDemo")) {
    #on vide tout au cas ou
    RsmlxDemo1.project <- RsmlxDemo2.project <- warfarin.data <- resMonolix <- NULL
    #on les supprime
    rm(RsmlxDemo1.project, RsmlxDemo2.project, warfarin.data,resMonolix)
    #parse(text = "data(RsmlxDemo)") --> transforme la chaîne "data(RsmlxDemo)" en une expression R (de type expression())
    #eval(...)-->evalue cette expression ie l'excecute
    eval(parse(text = "data(RsmlxDemo)"))

    tmp.dir <- tempdir() #renvoie le chemin du repertoire temporaire de l’environnement de travail
    #on creer les fichiers mlxtran "proprement"
    write(RsmlxDemo1.project, file = file.path(tmp.dir,
                                               "RsmlxDemo1.mlxtran"))
    write(RsmlxDemo2.project, file = file.path(tmp.dir,
                                               "RsmlxDemo2.mlxtran"))
    write.csv(warfarin.data, file = file.path(tmp.dir, "warfarin_data.csv"),
              quote = FALSE, row.names = FALSE)
    #on garde leur chemin
    project <- file.path(tmp.dir, project)

    demo <- TRUE
    #si il y a un parametre f non NULL donne par l'utilisateur
    if (!is.null(f)) {
      if (f == "boot") { #bootstrap
        if (is.null(settings))
          res <- resMonolix$r1.boot
        else if (!is.null(settings$N) & is.null(settings$covStrat))
          res <- resMonolix$r2.boot
        else res <- resMonolix$r3.boot
      }
      else if (f == "build") {
        if (identical(model, "all") & identical(paramToUse,
                                                "all"))
          res <- resMonolix$r1.build
        else if (identical(model, "all"))
          res <- resMonolix$r2.build
        else res <- resMonolix$r3.build
      }
      else if (f == "conf") {
        if (method == "fim" & level == 0.9)
          res <- resMonolix$r1.conf
        else if (method == "fim" & level == 0.95)
          res <- resMonolix$r2.conf
        else if (method == "proflike")
          res <- resMonolix$r3.conf
        else res <- resMonolix$r4.conf
      }
      else if (f == "cov") {
        if (identical(method, "COSSAC") & identical(paramToUse,
                                                    "all"))
          res <- resMonolix$r1.cov
        else if (identical(method, "SCM"))
          res <- resMonolix$r2.cov
        else res <- resMonolix$r3.cov
      }
      else if (f == "test") {
        if (length(tests) == 4)
          res <- resMonolix$r1.test
        else res <- resMonolix$r2.test
      }
      else if (f == "set")
        res = "foo"
    }
  }
  #si on est pas dans un projet de demo
  else {
    #si la version de Monolix est celle de 2019 ou 2020
    if (grepl("2020", Rsmlx::initRsmlx()$version) | grepl("2019",
                                                          Rsmlx::initRsmlx()$version))
      stop("Rsmlx versions above 4.0 are compatible only with MonolixSuite >= 2021R1",
           call. = FALSE)
    #sortie de initRsmlx : status: boolean equaling TRUE if the initialization has been successful
    #si Rsmlx n'est pas initialis comme il faut
    if (!Rsmlx::initRsmlx()$status)
      return()
    #si le projet n'a pas '.mlxtran'
    if (!grepl("\\.", project))
      project <- paste0(project, ".mlxtran") #on lui rajoute
    #si le projet n'existe pas
    if (!file.exists(project))
      stop(paste0("Project '", project, "' does not exist"),
           call. = FALSE)
    #on charge le projet
    lp <- lixoftConnectors::loadProject(project)
    if (!lp) #si ça ne fonctionne pas
      stop(paste0("Could not load project '", project,
                  "'"), call. = FALSE) #on affiche qu'on n'a pas reussi

    demo <- FALSE
    res <- NULL
  }
  return(list(project = project, demo = demo, res = res))
}
