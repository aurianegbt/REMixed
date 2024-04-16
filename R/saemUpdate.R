#' SAEM update of parameters
#'
#' Compute iterations of SAEM using monolix to estimates population parameters with fixed \eqn{alpha_1} parameters.
#'
#' @details
#'
#' We have a differential system of equations containing variables \eqn{(S_{p})_{p\leq P}} and \eqn{R}, that depends on some parameters. We define a non-linear mixed effects model depending on individual parameters \eqn{(\psi_i)_{i\leq N}} that defines some of the paramters from the structural model as individuals using a generalized linear model for each parameter \eqn{l\leq m}, and individual \eqn{i\leq N}: \deqn{h_l(\psi_{li}) = h_l(\psi_{lpop})+X_i\beta_l + \eta_{li}} with the covariates of individual \eqn{i\leq N}, \eqn{(X_i)_{i\leq N}}, random effects \eqn{\eta_i=(\eta_{li})_{l\leq m}\overset{iid}{\sim}\mathcal N(0,\Omega)}, the population parameters \eqn{\psi_{pop}=(\psi_{lpop})_{l\leq m}} and \eqn{\beta=(\beta_l)_{l\leq m}} is the vector of covariates effects on parameters. The rest of the population parameters of the structural model, that hasn't random effetcs, are denoted by \eqn{\phi_i}, and are defined as \deqn{\phi_i=\phi_{pop} + X_i \gamma }, they can depends on covariates effects or be constant over the all population. \cr To simplify formula, we write the individual process over the time, resulting from the differential system for a set of parameters \eqn{\phi_i,(\psi_{li})_{l\leq m}} for individual \eqn{i\leq N}, as \deqn{S_{p}(\cdot,\phi_i,(\psi_{li})_{l\leq m})=S_{pi}(\cdot)}, \eqn{p\leq P} and \deqn{R(\cdot,\phi_i,(\psi_{li})_{l\leq m})=R_i(\cdot)}. We assume that individual trajectories \eqn{(S_{pi})_{p\leq P,i\leq N}} are observed through a direct observation model, up to a transformation \eqn{g_p}, \eqn{p\leq P}, at differents times \eqn{(t_{pij})_{i\leq N,p\leq P,j\leq n_{ip}}} : \deqn{ Y_{pij}=g_p(S_{pi}(t_{pij}))+\epsilon_{pij} } with error \eqn{\epsilon_p=(\epsilon_{pij})\overset{iid}{\sim}\mathcal N(0,\varsigma_p^2)} for \eqn{p\leq P}. \cr The individual trajectories \eqn{(R_{i})_{i\leq N}} are observed through the latent processes, up to a transfromation \eqn{s_k}, \eqn{k\leq K}, observed in \eqn{(t_{kij})_{k\leq K,i\leq N,j\leq n_{kij}}} : \deqn{ Z_{kij}=\alpha_{k0}+\alpha_{k1} s_k(R_i(t_{kij}))+\varepsilon_{kij}} where \eqn{\varepsilon_k\overset{iid}{\sim} \mathcal N(0,\sigma_k^2)}. The parameters of the model are then \eqn{\theta=(\phi_{pop},\psi_{pop},B,\beta,\Omega,(\sigma^2_k)_{k\leq K},(\varsigma_p^2)_{p\leq P},(\alpha_{0k})_{k\leq K})} and \eqn{\alpha=(\alpha_{1k})_{k\leq K}}.
#'
#' Here we estimates the parameters \eqn{\theta} of the monolix project loaded or given while fixing regularization parameters \eqn{\alpha_1} parameters to the values contained in `a.final` arguments.
#'
#' If the pop.set argument is set to NULL, then the settings of the population parameters algorithm is set to 50 iterations for smoothing and 50 for exploratory phases, no simulated annealing, no auto stop for smoothing and exploratory phases.
#'
#' @param project string with the initial monolix project (default to NULL, and if so, used the currently loaded project) ;
#' @param final.project string with the final Monolix project (default adds "_built" to the original project) ;
#' @param alpha list of named vector "alpha0", "alpha1" (in good order), alpha1 mandatory even if 1 ;
#' @param a.final values of alpha1 parameters to fixed, in the same order as in alpha$alpha1 ;
#' @param iter current iteration (default NULL)
#' @param pop.set monolix population setting for SAEM, (default NULL).
#'
#' @return a list with the SAEM iterations value and the final estimates.
#' @export
#'
#' @examples
#' #TODO
saemUpdate <- function(project = NULL,final.project=NULL,
                       alpha, a.final,
                       iter=NULL,
                       pop.set=NULL){

  if(is.null(project)){
    project <- Rsmlx:::mlx.getProjectSettings()$project
  }else{
    Rsmlx:::mlx.loadProject(project)
  }

  if (is.null(final.project)){
    if(!is.null(iter))
      final.project <- paste0(sub(pattern = "(.*)\\..*$",
                                  replacement = "\\1", project), "_iter",iter,".mlxtran")
    else
      final.project <- paste0(sub(pattern = "(.*)\\..*$",
                                replacement = "\\1", project), "_upd.mlxtran")}
  if (!grepl("\\.", final.project))
    final.project <- paste0(final.project, ".mlxtran")
  if (!grepl("\\.mlxtran", final.project))
    stop(paste0(final.project, " is not a valid name for a Monolix project (use the .mlxtran extension)"),
         call. = FALSE)
  project.temp <- gsub(".mlxtran", "_temp.mlxtran", project)
  project.built <- final.project

  Rsmlx:::mlx.saveProject(project.temp)

  pset <- list(nbsmoothingiterations=50,nbexploratoryiterations=50,
               simulatedannealing=F, smoothingautostop=F, exploratoryautostop=F)
  pop.set <- Rsmlx:::mlx.getPopulationParameterEstimationSettings()
  pop.set <- modifyList(pop.set, pset[intersect(names(pset), names(pop.set))])

  Rsmlx:::mlx.setInitialEstimatesToLastEstimates(fixedEffectsOnly = FALSE)

  for(k in 1:length(alpha$alpha1)){
    eval(parse(text=paste0("lixoftConnectors::setPopulationParameterInformation(",alpha$alpha1[k],"_pop=list(initialValue=",a.final[k],",method='FIXED'))")))
  }

  Rsmlx:::mlx.setPopulationParameterEstimationSettings(pop.set)

  Rsmlx:::mlx.saveProject(project.temp)

  Rsmlx:::mlx.runPopulationParameterEstimation()

  Rsmlx:::mlx.saveProject(final.project)

  re = list(SAEMiterations = Rsmlx:::mlx.getSAEMiterations(),
            estimates = Rsmlx:::mlx.getEstimatedPopulationParameters())
  if(!is.null(iter))
    re <- append(re,list(iter=iter))

  return(re)
}
