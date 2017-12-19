#' Bootstrap
#'
#' This function generates bootstrapped betas.
#' @param formu formula to be used in the model
#' @param df data frame to be used in the model
#' @param model type of the model, "glm" and "glmm" currently available, "bayes" coming soon
#' @param n_bootstrap number of bootstraps
#' @param pop.size number of individuals
#' @param stratified logical wether a stratified or a pooled bootstrap should be performed
#' @param bayesianInput_ a list containing a Bayesian Setup, settings, and a character specifying the sampler created with Bayesian Tools
#'
#' @return returns betas of the bootstrapped models
#'
#' @keywords bootstrap
#'
#' @examples
#' bootstrap()
#'
#' @export

bootstrap <- function(formu = f.species, df = track, model = "glm", n_bootstrap = 100, pop_size = pop.size, stratified = F, bayesianInput_ = NA){
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  doSNOW::registerDoSNOW(cl)

  betas <- foreach(i = 1:n_bootstrap, .combine = rbind, .packages = c("lme4", "BayesianTools"))%dopar%{
    bootstrap = i

    if(stratified == T){
      ind.samples <- vector("list", pop_size)
      for(ind in 1:pop_size){
        ind.samples[[ind]] <- sample(which(df$ind == ind), size = sum(df$ind == ind), replace = T)
      }
      boot.sample <- unlist(ind.samples)
    }else{
      boot.sample <- sample(1:nrow(df), size = nrow(df), replace = T)
    }

    if(model == "glm") return(glm(formu, data = df[boot.sample, ], family = "binomial")$coefficients)
    
    if(model == "glmm") return(glmer(formu, data = df[boot.sample, ], family = "binomial")@beta)
    
    if(model == "bayes"){
      
        obs <- bayesianInput_[[1]][boot.sample]
      
        likelihood <- function(x, sum = TRUE){
        predicted <- plogis(bayesianInput_[[5]][boot.sample, ] %*% x)
        llValues <- dbinom(obs, size=1, predicted, log=T)
        if (sum == FALSE)return(llValues)
        else return(sum(llValues))
      }
      
      
      bayesianSetup_ <- createBayesianSetup(likelihood = likelihood, 
                                            prior = bayesianInput_[[4]], 
                                            names = c("Intercept", colnames(bayesianInput_[[5]])[-1]))
      
      return(MAP(runMCMC(bayesianSetup = bayesianSetup_, 
                         settings = bayesianInput_[[2]], 
                         sampler = bayesianInput_[[3]]))$parametersMAP)
      }
    
  }

  parallel::stopCluster(cl)
  return(betas)
}
