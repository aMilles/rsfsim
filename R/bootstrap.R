#' Bootstrap
#'
#' This function generates bootstrapped betas.
#' @param formu formula to be used in the model
#' @param df data frame to be used in the model
#' @param model type of the model, "glm" and "glmm" currently available, "bayes" coming soon
#' @param n_bootstrap number of bootstraps
#' @param pop.size number of individuals
#' @param stratified logical wether a stratified or a pooled bootstrap should be performed
#'
#' @return returns betas of the bootstrapped models
#'
#' @keywords bootstrap
#'
#' @examples
#' bootstrap()
#'
#' @export

bootstrap <- function(formu = f.species, df = track, model = "glm", n_bootstrap = 100, pop_size = pop.size, stratified = F){
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  doSNOW::registerDoSNOW(cl)

  betas <- foreach(i = 1:n_bootstrap, .combine = rbind)%dopar%{
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
    if(model == "bayes") return("bayes not implemented")
  }

  parallel::stopCluster(cl)
  beepr::beep(1)
  return(betas)
}
