#' Bootstrap
#'
#' This function generates bootstrapped models with glm/glmm .
#' @param f_species formula to be used in the analysis
#' @param df data frame to be used in the analysis
#' @param analysis type of the analysis, "glm" and "glmm" currently available, "bayes" coming soon
#' @param n_bootstrap number of bootstraps
#' @param pop_size number of individuals
#' @param stratified logical wether a stratified or a pooled bootstrap should be performed
#'
#' @return returns betas of the bootstrapped analysis
#'
#' @keywords bootstrap
#'
#' @examples
#' bootstrap()
#'
#' @export

bootstrap <- function(f_species = NA, df = NA, analysis = "glm", n_bootstrap = NA, pop_size = NA, stratified = F){
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  doSNOW::registerDoSNOW(cl)

  betas <- foreach(i = 1:n_bootstrap, .combine = rbind, .packages = c("lme4", "nlme"))%dopar%{
    bootstrap = i
    #Bootstrap sampling
    if(stratified == T){
      ind.samples <- vector("list", pop_size)

      for(ind in 1:pop_size){
        selection_ind_df <- seq.int(NROW(df))[which(df$ind == ind & df$presence == 1)]
        nonselection_ind_df <- seq.int(NROW(df))[which(df$ind == ind & df$presence == 0)]
        ind.samples[[ind]] <- c(sample(nonselection_ind_df, length(nonselection_ind_df), replace = T),
                                sample(selection_ind_df, length(selection_ind_df), replace = T))
      }
      boot.sample <- unlist(ind.samples)
    }else{
      boot.sample <- sample(1:nrow(df), size = nrow(df), replace = T)
    }

    #Beta estimation
    if(analysis == "glm") return(glm(f_species, data = df[boot.sample, ], family = "binomial")$coefficients)


    if(analysis == "glmm"){
      m_glmm <- lme4::glmer(f_species, data = df[boot.sample, ], family = "binomial")
      n_preds <- length(nlme::fixed.effects(m_glmm)) - 1
      GLMM_coef <- matrix(nrow = pop_size, ncol = n_preds)
      for(i in 1:pop_size) GLMM_coef[i,] <-as.vector(as.matrix(nlme::fixed.effects(m_glmm)[-1] + nlme::random.effects(m_glmm)$ind[i,]))
      return(GLMM_coef)
    }

  }
  parallel::stopCluster(cl)
  return(betas)
}

