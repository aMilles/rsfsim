#' Bootstrap
#'
#' This function generates bootstrapped betas.
#' @param formu formula to be used in the analysis
#' @param df data frame to be used in the analysis
#' @param analysis type of the analysis, "glm" and "glmm" currently available, "bayes" coming soon
#' @param n_bootstrap number of bootstraps
#' @param pop_size number of individuals
#' @param stratified logical wether a stratified or a pooled bootstrap should be performed
#' @param bayesianInput_ a list containing a Bayesian Setup, settings, and a character specifying the sampler created with Bayesian Tools
#'
#' @return returns betas of the bootstrapped analysis
#'
#' @keywords bootstrap
#'
#' @examples
#' bootstrap()
#'
#' @export

bootstrap <- function(formu = f.species, df = track, analysis = "glm", n_bootstrap = 100, pop_size = pop.size, stratified = F, bayesianInput_ = NA){
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  doSNOW::registerDoSNOW(cl)

  betas <- foreach(i = 1:n_bootstrap, .combine = rbind, .packages = c("lme4", "BayesianTools", "nlme"))%dopar%{
    bootstrap = i

    #Bootstrap sampling
    if(stratified == T){
      ind.samples <- vector("list", pop_size)
      for(ind in 1:pop_size){
        ind.samples[[ind]] <- sample(which(df$ind == ind), size = sum(df$ind == ind), replace = T)
      }
      boot.sample <- unlist(ind.samples)
    }else{
      boot.sample <- sample(1:nrow(df), size = nrow(df), replace = T)
    }

    #Beta estimation
    if(analysis == "glm") return(glm(formu, data = df[boot.sample, ], family = "binomial")$coefficients)


    if(analysis == "glmm"){
      m_glmm <- glmer(formu, data = df[boot.sample, ], family = "binomial")
      n_preds <- length(fixed.effects(m_glmm)) - 1
      GLMM_coef <- matrix(nrow = pop_size, ncol = n_preds)
      for(i in 1:pop_size) GLMM_coef[i,] <-as.vector(as.matrix(fixed.effects(m_glmm)[-1] + random.effects(m_glmm)$ind[i,]))

      return(GLMM_coef)
    }


    if(analysis == "bayes_nornd"){

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


    if(analysis == "bayes_rnd"){
      obs <- bayesianInput_[[1]][boot.sample]
      mm_ind <- bayesianInput_[[5]][boot.sample, ]
      ind_col <- match("ind", colnames(mm_ind))
      n_preds <- ncol(mm_ind)-2

      likelihood_rnd_Slope <- function(x, sum = TRUE){
        fix_eff <- x[1:(n_preds+1)]
        rnd_slope <- x[(n_preds+2) : (length(x)-1)]
        rnd_slope_sd <- x[length(x)]
        ind_col <- match("ind", colnames(mm_ind))
        predicted <- vector(length = nrow(mm_ind))
        for(ind in 1:pop_size){
          ind_rows <- mm_ind[, ind_col] == ind
          mixed_effect <- c(fix_eff[1], fix_eff[-1] +rnd_slope[(ind - 1) * n_preds + (1:n_preds)])
          predicted[ind_rows] <- plogis(mm_ind[ind_rows, - ind_col] %*% mixed_effect)
        }
        predicted <- plogis((mm_ind[, - ind_col] %*% fix_eff) + apply((mm_ind[, -c(1, ind_col)] * rnd_slope[(mm_ind[,ind_col] - 1) * n_preds + (1:n_preds)]), 1, sum))
        llValues <- dbinom(obs, size=1, predicted, log=T)
        llRandom <- dnorm(rnd_slope, sd = rnd_slope_sd, log = T)

        return(sum(llValues) + sum(llRandom))
      }
      rnd_names <- vector("list", length = pop_size)
      for(i in 1:pop_size) rnd_names[[i]] <- paste0(colnames(mm_ind)[-c(1, ind_col)], "_ind_", i)
      bayesianSetup_rnd_Slope <- createBayesianSetup(likelihood_rnd_Slope,
                                                     bayesianInput_[[4]],
                                                     names = c("Intercept", colnames(mm_ind)[-c(1, ind_col)], unlist(rnd_names), "rndsigma"))
      MCMC_run <- runMCMC(bayesianSetup = bayesianSetup_rnd_Slope,
                          settings = bayesianInput_[[2]],
                          sampler = bayesianInput_[[3]])

      Bayes_coef <- matrix(nrow = pop_size, ncol = n_preds)
      for(i in seq.int(pop_size)) Bayes_coef[i, ] <- MAP(MCMC_run)$parametersMAP[1:n_preds + 1] + MAP(MCMC_run)$parametersMAP[(2+i*n_preds):(n_preds+ 1 + i*n_preds)]

      return(Bayes_coef)
    }

  }

  parallel::stopCluster(cl)
  return(betas)
}

