#' A rsf score function
#'
#' This function generates rsf scores/effect plots
#'
#' @param formu  formula without an response such as ~ x1 + x2
#' @param pop_size integer defining the number of individuals
#' @param base_ vector of numerics as the parameters coefficients
#' @param pred_choice vector  of integers defining the set of predictors originally selected
#' @param n_preds integer defining total number of predictors, quadratic conditions and interactions
#' @param n_preds_original integer defining the number of predictors
#' @param gridsize integer determining size of the map
#' @param predictors_ matrix of predictors
#' @param ind_sd integer as the standard deviation of the normal distribution with a mean of 1 that is multiplied with the parameter coefficients
#' @param n_values integer as the length of the min - max - sequence of a predictor (x-axis of an effect plot)
#' @param uc_type character defining the type of uncertainty, "both", "marg" and "cond" are possible
#' @param mode character defining wether effect plots should be generated or if the function is used inside a bootstrap. "plot" and "bootstrap" are possible
#' @param param_list list of parameters for shorter lines, see code for details
#'
#' @return In plot mode it returns a data.frame with quantiles suitable for ggplot2 and the plot itself, in  bootstrap mode it returns an array with non-aggregated rsf-scores.
#'
#' @keywords rsf
#'
#' @examples
#' plot.truth()
#'
#' @export
plot.truth <- function(formu = NA, pop_size = NA, base_ = NA, pred_choice = NA, n_preds = NA, n_preds_original = NA, gridsize = NA, predictors_ = NA, ind_sd = NA, n_values = NA, uc_type = NA, mode = NA, param_list = list(NA), silent = T){


    if(!is.na(param_list)[[1]]) for(param in names(param_list[names(param_list) %in% names(as.list(environment(), all = TRUE))])) if(is.na(get(param))) assign(param, param_list[[param]])


  if(mode == "bootstrap" & uc_type == "both" | uc_type == "marg"){
    stop("marginal bootstrap currently not available")
  }

  #create matrices containing sequences of the median and the range of possible values
  cond.matrix <- matrix(ncol = n_preds_original, nrow = n_values)
  colnames(cond.matrix) <- colnames(predictors_)[pred_choice]
  med.matrix <- cond.matrix
  for(pred in pred_choice) cond.matrix[,pred] <- seq(min(predictors_[,pred]), max(predictors_[,pred]), len = n_values)
  for(pred in pred_choice) med.matrix[,pred] <- rep(median(predictors_[,pred]), len = n_values)

  #array for conditional and marginal rsf-scores
  cond.scores <- array(dim = c(n_values, pop_size))
  marg.scores <- array(dim = c(nrow(predictors_), n_values, pop_size))

  #initialize parallel computing in a for-loop to calculate marginal and conditional scores
  cl <- parallel::makeCluster(parallel::detectCores())
  doSNOW::registerDoSNOW(cl)

  if(uc_type == "both" | uc_type == "marg"){
    out <- foreach(i = 1:n_preds_original, inputs = vector("list", n_preds_original)) %dopar% {
      pred <- i
      for(ind in 1:pop_size){
        set.seed(ind)
        ind_pref <- base_*rnorm(mean = 1, sd = ind_sd, n = n_preds)

        for(marg.value in 1:n_values){
          temp.marg.matrix <- predictors_[,1:n_preds_original]
          temp.marg.matrix[, pred] <- rep(cond.matrix[marg.value, pred], nrow(predictors_))
          temp.marg.matrix <- model.matrix(formu, data.frame(temp.marg.matrix))[, -1]
          marg.scores[, marg.value, ind] <- apply(temp.marg.matrix, 1, function(x) exp(sum(ind_pref*x)))
        }
        
        temp.cond.matrix <- med.matrix
        temp.cond.matrix[, pred] <- cond.matrix[, pred]
        temp.cond.matrix <- model.matrix(formu, data.frame(temp.cond.matrix))[, -1]
        cond.scores[, ind] <- apply(temp.cond.matrix, 1, function(x) exp(sum(ind_pref*x)))
      }
      return(list(marg.scores, cond.scores))
    }
  }

  if(uc_type == "cond"){
    out <- foreach(i = 1:n_preds_original, inputs = vector("list", n_preds_original)) %dopar% {
      pred <- i
      for(ind in 1:pop_size){
        ind_pref <- base_*rnorm(mean = 1, sd = ind_sd, n = n_preds)
        temp.cond.matrix <- med.matrix
        temp.cond.matrix[, pred] <- cond.matrix[, pred]
        temp.cond.matrix <- model.matrix(formu, data.frame(temp.cond.matrix))[, -1]
        cond.scores[, ind] <- apply(temp.cond.matrix, 1, function(x) exp(sum(ind_pref*x)))
      }
      return(list(marg.scores, cond.scores))
    }
  }
  parallel::stopCluster(cl)

  if(mode == "plot"){
    if(uc_type == "both" | uc_type == "marg"){
      scores3 <- matrix(nrow = 3, ncol = n_values*n_preds_original)
      for(i in 1:n_preds_original) scores3[,(1+(i-1)*n_values):(i*n_values)]<- apply(out[[i]][[1]], 2, quantile, c(0.025, 0.5, 0.975))
      scores4 <- cbind(data.frame(t(scores3)), score = rep(colnames(med.matrix), each = n_values), pred = as.vector(cond.matrix))
      names(scores4)[1:3] <- c("lwr", "med", "upr")
      marg.scores5 <- reshape::melt(scores4, id.vars = c("score", "pred"))
    }

    if(mode == "plot")
    if(uc_type == "both" | uc_type == "cond"){
      scores3 <- matrix(nrow = 3, ncol = n_values*n_preds_original)
      for(i in 1:n_preds_original) scores3[,(1+(i-1)*n_values):(i*n_values)]<- apply(out[[i]][[2]], 1, quantile, c(0.025, 0.5, 0.975))
      scores4 <- cbind(data.frame(t(scores3)), score = rep(colnames(med.matrix), each = n_values), pred = as.vector(cond.matrix))
      names(scores4)[1:3] <- c("lwr", "med", "upr")
      cond.scores5 <- reshape::melt(scores4, id.vars = c("score", "pred"))
    }

    gg.condmarg <-ggplot()+
      theme_minimal()+
      scale_linetype_manual(values = c("dashed", "solid", "dashed"))+
      guides(color = guide_legend(title = "model"), linetype = guide_legend(title = "CI"))


    if(uc_type == "both"){
      gg.condmarg <- gg.condmarg +
      geom_line(data = cond.scores5, aes(x = pred, y = value, linetype = variable, col = "cond_truth"))+
      geom_line(data = marg.scores5, aes(x = pred, y = value, linetype = variable, col = "marg_truth"))+
      facet_wrap(~score)

      if(!silent) beepr::beep(1)
      return(list(gg.condmarg, cond.scores5, marg.scores5))
    }

    if(uc_type == "cond"){
      gg.condmarg <- gg.condmarg +
        geom_line(data = cond.scores5, aes(x = pred, y = value, linetype = variable, col = "cond_truth"))+
        facet_wrap(~score)

      if(!silent) beepr::beep(1)
      return(list(gg.condmarg, cond.scores5))
    }

    if(uc_type == "marg"){
      gg.condmarg <- gg.condmarg +
        geom_line(data = marg.scores5, aes(x = pred, y = value, linetype = variable, col = "marg_truth"))+
        facet_wrap(~score)

      if(!silent) beepr::beep(1)
      return(list(gg.condmarg, marg.scores5))
    }
  }

  if(mode == "bootstrap"){
    if(!silent) beepr::beep(1)
    if(uc_type == "cond") return(lapply(out, function(x) x[[2]]))
    if(uc_type == "marg") return(lapply(data, function(x) x[[1]]))
    if(uc_type == "both") return(out)
  }
}
