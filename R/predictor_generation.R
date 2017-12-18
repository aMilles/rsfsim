#'Predictor Generation
#'
#' Generate predictors_ from a landscape
#'
#' @param formu  formula without an response such as ~ x1 + x2
#' @param base vector of numerics as the parameters coefficients
#' @param pred_choice vector  of integers defining the set of predictors originally selected
#' @param n_preds integer defining total number of predictors, quadratic conditions and interactions
#' @param n_preds_original integer defining the number of predictors
#' @param gridsize integer determining size of the map
#' @param code code (e.g. "110") to be used in simData, check FReibier::simData documentation for further help
#' @param effect.size max distance of the base coefficients to zero
#'
#' @keywords predictors
#'
#' @return matrix of predicotrs and a map suitable for plotting with ggplot2
#'
#' @examples
#' preds.from.ls()
#'
#' @export
preds.from.ls <- function(pred_choice, gridsize, code, n_preds_original, n_preds, base_ = NA, formu, effect.size = 1){
  set.seed(123)
  if(is.na(base_)) base_ <- runif(min = -effect.size, max = effect.size, n = n_preds)

  if(!paste0(code, gridsize, ".nc")[1] %in% list.files()){
    #source("simData.R")
    FReibier::simData(code, filename = paste0(code, gridsize)[1], gridsize = rep(gridsize, 2))
  }

  ls <- ncdf4::nc_open(paste0(code, gridsize, ".nc")[1])

  pred_choice = c(2, 3, pred_choice + 4) # add 2, 3 to pred_choice for Lon, Lat
  map <- matrix(ncol = length(names(ls$var)[pred_choice]), nrow = gridsize^2)
  colnames(map) <- names(ls$var)[pred_choice]
  for(i in names(ls$var)[pred_choice]) map[,i] <- ncvar_get(ls, varid = i)
  map <- melt(data.frame(map), id.vars = c("Lat", "Lon"))
  predictors_ <- matrix(map[,4], byrow = F, ncol = n_preds_original)
  colnames(predictors_) <- levels(map$variable)
  predictors_ <- model.matrix(formu, data.frame(predictors_))
  pred.names <- colnames(predictors_)[-1]
  predictors_ <- matrix(predictors_[,-1], ncol = n.preds)
  colnames(predictors_) <- pred.names 
  
  probs <- plogis(predictors_ %*%base_)
  
  map$probs <- rep(probs, nrow(map)/gridsize^2)

  return(list(predictors_, map, base_))
}

# base_ = NA
# pred_choice = pred.choice
# gridsize = gridsize
# code = code
# n_preds_original = n.preds.original
# n_preds = n.preds
# formu = f.species
