#' A resource selection simulator
#'
#' This function generates a dataset usable in resource selection functions
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
#' @param steps_ number of relocations/steps
#' @param non_steps number of absences
#' @param absence_sampling default is "full" , meaning all cells are used as availability background
#' @param stepwise logical, if relocations are only dependent of the set of predictors  or also dependent of the location of the last step
#' @param dist_effect numeric defining the impact of the distance on relocation as plogis(1/dist^dist_effect)
#'
#' @return returns a dataset usable in resource selection functions and a map suitable for plotting with ggplot2
#'
#' @keywords movement
#'
#' @examples
#' move.from.preds.stepwise()
#'
#' @export
move.from.preds.stepwise <- function(formu, pop_size = pop.size, steps_ = steps, non_steps = non.steps, base_, ind_pref_mat = NA, n_preds, n_preds_original, gridsize = sqrt(nrow(predictors)), predictors_, ind_sd = NA, dist_effect = NA, stepwise = F, absence_sampling = NA, param_list = list(NA)){

  if(absence_sampling == "full" & non_steps != gridsize^2){
    non_steps = gridsize^2
    warning(paste("full absence sampling is used, non_steps is set to",gridsize^2))
    }

  #create a distance matrix to sample from later on if a stepwise approach is used
  if(stepwise == T)  dist_matrix_ext <- matrix(apply(
     X = gtools::permutations(gridsize * 3, 2,1:(gridsize * 3), repeats.allowed = T),                            MARGIN = 1,
     FUN = function(x) sqrt((x[1]-gridsize*3/2)^2 + (x[2]-(gridsize*3/2))^2)), 
    ncol = gridsize * 3)


  #initialize simulation
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  doSNOW::registerDoSNOW(cl)
  out<-foreach::foreach(input = vector("list", pop_size), i = 1:pop_size)%dopar%{
    set.seed(i)
    if(is.na(ind_pref_mat)){
      ind_pref <- base_*rnorm(mean = 1, sd = ind_sd , n = n_preds)
    }else{
      ind_pref <- ind_pref_mat[i,]
    }

    #sample presence - sample stepwise
    if(stepwise == F){
      present_cells <- sample(1:(gridsize^2), size = steps_, prob =  plogis(predictors_%*%ind_pref), replace = T)
    }else{
      initial <- sample(1:(gridsize^2), size = 1, prob =  plogis(predictors_%*%ind_pref))
      col <- ceiling(initial/gridsize)
      row <- initial - (col-1)*gridsize
      present_cells <- vector(length = steps_)
      present_cells[1] <- initial

      for(step in 2:steps_){
        dist_matrix <- dist_matrix_ext[(gridsize*1.5):(gridsize*1.5+gridsize-1) - (row - 1), (gridsize*1.5):(gridsize*1.5+gridsize-1) - (col - 1)]
        present_cells[step] <- sample(1:(gridsize^2), size = 1, prob = (plogis(1/dist_matrix^dist_effect)-0.5)*2*matrix(plogis(predictors_%*%ind_pref), ncol = gridsize))
        col <- ceiling(present_cells[step]/gridsize)
        row <- present_cells[step] - (col-1)*gridsize
      }
    }


    # sample absence - use all cells or sample a certain amount of cells where no individual was sampled
    if(absence_sampling == "not_full") absence_cells <- sample((1:(gridsize^2))[-unique(present_cells)], size = non_steps, replace = T)
    if(absence_sampling == "full") absence_cells <- 1:(gridsize^2)
      
    return(list(present_cells, absence_cells, plogis(predictors_%*%ind_pref)))
  }

  parallel::stopCluster(cl)

  # process output from parallel computing
  present_cells <- unlist(lapply(out, function(x) x[1][[1]]))
  absence_cells <- unlist(lapply(out, function(x) x[2][[1]]))
  probs <- lapply(out, function(x) x[3][[1]])

  #process further output, track contains predictors, 0/1 values, individual id and step numbers, indvidual map is used for ggplots later on
    track <- data.frame(ind = as.factor(c(rep(1:pop_size, each = steps_), rep(1:pop_size, each = non_steps))),
                      presence = c(rep(1, length(present_cells)),rep(0, pop_size * non_steps)),
                      step_nr = c(rep(1:steps_, len = length(present_cells)), rep(NA, len = length(absence_cells))))
  
  track_preds_present <- matrix(predictors_[present_cells,], ncol = n_preds)
  track_preds_absent <- matrix(predictors_[absence_cells,], ncol = n_preds)
  track_preds <- rbind(track_preds_present, track_preds_absent)
  if(n_preds == 1) track_preds <- c(track_preds_present, track_preds_absent)

  track <- data.frame(track, track_preds)
  names(track)[4:ncol(track)] <- colnames(predictors_)
  #track[, colnames(predictors_)] <- apply(track[,colnames(predictors_)], 2, scale)
  track$presence <- factor(track$presence, levels = c(0, 1))
  individual_map <- map[present_cells, ]
  individual_map <- cbind(individual_map, ind = track$ind[track$presence == 1], step_nr = rep(1:steps_, pop_size))

  return(list(individual_map, track, probs))
}
