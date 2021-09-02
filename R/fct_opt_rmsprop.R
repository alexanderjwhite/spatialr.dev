#' Internal RMSprop Function
#'
#' @inheritParams fct_opt_amsgrad
#'
#' @return updates
#' @export
#'
#' @examples fct_opt_rmsprop(matrix(rnorm(20), nrow = 10, ncol = 2), state = NULL)
fct_opt_rmsprop <- function(gradients, state){
  #%RMSPROP rmsprop optimization
  #%   Detailed explanation goes here
  
  if(length(state)==0){
    state$alpha = 1e-3;
    state$rho = 0.9;
    state$epsilon = 1e-8;
    state$iteration = 1;
    state$history = 0*gradients;
  }
  
  state$history = state$rho * state$history + (1 - state$rho) * gradients^2;
  
  #% update parameters
  updates = gradients * state$alpha / sqrt(state$history + state$epsilon);
  
  #% update iteration number
  state$iteration = state$iteration + 1;
  
  state$alpha = state$alpha;
  state$rho = state$rho;
  state$epsilon = state$epsilon;
  state$iteration = state$iteration ;
  state$history = state$history;
  return(list(updates=updates, state=state))
}

