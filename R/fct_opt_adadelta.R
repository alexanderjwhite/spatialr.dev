#' Internal Adadelta Function
#'
#' @inheritParams fct_opt_amsgrad
#'
#' @return updates
#' @export
#'
#' @examples
fct_opt_adadelta <- function(gradients, state){
  #%ADADELTA optimization
  #%   Detailed explanation goes here
  
  if(length(state)==0){
    state$epsilon = 1e-6;
    state$rho = .95;
    state$iteration = 1;
    state$history = 0*gradients;
    state$u = 0*gradients;
  }
  
  #% accumulate gradient
  state$history = state$rho * state$history + (1 - state$rho) * gradients^2;
  
  #% update parameters
  updates = gradients * sqrt((state$u + state$epsilon) / (state$history + state$epsilon));
  
  #% accumulate updates
  state$u = state$rho * state$u + (1 - state$rho) * updates^2;
  
  #% update iteration number
  state$iteration = state$iteration + 1;
  
  state$epsilon = state$epsilon;
  state$rho = state$rho;
  state$iteration = state$iteration;
  state$history = state$history;
  state$u = state$u;
  
  return(list(updates=updates, state=state))
}

