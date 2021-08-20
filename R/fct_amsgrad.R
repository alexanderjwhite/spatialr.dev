#' Internal AMSGrad Function
#'
#' @param gradients matrix; gradients
#' @param state object; state information
#'
#' @return updates
#' @export
#'
#' @examples fct_amsgrad(matrix(rnorm(20), nrow = 10, ncol = 2), state = NULL)
fct_amsgrad <- function(gradients, state){
  
  if(length(state)==0){
    
    state$beta1 = 0.9;
    state$beta2 = 0.999;
    state$epsilon = 1e-8;
    state$iteration = 1;
    state$m = 0*gradients;
    state$v = 0*gradients;
    state$vhat = 0*gradients;
    state$alpha = 1e-2;
  }
  
  #% update biased first moment estimate
  state$m = state$beta1 * state$m + (1 - state$beta1) * gradients;
  
  #% update biased second raw moment estimate
  state$v = state$beta2 * state$v + (1 - state$beta2) * gradients^2;
  
  #% non-decreasing
  state$vhat = pmax(state$vhat, state$v);
  
  #% update parameters
  updates = state$alpha * state$m / (sqrt(state$vhat) + state$epsilon);
  
  #% update iteration number
  state$iteration = state$iteration + 1;
  
  
  state$beta1 = state$beta1;
  state$beta2 = state$beta2;
  state$epsilon = state$epsilon;
  state$iteration = state$iteration;
  state$m = state$m;
  state$v = state$v;
  state$vhat = state$vhat;
  state$alpha = state$alpha ;
  
  return(list(updates=updates, state=state))
  
}
