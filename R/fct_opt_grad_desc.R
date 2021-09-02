#' Internal Gradient Descent Function
#'
#' @inheritParams fct_opt_amsgrad
#'
#' @return updates
#' @export
#'
#' @examples fct_opt_grad_desc(matrix(rnorm(20), nrow = 10, ncol = 2), state = NULL)
fct_opt_grad_desc<-function(gradients,state){
  if(length(state)==0){
    
    state$epsilon = 1e-8;
  }
  updates = state$epsilon*gradients
  state$epsilon = 1e-8;
  return(list(updates=updates,state=state))
}
