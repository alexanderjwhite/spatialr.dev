fct_v_grad <- function(x, u, v){
  t(t(x)-exp(u%*%t(v)))%*%u
}
