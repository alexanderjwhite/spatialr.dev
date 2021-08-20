fct_u_grad <- function(x, u, v, w, j, lambda){
  (t(x)-exp(u%*%t(v)))%*%v-2*lambda*(sum(diag(j%*%t(w)%*%u))-t(w)%*%u - w%*%u)
}
