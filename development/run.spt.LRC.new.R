rm(list=ls())
library(dplyr)
library(ggplot2)
load("development/spatialLIBD.spe.RData")

# if(inherits(try(expr_counts[["counts"]][1:5000,1:5000],silent = TRUE),'try-error')){counts=expr_counts[["counts"]][1:5000,1:5000]}else{counts=expr_counts[["counts"]][1:5000,1:5000]}
if(inherits(try(expr_counts[["counts"]][1:100,1:100],silent = TRUE),'try-error')){counts=expr_counts[["counts"]][1:100,1:100]}else{counts=expr_counts[["counts"]][1:100,1:100]}
x = counts
x = t(counts)
coords=coords[rownames(x),]
w <- fct_dist_matrix(coords)


library(sparsesvd)
sp_svd = sparsesvd(log(X+1),rank=20)

U0=sp_svd$u %*% diag(sqrt(sp_svd$d))
V0=sp_svd$v %*% diag(sqrt(sp_svd$d))

#%% Initialize Optimizer
#source("Adam.R")
#source("Nadam.R")
#source("AMSgrad.R")
#source("AdaMax.R")
#source("RMSprop.R")
#source("Adam2.R")
#source("Adadelta.R")
source("development/Adam.new.R")
source("development/Nadam.new.R")
source("development/AMSgrad.new.R")
source("development/AdaMax.new.R")
source("development/RMSprop.new.R")
source("development/Adam2.new.R")
source("development/Adadelta.new.R")
source("development/GD.new.R")
f=Adam
f=RMSprop
objs_list = vector("list",8)
names(objs_list)=c("Adam","Nadam","AMSgrad","AdaMax","RMSprop","Adam2","Adadelta","GD")
U_list=V_list = objs_list
ccc=0
# for(f in c(Adam,Nadam,AMSgrad,AdaMax,RMSprop,Adam2,Adadelta,GD)){
for(f in c(AMSgrad)){
	ccc=ccc+1
	print(names(objs_list)[ccc])
	i = 0;
	U = U0;
	V = V0;
	PP = U %*% t(V);
	EE = exp(PP);
	
	epsilon=1e-4
	max_iter=100
	objs=vector("numeric",max_iter)
	objs = rep(0,max_iter);
	obj = -Inf;
	
	flag = 0;
	#state.m_U=state.m_V=NULL
	#state.v_U=state.v_V=NULL
	state_U=state_V=NULL
	while(flag == 0 & i < max_iter){
		print(paste("The",i,"th iteration"))
		i = i + 1;
		U_old = U;
		V_old = V;
		EE_old = EE;
		TTT = t(X) - EE_old;
		
		 #   % Gradients
		update_U = TTT %*% V_old;
		update_V = t(TTT) %*% U_old;
		    
		#    % Gradient Descent Optimizer
		sgd_res_U = f(update_U,state_U)
		sgd_res_V = f(update_V,state_V)
		state_U=sgd_res_U[["state"]]
		updates_U=sgd_res_U[["updates"]]
		state_V=sgd_res_V[["state"]]
		updates_V=sgd_res_V[["updates"]]
		#    % Perform Updates
		U = U_old + updates_U;
		V = V_old + updates_V;
		PP = U %*% t(V);
		
		EE = exp(PP);
		u_diff = sum(sum(abs(U_old - U))) / sum(sum(abs(U_old)));
		v_diff = sum(sum(abs(V_old - V))) / sum(sum(abs(V_old)));
		obj_old = obj;
		XPP = PP * t(X);
		obj = sum(XPP) - sum(EE);
		objs[i] = obj;
		o_diff = abs(obj - obj_old) / abs(obj_old);
		if(u_diff < epsilon & v_diff < epsilon & o_diff < epsilon  & i > max_iter){ ##It looks like that the u_diff doesn't do much good, as the changes are already very small in the beginning....
		    flag = 1;
		}
		print(paste("The convergence of U and V, and objective value is are: ", u_diff, v_diff, obj, o_diff))
	}
	objs_list[[ccc]]=objs
	U_list[[ccc]]=U
	V_list[[ccc]]=V
}

save(objs_list,U_list,V_list, file="objs-U-V_list.RData")
load("objs-U-V_list.RData")

par(mfrow=c(2,5))
for(jj in 1:length(objs_list)){
	plot(objs_list[[jj]], main=names(objs_list)[jj])
}

plot_data <- tibble(method = "Adam", lik = objs_list$Adam, index = 1:100) %>% 
  bind_rows(tibble(method = "Nadam", lik = objs_list$Nadam, index = 1:100)) %>% 
  bind_rows(tibble(method = "AMSgrad", lik = objs_list$AMSgrad, index = 1:100)) %>% 
  bind_rows(tibble(method = "AdaMax", lik = objs_list$AdaMax, index = 1:100)) %>% 
  bind_rows(tibble(method = "RMSprop", lik = objs_list$RMSprop, index = 1:100)) %>% 
  bind_rows(tibble(method = "Adam2", lik = objs_list$Adam2, index = 1:100)) %>% 
  bind_rows(tibble(method = "Adadelta", lik = objs_list$Adadelta, index = 1:100)) %>% 
  bind_rows(tibble(method = "GD", lik = objs_list$GD, index = 1:100))

(plot_data %>% 
  ggplot() +
  geom_line(aes(x = index, y = lik, color = method))) %>% 
  plotly::ggplotly()


# USE A<S GRAD
