#function [updates, state] = Adam(gradients, state)
GD<-function(gradients,state){
if(length(state)==0){

    state$epsilon = 1e-8;
}
updates = state$epsilon*gradients
state$epsilon = 1e-8;
	return(list(updates=updates,state=state))
}
