#function [updates, state] = AdaMax(gradients, state)
AdaMax<-function(gradients, state){
#%ADAMAX Summary of this function goes here
#%   Detailed explanation goes here
if(length(state)==0){
    state$beta1 = 0.9;
    state$beta2 = 0.999;
    state$epsilon = 1e-8;
    state$iteration = 1;
    state$m = 0*gradients;
    state$u = 0*gradients;
    state$alpha = 1e-2;

}

#% update biased first moment estimate
state$m = state$beta1 * state$m + (1 - state$beta1) * gradients;
    
#% update biased second raw moment estimate
state$u = pmax(state$beta2 * state$u, abs(gradients));
    
#% compute bias-corrected first moment estimate
mhat = state$m / (1 - (state$beta1)^(state$iteration));
    
#% update parameters
updates = state$alpha * mhat / (state$u + state$epsilon);

#% update iteration number
state$iteration = state$iteration + 1;
return(list(updates=updates, state=state))

}

