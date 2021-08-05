#function [updates, state] = Nadam(gradients, state)
Nadam<-function(gradients, state){
#%NADAM Summary of this function goes here
#%   Detailed explanation goes here

if(length(state)==0){

    state$beta1 = 0.9;
    state$beta2 = 0.999;
    state$epsilon = 1e-8;
    state$iteration = 1;
    state$m = 0*gradients;
    state$v = 0*gradients;
    state$alpha = 1e-2;

}

#% update biased first moment estimate
state$m = state$beta1 * state$m + (1 - state$beta1) * gradients;
    
#% update biased second raw moment estimate
state$v = state$beta2 * state$v + (1 - state$beta2) * gradients^2;
    
#% compute bias-corrected first moment estimate
mhat = state$m / (1 - state$beta1^(state$iteration + 1));
    
#% compute bias-corrected second raw moment estimate
vhat = state$v / (1 - state$beta2^state$iteration);

#% nadam
mhat = state$beta1 * mhat + (((1 - state$beta1) * gradients) / (1 - state$beta1^state$iteration));

#% update parameters
updates = (state$alpha * mhat) / (sqrt(vhat) + state$epsilon);

#% update iteration number
state$iteration = state$iteration + 1;


    state$beta1 = state$beta1;
    state$beta2 = state$beta2;
    state$epsilon = state$epsilon;
    state$iteration = state$iteration;
    state$m = state$m;
    state$v = state$v;
    state$alpha = state$alpha;

return(list(updates=updates, state=state))
}

