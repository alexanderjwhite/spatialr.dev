#function [updates, state] = Adam2(gradients, state)
Adam2<-function(gradients, state){
#%ADAM2 Summary of this function goes here
#%   Detailed explanation goes here
if(length(state)==0){
    state$beta1 = 0.9;
    state$beta2 = 0.999;
    state$beta1t = state$beta1;
    state$beta2t = state$beta2;
    state$epsilon = 1e-8;
    state$m = 0*gradients;
    state$v = 0*gradients;
    state$vhat = state$v;
    state$alpha = 1e-2;

}

#% update biased first moment estimate
state$m = state$beta1 * state$m + (1 - state$beta1) * gradients;

#% update biased second raw moment estimate
state$v = state$beta2 * state$v + (1 - state$beta2) * gradients^2;

#% bias correction
bc = sqrt(1 - state$beta2t) / (1 - state$beta1t);

#% compute bias-corrected second raw moment estimate
state$vhat = pmax(state$vhat, state$v);

#% update parameters
updates = state$alpha * state$m * bc / (sqrt(state$vhat) + state$epsilon);
state$beta1t = state$beta1 * state$beta1t;
state$beta2t = state$beta2 * state$beta2t;

    state$beta1 = state$beta1;
    state$beta2 = state$beta2;
    state$beta1t = state$beta1t;
    state$beta2t = state$beta2t;
    state$epsilon = state$epsilon;
    state$m = state$m;
    state$v = state$v;
    state$vhat = state$vhat;
    state$alpha = state$alpha;
return(list(updates=updates, state=state))
}
