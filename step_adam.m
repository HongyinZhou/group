function [grad_adam, mopt, vopt] = step_adam(gradient, mopt_old, vopt_old, iteration)
%Performs one step of adam optimization
%   此处显示详细说明
beta1=0.9; beta2=0.999;
epsilon=1e-8;
    mopt = beta1 * mopt_old + (1 - beta1) * gradient;
    mopt_t = 1/ (1 - (beta1)^(iteration + 1))*mopt;
    vopt = beta2 * vopt_old + (1 - beta2) * (gradient.^2);
    vopt_t =1 / (1 - (beta2)^(iteration + 1))* vopt;
    
    grad_adam=mopt_t./(sqrt(vopt_t)+epsilon); 
end

