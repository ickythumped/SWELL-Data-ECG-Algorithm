function [f_val] = f(j, theta_p_plus1, delta, u_val, mu_val)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

lambda_idg = theta_p_plus1;
x_idg = j*delta - u_val;
f_val = mean(pdf('InverseGaussian', x_idg, 'mu', mu_val, 'lambda', lambda_idg));

