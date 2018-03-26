function [f_val] = f(j, k, theta_p_plus1, delta, u, mu)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

lambda_idg = theta_p_plus1;
x_idg = j*delta - u(k);
mu_idg = mu(j);
f_val = mean(pdf('InverseGaussian', x_idg, 'mu', mu_idg, 'lambda', lambda_idg));

