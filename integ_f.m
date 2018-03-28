function [integ_f_val] = integ_f(j, delta, u_val, mu_val, theta_p_plus1)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

integ_f_val = integral(@(t) sym_f(t, theta_p_plus1, u_val, mu_val), u_val, j*delta);
end

