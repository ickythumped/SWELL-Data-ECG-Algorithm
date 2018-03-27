function [dsquaref_reg3_val] = dsquaref_reg3(j, delta, f_val, df_vec, theta_p_plus1)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

expr1 = j*delta - u_val;
expr2 = ((expr1 - u_val)^2)./(2*(mu^2)*expr1);
expr3 = df_vec(end) - f_val./theta_p_plus1;

dsquaref_reg3_val = (df_vec(end)*(-expr2)) + ((1/(2*theta_p_plus1))*expr3);
end

