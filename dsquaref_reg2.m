function [dsquaref_reg2_vec] = dsquaref_reg2(j, k, delta, u_val,...
    H, mu_val, f_val, df_vec, theta_p_plus1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dsquaref_reg2_vec = zeros(1:nparams+1);
expr1 = j*delta - u_val;
expr2 = ((expr1 - u_val)^2)./(2*(mu^2)*expr1);
expr3 = (f_val*(expr1 - mu_val))./mu_val^3;

for i = 1:nparams+1
    dsquaref_reg2_vec(i) = ((1/(2*theta_p_plus1))*df_vec(i)) - (df_vec(i)*expr2) + (expr3*H(k-i+2)); 
        
end
end

