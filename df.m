function [df_vec] = df(j, k, nparams, delta, u_val, H, mu_val, f_val, theta_p_plus1)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
df_vec = zeros(1, nparams+2);

expr1 = (f_val*theta_p_plus1)./(mu_val^3);
expr2 = (j*delta - u_val - mu_val);
for i = 1:nparams+1
    df_vec(i) = H(k-i+2)*expr1*expr2;
end

expr3 = (-(expr2)^2)./(2*(mu_val^2)*(j*delta - u_val));
expr4 = f_val./(2*theta_p_plus1);
df_vec(nparams+2) = f_val*expr3*expr4;
end

