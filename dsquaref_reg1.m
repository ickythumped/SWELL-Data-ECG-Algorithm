function [dsquaref_reg1_mat] = dsquaref_reg1(j, k, delta, nparams, u_val,...
    H, mu_val, f_val, theta_p_plus1)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
dsquaref_reg1_mat = zeros(nparams+1, nparams+1);
expr1 = j*delta - u_val;
expr2 = (expr1 - mu_val)./(mu_val^3);
expr3 = (2*mu_val - 3*expr1)./(mu_val^4);
df_dmu = (f_val*theta_p_plus1*(expr1 - mu_val))./(mu_val^3);

dsquaref_dmusquare = theta_p_plus1.*(df_dmu*expr2 + f_val*expr3);

for a = 1:nparams+1
    for b = 1:nparams+1
        dsquaref_reg1_mat(a,b) = dsquaref_dmusquare*H(k-a+2)*H(k-b+2);
    end
end
end

