function [integ_dfsquare_mat] = integ_dfsquare(j, delta, k, nparams, u_val, H,...
    mu_val, f_val, df_vec, theta_p_plus1)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

integ_dfsquare_mat = zeros(nparams+2, nparams+2);
for a = 1:nparams+2
    for b = 1:nparams+2
        f_handle = @(t) sym_dsquaref(a, b, t, k, nparams, u_val, H,...
        mu_val, f_val, df_vec, theta_p_plus1);
    
        integ_dfsquare_mat(a, b) = integral(f_handle, (j-1)*delta, j*delta);
    end
end
end

