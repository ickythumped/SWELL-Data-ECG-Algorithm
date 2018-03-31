function [integ_df_vec] = integ_df(j, k, delta, nparams,...
    u_val, H, mu_val, f_val, theta_p_plus1)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

integ_df_vec = zeros(1, nparams+2);
for i = 1:nparams+2
    f_handle = @(t) sym_df(i, t, k, nparams,...
        u_val, H, mu_val, f_val, theta_p_plus1);
    
    integ_df_vec(i) = integral(f_handle, (j-1)*delta, j*delta);
end

end

