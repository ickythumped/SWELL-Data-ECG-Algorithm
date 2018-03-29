function [sym_df_val] = sym_df(i, t, k, nparams,...
    u_val, H, mu_val, f_val, theta_p_plus1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

expr1 = (f_val*theta_p_plus1)./(mu_val^3);
expr2 = (t - u_val - mu_val);
    
expr3 = (-(expr2).^2)./(2*(mu_val^2).*(t - u_val));
expr4 = f_val./(2*theta_p_plus1);

if (i~=nparams+2)
    sym_df_val = H(k-i+2).*expr1.*expr2;
else
    sym_df_val = f_val.*expr3 + expr4;
end

end

