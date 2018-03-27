function [d_lambda_vec] = d_lambda(j, nparams, integ_f_val, df_vec, integ_df_vec, f_vec)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

expr1 = (1 - integ_f_val);
d_lambda_vec = zeros(1, nparams+2);
for i = 1:nparams+2
    d_lambda_vec(i) = (expr1*df_vec(i) +  f_vec(j)*integ_df_vec(i))./(expr1^2);
end

end

