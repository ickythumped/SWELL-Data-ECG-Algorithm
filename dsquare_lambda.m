function [dsquare_lambda_mat] = dsquare_lambda(nparams, f_val, integ_f_val, ...
    d_lambda_vec, df_vec, integ_df_vec, dsquaref_mat, integ_dfsquare_mat) 
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

dsquare_lambda_mat = zeros(nparams+2, nparams+2);
expr1 = 1/(1 - integ_f_val);

for a = 1:nparams+2
    for b = 1:nparams+2
        term1 = (2.*expr1).*(d_lambda_vec(b).*integ_df_vec(a));
        term2 = dsquaref_mat(a,b).*expr1;
        term3 = (expr1^2).*df_vec(b).*integ_df_vec(a);
        term4 = (expr1^2).*df_vec(a).*integ_df_vec(b);
        term5 = f_val.*(expr1^2).*integ_dfsquare_mat(a,b);
        
        dsquare_lambda_mat(a, b) = term1 + term2 - term3 + term4 + term5;
    end
end


end

