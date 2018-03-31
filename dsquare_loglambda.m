function [dsquare_loglambda_mat] = dsquare_loglambda(nparams, lambda_val, dlambda_vec,...
    dsquare_lambda_mat)
% Calculating 2nd derivative of log(lambda) with respect to theta 

dsquare_loglambda_mat = zeros(nparams+2, nparams+2);
for i = 1:nparams+2
    for j = 1:nparams+2
        dsquare_loglambda_mat(i, j) = dsquare_lambda_mat(i,j).*(1/lambda_val)...
            - dlambda_vec(i).*dlambda_vec(j).*(1./(lambda_val)^2);
    end
end
end

