function [d_loglambda_vec] = d_loglambda(lambda_val, d_lambda_vec)
%Calculate d_loglambda

d_loglambda_vec = d_lambda_vec./lambda_val;
end

