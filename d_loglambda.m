function [d_loglambda_vec] = d_loglambda(lambda_val, d_lambda_vec)
% Calculating 1st derivative of log(lambda) with respect to theta 

d_loglambda_vec = d_lambda_vec./lambda_val;
end

