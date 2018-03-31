function [lambda_val] = cif(f_val, integ_f_val)
%Calculating lambda using the conditional intensity function

lambda_val = f_val./(1 - integ_f_val);
end

