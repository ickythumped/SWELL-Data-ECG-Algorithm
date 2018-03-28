function [integ_f_val] = integ_f(j, delta, u_val, sym_f_val)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
sym_f_val = double(sym_f_val);
func_handle = @(t) sym_f_val;
integ_f_val = integral(func_handle, u_val, j*delta);
end

