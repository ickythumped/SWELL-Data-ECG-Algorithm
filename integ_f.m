function [integ_f_val] = integ_f(j, delta, u_val, mu_val, theta_p_plus1, integ_f_value)
%Calculating integral of f()

f_handle = @(t) sym_f(t, theta_p_plus1, u_val, mu_val);

if(integ_f_value == 0)
    integ_f_val = integral(f_handle, u_val, j*delta);
else
    integ_f_val = integral(f_handle, (j-1)*delta, j*delta);
    
%integ_f_val = integral(@(t) sym_f(t, theta_p_plus1, u_val, mu_val), u_val, j*delta);
end

