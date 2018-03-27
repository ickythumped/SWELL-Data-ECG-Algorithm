function [f_eq] = sym_f(t, theta_p_plus1, u_val, mu_val)
%Symbolic equation for f() in order to compute integ_f()

expr1 = theta_p_plus1;
expr2 = t - u_val;
expr3 = exp((-0.5*expr1*(expr2 - mu_val)^2)./(mu_val^2 * expr2));

f_eq = sqrt(expr1./(2*pi*(expr2)^3)) * expr3;


end

