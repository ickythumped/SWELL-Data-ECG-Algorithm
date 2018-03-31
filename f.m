function [f_val] = f(j, theta_p_plus1, delta, u_val, mu_val)
%Calculating f() for History dependent Inverse Gaussian Distribution

%% PDF approach
% lambda_idg = theta_p_plus1;
% x_idg = j*delta - u_val;
% f_val = mean(pdf('InverseGaussian', x_idg, 'mu', mu_val, 'lambda', lambda_idg));


%% Traditional Approach
expr1 = theta_p_plus1;
expr2 = j*delta - u_val;
expr3 = exp((-0.5*expr1*(expr2 - mu_val)^2)./(mu_val^2 * expr2));

f_val = sqrt(expr1./(2*pi*(expr2)^3)) * expr3;
end

