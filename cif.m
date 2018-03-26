function [lambda_val, integ_term] = cif(j, k, u, theta_p_plus1, delta, cif_val, mu)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%with adaptive filter
% u_limit = cdf('InverseGaussian', j*delta, mu(j), theta_p_plus1(j));

%without adaptive filter
u_limit = cdf('InverseGaussian', j*delta, mu(j), theta_p_plus1(end));


lastspike_bin = ceil(u(k)/delta);

%with adaptive filter
% l_limit = cdf('InverseGaussian', lastspike_bin*delta,...
%         mu(lastspike_bin), theta_p_plus1(lastspike_bin));

%without adaptive filter
l_limit = cdf('InverseGaussian', lastspike_bin*delta,...
        mu(lastspike_bin), theta_p_plus1(end));

integ_term = u_limit - l_limit;
lambda_val = cif_val/(1 - integ_term);
end

