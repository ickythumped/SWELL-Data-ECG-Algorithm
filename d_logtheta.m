function [d_logval, d_den_square] = d_logtheta(theta_vec, f_vec)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
f_vec(f_vec == 0) = 1e-2;
f_vec = log(f_vec);

theta_vec(theta_vec == 0) = 1e-2;

d_num = diff(f_vec, 1, 2);
d_num(d_num == 0) = 1e-2;

d_den = diff(theta_vec, 1, 2);
d_den(d_den == 0) = 1e-2;

d_den_square = diff(theta_vec, 2, 2);
d_den_square(d_den_square == 0) = 1e-2;
    
d_logval = d_num./d_den;


end

