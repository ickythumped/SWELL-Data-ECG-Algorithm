function [mean_val] = mean_rate(k, history_rate, num_params, theta_vec)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

mean_val = theta_vec(1);
for i = 2:num_params+1
    if k - i + 2 > 0
        mean_val = mean_val + theta_vec(i).*history_rate(k-i+2);
    end
end
clear i
end

