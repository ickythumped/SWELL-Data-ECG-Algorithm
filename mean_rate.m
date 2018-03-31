function [mean_val] = mean_rate(k, history_vec, num_params, theta_j_vec)
%Calculating mean heart rate

mean_val = theta_j_vec(1);
for i = 2:num_params+1
    if k - i + 2 > 0
        mean_val = mean_val + theta_j_vec(i).*history_vec(k-i+2);
    end
end
clear i
end

