function [dsquare_val] = dsquare_theta(del_lambda, diffsquare_theta)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
dsquare_val = zeros(length(diffsquare_theta));
w = 1;
for z = 1:numel(diffsquare_theta)
    for y = 1:numel(del_lambda)
        dsquare_val(w) = del_lambda(y)/diffsquare_theta(z);
        w = w+1;
    end
end

end

