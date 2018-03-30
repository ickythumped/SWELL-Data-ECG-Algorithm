function [J, delta, n] = del_parts(T, u)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Computing J, delta
%Divide [0,T] into J equal parts
delta = 0.005;
J = floor(T/delta); 

%% Computing n_j

% n_j -> binary value for peak in each delta
n = zeros(1,J); 
for j = 1:J
    if(u(1) > (j-1)*delta && u(1) <= j*delta)
        n(j) = 1;
        u = circshift(u,1,2); 
    end
end
clear j

end

