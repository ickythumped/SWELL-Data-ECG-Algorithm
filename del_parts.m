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
    for k = 1:length(u)
        if (u(k) >
end
clear j

end

