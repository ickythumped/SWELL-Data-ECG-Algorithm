function [u, K, T, N, H] = init(r_peaks)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%R-R peaks in seconds
u = r_peaks(2:end);

%Number of R peaks
K = length(u); 

% Length of data [0, T]
t = 0:0.005:660;
T = length(t); 

% R-peaks vector in samples
N = 0:1:660/0.05; 
for i = 1:length(u)
    N(u(i)) = 1;
end
clear i

% History dependence vector
H = [u(1) diff(u)]; 

end

