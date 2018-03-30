function [u, K, T, N, H] = init(r_peaks)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%R-R peaks in seconds
u = r_peaks;

%Number of R peaks
K = length(u); 

% Length of data [0, T]
t = 1:0.005:160;
T = length(t); 

% R-peaks vector in samples
N = zeros(1,T); 
for i = 1:length(u)
    N(u(i)) = 1;
end
clear i

% History dependence vector
H = [u(1) diff(u)]; 

end

