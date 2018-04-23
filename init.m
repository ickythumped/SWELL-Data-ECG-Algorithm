function [u, K, T, N, H] = init(HR_data, fs, U)
%Initializing required vectors

%R-R peaks in seconds
u = U./fs;

%Number of R peaks
K = length(U); 

% Length of data [0, T]
T = length(HR_data); 

% R-peaks vector in samples
N = zeros(1,length(HR_data)); 
for i = 1:length(U)
    N(U(i)) = 1;
end
clear i

% History dependence vector
H = [u(1) diff(u)]; 
%H(H > 1.5) = 1;

end

