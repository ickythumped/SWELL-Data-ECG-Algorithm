function [J, delta, n] = del_parts(T, fs, N)
%Dividing [0,T] into J intervals

%% Computing J, delta
%Divide [0,T] into J equal parts
J = floor(T/(fs/128)); 

%Length of each part
delta_samples = T/J; 

%delta in seconds
delta = delta_samples/fs;

%% Computing n_j

% n_j -> binary value for peak in each delta
n = zeros(1,J); 
for j = 1:J
    for c = ((j-1)*delta_samples)+1:j*delta_samples
        if (N(c) == 1)
            n(j) = 1;
            break;
        end
    end
    clear c
end
clear j

end

