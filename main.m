clear all; %#ok<CLALL>
close all;
clc;
[totalData1, HR_data1, ~] = myDoReadData("D:\4th sem\Physiological Signal Processing\SWELL Dataset\Data\pp1_18-9-2012_c1.S00");

%% Downsampling
r = 1;
fs = totalData1.fs/r;
%HR_data = decimate(HR_data1, r);
HR_data = HR_data1(1 : (floor(length(HR_data1)/fs))*fs);

%% Locating and Plotting R-R peaks
% U -> Location of R-peaks in samples
% For accurate identification of peaks please select appropriate value for
% ...MinPeakHeight & MinPeakDistance
MinPeakHeight = 500;
MinPeakDistance = 1/4;
[~, U] = peakfinder(HR_data, fs, MinPeakHeight, MinPeakDistance);

%% Declarations
% u -> Location of R-peaks in seconds
% K -> Number of R-peaks
% T -> Length of data [0, T]
% N -> R-peaks vector in samples
% H -> History vector
[u, K, T, N, H] = init(HR_data, fs, U);
syms t; %variable for computing integ_f & integ_df

% R-R plot
figure
scatter(u, H, '*' , 'r')
xlabel('time in seconds')
ylabel('R-R intervals')
title('R-R Plot')
%% j division
% J -> Divide [0,T] into J equal parts
% delta -> Length of each part in seconds
% n_j -> binary value for peak in each delta
[J, delta, n] = del_parts(T, fs, N);

%% Parameter initialization
nparams = 2; %Number of parameters in theta

%wait for 'p' previous spikes
k = 0;
start_iter = 1;
while(k ~= nparams)
    if(u(k+1) <= (start_iter-1)*delta)
        k = k+1;
    end
    start_iter = start_iter+1;
end

theta_predict = zeros(nparams + 2, J); %Model parameter vector (j|j-1)
Wvar_predict = ones(nparams+2, nparams+2); %Model parameter vector
Wvar_update = zeros(nparams+2, nparams+2); %Model parameter vector
mu = zeros(1, J); %mean of each interval
sigma_square = zeros(1, J); %variance of each interval
lambda_vec = zeros(1, J); % lambda vector
integ_term = zeros(1, J);
f_vec = zeros(1, J); %f() vector

%% Initializations (Note: work to be done)
theta_update = [0.834; -0.15; -0.25; 0.396]; %intitalizing theta(j|j)
theta_predict(:, start_iter) = theta_update;
mu(start_iter) = mean_rate(k, H, nparams, theta_update);
sigma_square(start_iter) = 0.75;
theta_predict(nparams+2, start_iter) = mu(start_iter)^3./sigma_square(start_iter)^2;
theta_update(nparams+2) = mu(start_iter)^3./sigma_square(start_iter)^2;
covar_matrix = diag(theta_update);

% del initializations

%% Algorithm without adaptive filter

for j = start_iter+1:J
   if(u(k+1) <= (j-1)*delta)
       k = k+1;
   end
   
   % Compute mean
   mu(j) = mean_rate(k, H, nparams, theta_update);
   
   % Compute f()
   f_vec(j) = f(j, theta_update(end), delta, u(k), mu(j));
   
   % Compute sym_f() and integ_f()
   f_eq = sym_f(t, theta_update(end), u(k), mu(j));
   integ_f_val = integ_f(j, delta, u_(k), f_eq); %have to complete this
   
   % Compute cif
   lambda_vec(j) = cif(f_vec(j), integ_f_val);
   
   % Compute del
   df_vec = df(j, k, nparams, delta, u(j), H, mu(j), f_vec(j), theta_update(end));
   integ_df_vec = integ_df(df_vec); %have to code this
   
   d_lambda_vec = d_lambda(j, nparams, integ_f_val, df_vec, integ_df_vec, f_vec);
   
   d_loglambda_vec = d_loglambda(lambda_val, d_lambda_vec);
   
   % Compute del square
   
end

%% Algorithm - with adaptive filter

% for j = start_iter+1:J
%    if(u(k+1) <= (j-1)*delta)
%        k = k+1;
%    end
%    
%    % Prediction Step
%    theta_predict(:, j)  = theta_update;
%    Wvar_predict = Wvar_update + covar_matrix;
%   
%    % Compute mean
%    mu(j) = mean_rate(k, H, nparams, theta_predict(:, j));
%    
%    % Compute f()
%    f_idg(j) = f(j, theta_predict(nparams+2, j), delta, u(k), mu(j));
%    
%    % Compute lambda
%    [lambda_j(:, j), integ_term] = cif(j,k, u, theta_predict(nparams+2, :), delta, f_idg(j), mu);
%    
%    % Compute derivates wrt theta
%    
%
%    %update posterior parameters
%    theta_update = theta_predict(:, j) + (Wvar_predict * del_loglambda(:,end) .* (n(j) - (lambda_j(j)*delta)));
%    
%    Wvar_update = inv(inv(Wvar_predict) - delsquare_loglambda.*(n(j) - (lambda_j(j)*delta))...
%        - del_loglambda(:,end) * (del_lambda.*delta).'); 
%     
% end

