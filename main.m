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
lambda_j = zeros(1, J); % lambda vector
integ_term = zeros(1, J);
f_idg = zeros(1, J); %f() vector

%% Initializations (Note: work to be done)
theta_update = [0.834; -0.15; -0.25; 0.396]; %intitalizing theta(j|j)
theta_predict(:, start_iter) = theta_update;
mu(start_iter) = mean_rate(k, H, nparams, theta_update);
sigma_square(start_iter) = 0.75;
theta_predict(nparams+2, start_iter) = mu(start_iter)^3./sigma_square(start_iter)^2;
theta_update(nparams+2) = mu(start_iter)^3./sigma_square(start_iter)^2;
covar_matrix = diag(theta_update);

%% Algorithm without adaptive filter

for j = start_iter+1:J
   if(u(k+1) <= (j-1)*delta)
       k = k+1;
   end   
   % Compute mean
   mu(j) = mean_rate(k, H, nparams, theta_update);
   % Compute f()
   f_idg(j) = f(j, k, theta_update, delta, u, mu);
   % Compute lambda
   [lambda_j(j), integ_term(j)] = cif(j,k, u, theta_update, delta, f_idg(j), mu, integ_term(j-1));
end

%% Algorithm - with adptive filter

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
%    f_idg(j) = f(j, k, theta_predict(nparams+2, j), delta, u, mu);
%    
%    % Compute lambda
%    [lambda_j(:, j), integ_term] = cif(j,k, u, theta_predict(nparams+2, :), delta, f_idg(j), mu);
%    
%    % Compute derivates wrt theta
%    
%    [del_loglambda, diffsquare_theta] = d_logtheta(theta_predict(:, j-2:j),...
%        lambda_j(:, j-2:j)); %d(log(lambda))/d(theta)
%     
%    delsquare_loglambda = dsquare_theta(del_loglambda(:, 2), diffsquare_theta);
%    
%    del_lambda = d_theta(theta_predict(:, j-1:j), lambda_j(:, j-1:j));
%    
%    %update posterior parameters
%    theta_update = theta_predict(:, j) + (Wvar_predict * del_loglambda(:,end) .* (n(j) - (lambda_j(j)*delta)));
%    
%    Wvar_update = inv(inv(Wvar_predict) - delsquare_loglambda.*(n(j) - (lambda_j(j)*delta))...
%        - del_loglambda(:,end) * (del_lambda.*delta).'); 
%     
% end

