%%Run this file.

clear all; %#ok<CLALL>
close all;
clc;

%Set the address of file
[totalData1, HR_data1, ~] = myDoReadData...
    ("D:\4th sem\Physiological Signal Processing\SWELL Dataset\Data\pp2_c1.S00");

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
% syms t; %variable for computing integ_f & integ_df

%% j division
% J -> Divide [0,T] into J equal parts
% delta -> Length of each part in seconds
% n_j -> binary value for peak in each delta
[J, delta, n] = del_parts(T, fs, N);

%% Parameter initialization
nparams = 9; %Number of parameters in theta

%wait for 'p' previous spikes
k = 0;
start_iter = 0;
while(k ~= nparams)
    start_iter = start_iter+1;
    if(u(k+1) <= (start_iter)*delta)
        k = k+1;
    end
end

theta_predict = zeros(nparams + 2, J); %Model parameter vector (j|j-1)
mu = zeros(1, J); %mean of each interval
std_vec = zeros(1, J); %variance of each interval

theta_update = [0.14937830871589303; 1.1564059506167039; -0.39835826229059518; ...
    0.83088601289788633; -1.023213140764317; 0.41452982839565033; -0.26784526179216739; ...
    0.6149154253134137; -0.60566867990517192; 0.11520209772910678; 1507.727429652082];
theta_predict(:, start_iter) = theta_update;
mu(start_iter) = mean_rate(k, H, nparams, theta_update);


%std_vec(start_iter) = 0.04;
covar = [3e-7; 4e-13; 4e-13; 4e-13; 4e-13; 4e-13; 4e-13; 4e-13; 4e-13; 4e-13; 0.0039];
covar_matrix = diag(covar);
Wvar_predict = ones(nparams+2, nparams+2); %Model parameter vector
Wvar_update = covar_matrix; %Model parameter vector

% Integral value initializations
integ_f_value = 0;
integ_df_vector = zeros(1, nparams+2);
integ_dfsquare_matrix = zeros(nparams+2, nparams+2);

% Paramter initializations
lambda_vec = zeros(1, J); % lambda vector
f_vec = zeros(1, J); %f() vector

%% Algorithm with adaptive filter
for j = start_iter:J
  
   if(k == length(u))
       continue
   end
   % Check for previous spike
   if(u(k+1) < j*delta)
       k = k+1;
       integ_f_value = 0;
       integ_df_vector = zeros(1, nparams+2);
       integ_dfsquare_matrix = zeros(nparams+2, nparams+2);
   end
   
   % Prediction Step
   theta_predict(:, j) = theta_update;
   Wvar_predict = Wvar_update + covar_matrix;
   
   % Compute mean
   mu(j) = mean_rate(k, H, nparams, theta_predict(:, j));
   if (mu(j) < 0), mu(j) = 1; end
   
   % Compute standard deviation
   std_vec(j) = ((mu(j).^3).*(theta_predict(nparams+2,j).^-1)).^0.5;
   
   % Compute f()
   f_vec(j) = f(j, theta_predict(end, j), delta, u(k), mu(j));
   if (f_vec(j) <= 1e-18)
       continue
   end
   
   % Compute sym_f() and integ_f()
   integ_f_temp = integ_f(j, delta, ...
       u(k), mu(j), theta_predict(end, j), integ_f_value);
   integ_f_value = integ_f_value + integ_f_temp; 
   if(integ_f_value <= 1e-28)
       continue
   end
   
   % Compute cif
   lambda_vec(j) = cif(f_vec(j), integ_f_value);
   
   % Compute del
   df_vec = df(j, k, nparams, delta, u(k), H, mu(j), f_vec(j), theta_predict(end, j));
   
   integ_df_vector = integ_df_vector + integ_df(j, k, delta, nparams,...
    u(k), H, mu(j), f_vec(j), theta_predict(end, j));
   
   d_lambda_vec = d_lambda(j, nparams, integ_f_value,...
       df_vec, integ_df_vector, f_vec);
   
   d_loglambda_vec = d_loglambda(lambda_vec(j), d_lambda_vec);
   
   % Compute del square
   dsquaref_mat = dsquaref(j, delta, k, nparams, u(k), H,...
    mu(j), f_vec(j), df_vec, theta_predict(end, j));

   integ_dfsquare_matrix = integ_dfsquare_matrix + integ_dfsquare(j,...
       delta, k, nparams, u(k), H, mu(j), f_vec(j), df_vec, theta_predict(end, j));

   dsquare_lambda_mat = dsquare_lambda(nparams, f_vec(j), integ_f_value, ...
    d_lambda_vec, df_vec, integ_df_vector, dsquaref_mat, integ_dfsquare_matrix);

   dsquare_loglambda_mat = dsquare_loglambda(nparams, lambda_vec(j), d_lambda_vec,...
    dsquare_lambda_mat);

   % Update Step
   theta_update = theta_predict(:, j)...
       + ((Wvar_predict * d_loglambda_vec') .* (n(j) - (lambda_vec(j)*delta)));
    
   Wvar_update = inv(inv(Wvar_predict) - dsquare_loglambda_mat.*(n(j)...
        - (lambda_vec(j)*delta)) - (d_loglambda_vec * (d_lambda_vec.*delta)')); 
end

%% R-R plots
% Mean R-R plot
t = (1:length(mu))*delta;
figure
hold on
plot(u, H, 'r*')
%line(u, H)
plot(t, mu)
xlabel('time in seconds')
ylabel('R-R intervals')
title('R-R intervals in Interruptions condition')
hold off


% Shifted r-r graph 
final_mean = mean(mu);
disp(final_mean);
disp(mean(H));
X = randn*0.0001;

temp1 = final_mean - mean(H);
temp2 = temp1; 
% shifted_mean = final_mean - temp2;
shifted_muvec = mu - temp2;
shifted_muvec(shifted_muvec < 0) = 0.1;
disp(mean(shifted_muvec));

figure
plot(t, shifted_muvec)
xlabel('time in seconds')
ylabel('mean R-R intervals')
title('mean of R-R intervals in Interruptions condition')
ax = gca;
ax.FontSize = 18;
%ax.YLim = [0 2];


% std plot
figure
plot(t, std_vec)
xlabel('time in seconds')
ylabel('standard deviation of R-R')
title('standard deviation in Interruptions condition')
ax = gca;
ax.FontSize = 16;
%ax.YLim = [0 0.1125];


%% Heart Rate 
mu_temp = shifted_muvec/60;
theta_temp = theta_predict(nparams+2, :)./60;
mu_temp(mu_temp < 0) = 0;
theta_temp(theta_temp < 0) = 0;

mu_hrv = mu_temp.^(-1) + theta_temp.^(-1);
mu_hrv(mu_hrv == Inf) = 0;

std_temp1 = 2.*mu_temp + theta_temp;
std_temp2 = mu_temp.*(theta_temp.^2);
std_hrv = (std_temp1./std_temp2).^0.5;


disp(mean(mu_hrv));

figure
plot(t, mu_hrv)
xlabel('time in seconds')
ylabel('mean HR in bpm')
title('Heart rate in Interruptions condition')
ax = gca;
ax.FontSize = 16;
ax.YLim = [0 150];

figure
plot(t, std_hrv)
xlabel('time in seconds')
ylabel('standard deviation of HR in bpm')
title('standard deviation of HR in Interruptions condition')
ax = gca;
ax.FontSize = 16;
ax.YLim = [0 3];
