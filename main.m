clear all; %#ok<CLALL>
close all;
clc;
load('r_peaks.mat');

%% Declarations
% u -> Location of R-peaks in seconds
% K -> Number of R-peaks
% T -> Length of data [0, T]
% H -> History vector
[u, K, T, H] = init(r_peaks);

%% j division
% J -> Divide [0,T] into J equal parts
% delta -> Length of each part in seconds
% n_j -> binary value for peak in each delta
[J, delta, n] = del_parts(T, u);

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
sigma_square = zeros(1, J); %variance of each interval

%% Initializations (Note: work to be done)
theta_update = [0.14937830871589303; 1.1564059506167039; -0.39835826229059518; ...
    0.83088601289788633; -1.023213140764317; 0.41452982839565033; -0.26784526179216739; ...
    0.6149154253134137; -0.60566867990517192; 0.11520209772910678; 1507.727429652082];
theta_predict(:, start_iter) = theta_update;
mu(start_iter) = mean_rate(k, H, nparams, theta_update);

% sigma_square(start_iter) = 0.00075;
% theta_predict(nparams+2, start_iter) = mu(start_iter)^3./sigma_square(start_iter);
% theta_update(nparams+2) = mu(start_iter)^3./sigma_square(start_iter);

covar_matrix = diag(theta_update);
Wvar_predict = ones(nparams+2, nparams+2); %Model parameter vector
Wvar_update = covar_matrix; %Model parameter vector

% Integral value initializations
integ_f_value = 0;
integ_df_vector = zeros(1, nparams+2);
integ_dfsquare_matrix = zeros(nparams+2, nparams+2);

% Paramter initializations
lambda_vec = zeros(1, J); % lambda vector
f_vec = zeros(1, J); %f() vector
%% Algorithm without adaptive filter

for j = start_iter:J
    
   % Check for previous spike
   if(k == length(u))
       break
   end
   
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
   
   % Compute f()
   f_vec(j) = f(j, theta_predict(end, j), delta, u(k), mu(j));
   if (f_vec(j) <= 1e-18)
       continue
   end
%    if (f_vec(j) == Inf)
%        continue
%    end
   
   % Compute sym_f() and integ_f()
   integ_f_value = integ_f_value + integ_f(j, delta, ...
       u(k), mu(j), theta_predict(end, j), integ_f_value); %check after completing del and update eqs
   if(integ_f_value <= 1e-18)
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

%% R-R plot
t = (1:length(mu))*delta;
figure
hold on
plot(u, H, 'r.')
%line(u, H)
plot(t, mu)
xlabel('time in seconds')
ylabel('R-R intervals')
title('R-R Plot')
hold off

%% Mu_value
final_mean = mean(mu);
disp(final_mean);
disp(mean(H));

