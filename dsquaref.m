function [dsquaref_mat] = dsquaref(j, delta, k, nparams, u_val, H,...
    mu_val, f_val, df_vec, theta_p_plus1)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% intializations
dsquaref_reg1_mat = zeros(nparams+1, nparams+1);
dsquaref_reg2_vec = zeros(1, nparams+1);

%% reg1
expr1 = j*delta - u_val;
expr2 = (expr1 - mu_val)./(mu_val^3);
expr3 = (2*mu_val - 3.*expr1)./(mu_val^4);
df_dmu = (f_val.*theta_p_plus1.*(expr1 - mu_val))./(mu_val.^3);
dsquaref_dmusquare = theta_p_plus1.*(df_dmu.*expr2 + f_val.*expr3);

%% reg3
expr4 = ((expr1 - u_val).^2)./(2.*(mu_val^2).*expr1);
expr5 = df_vec(end) - f_val./theta_p_plus1;

%% reg2
expr6 = (f_val*(expr1 - mu_val))./mu_val^3;

%% Computing all regions
for a = 1:nparams+1
    %reg2
    dsquaref_reg2_vec(a) = ((1/(2*theta_p_plus1))*df_vec(a)) -...
        (df_vec(a)*expr4) + (expr6*H(k-a+2));
    
    for b = 1:nparams+1
        %reg1
        dsquaref_reg1_mat(a,b) = dsquaref_dmusquare*H(k-a+2)*H(k-b+2);
    end
end

%reg3
dsquaref_reg3_val = (df_vec(end)*(-expr4)) + ((1/(2*theta_p_plus1))*expr5);

%% Creating Hessian matirx
dsquaref_temp1 = [dsquaref_reg1_mat dsquaref_reg2_vec'];
dsquaref_temp2 = [dsquaref_reg2_vec dsquaref_reg3_val];
dsquaref_mat = [dsquaref_temp1; dsquaref_temp2];

end

