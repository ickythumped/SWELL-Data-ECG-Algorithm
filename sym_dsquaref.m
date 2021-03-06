function [sym_dsquaref_val] = sym_dsquaref(a, b, t, k, nparams, u_val, H,...
    mu_val, f_val, df_vec, theta_p_plus1)
%handler function

%% reg1
expr1 = t - u_val;
expr2 = (expr1 - mu_val)./(mu_val^3);
expr3 = (2*mu_val - 3.*expr1)./(mu_val^4);
df_dmu = (f_val.*theta_p_plus1.*(expr1 - mu_val))./(mu_val.^3);
dsquaref_dmusquare = theta_p_plus1.*(df_dmu.*expr2 + f_val.*expr3);

%% reg3
expr4 = ((expr1 - u_val).^2)./(2.*(mu_val^2).*expr1);
expr5 = df_vec(end) - f_val./theta_p_plus1;

%% reg2
expr6 = (f_val*(expr1 - mu_val))./mu_val^3;

%%
if (a == 1 && b == 1)
    sym_dsquaref_val = dsquaref_dmusquare;
    
elseif (a == 1 && b <= nparams+1 && b > 1)
    sym_dsquaref_val = dsquaref_dmusquare.*H(k-b+2);
    
elseif (a > 1 && a <= nparams+1  && b == 1)
    sym_dsquaref_val = dsquaref_dmusquare.*H(k-a+2);
    
elseif (a <= nparams+1 && b <= nparams+1 && a > 1 && b > 1) %reg1
    sym_dsquaref_val = dsquaref_dmusquare.*H(k-a+2).*H(k-b+2);
    
elseif (a == 1 && b == nparams+2)
    sym_dsquaref_val = ((1/(2*theta_p_plus1))*df_vec(a))...
        - (df_vec(a)*expr4) + expr6;
    
elseif (a == nparams+2 && b == 1)
    sym_dsquaref_val = ((1/(2*theta_p_plus1))*df_vec(a))...
        - (df_vec(a)*expr4) + expr6;

elseif (a == nparams+2 && b == nparams+2) %reg3
    sym_dsquaref_val = (df_vec(end).*(-expr4))...
            + ((1/(2.*theta_p_plus1)).*expr5);

elseif (a <= nparams+1 && b == nparams+2)
    sym_dsquaref_val = ((1/(2*theta_p_plus1))*df_vec(a))...
        - (df_vec(a)*expr4) + (expr6*H(k-a+2));

elseif (a == nparams+2 && b <= nparams+1)
    sym_dsquaref_val = ((1/(2*theta_p_plus1))*df_vec(b))...
        - (df_vec(b)*expr4) + (expr6*H(k-b+2));
end

end
