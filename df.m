function [df_vec] = df(j, k, nparams, delta, u_val, H, mu_val, f_val, theta_p_plus1)
%Computing (df/dtheta) for all thetas

%  dtheta = [dtheta(1) dtheta(2)... dtheta(p+2)]
df_vec = zeros(1, nparams+2);

expr1 = (f_val*theta_p_plus1)./(mu_val^3);
expr2 = (j*delta - u_val - mu_val);

for i = 1:nparams+1
    if(i == 1)
        df_vec(1) = expr1*expr2;
    else
        df_vec(i) = H(k-i+2)*expr1*expr2;
    end
end

expr3 = (-(expr2)^2)./(2*(mu_val^2)*(j*delta - u_val));
expr4 = f_val./(2*theta_p_plus1);
df_vec(nparams+2) = f_val*expr3 + expr4;
end

