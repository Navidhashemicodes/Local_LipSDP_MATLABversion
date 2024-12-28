function [alpha_param,beta_param] =   piecewise_linear_slopes_ReLU(l,u)
    
    alpha_param  =  (l>0) .* ones(size(l));
    beta_param   =  (u>=0).* ones(size(u));
    
end