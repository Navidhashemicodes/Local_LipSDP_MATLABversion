function  [alpha_param,beta_param]  =   piecewise_linear_slopes_tanh (l,u)
    
    diffil=tanh_diff(l);
    diffiu=tanh_diff(u);
    diffi0= ((l.*u)<=0).*tanh_diff(zeros(size(l)));
    alpha_param  = min(diffil,diffiu);
    beta_param   = max(diffi0, max(diffil,diffiu));
    
    function   y_prim = tanh_diff(x)
        epsi=0.0001;
        y_prim= (tanh(x+epsi)-tanh(x-epsi)) /  (2*epsi)  ;
    end
end