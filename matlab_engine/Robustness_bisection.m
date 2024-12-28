%%%%% epsilon finder
function  [L,epsilon,rho]  =  Robustness_bisection(net, input_data, L_g)

dims_in = net.dims(1);
Output_data=net.eval(input_data);

options.language='yalmip';
options.solver='mosek';
options.verbose=1;
mode='upper';

i=find(Output_data==max(Output_data));
i=i(1);
disp(['The introduced data point is classified to be   ', num2str(i-1)]);
leng=length(Output_data);
I=eye(leng);
s=0;
for j=1:leng
    if i~=j
        s=s+1;
        Rho(s)=abs(   (I(:,i)-I(:,j)).'   * Output_data  /sqrt(2)   );
    end
end
rho=min(Rho);


epsi_lower=0.001;
epsi_upper=2;   %%% we assume L_min can not be smaller than 0.01*L_global
while epsi_upper-epsi_lower>=0.0001
    epsi=0.5*(epsi_lower+epsi_upper);
    epsi_vector=epsi*ones(size(input_data));   %%% this command surrounds the sphere with a cube.
    [alpha_param , beta_param] = net.slope_bounds(input_data, epsi_vector);
    bound =   lipsdp_local_lip(net, alpha_param, beta_param, options, mode);
    if bound*epsi*sqrt(dims_in)>=rho
        epsi_upper=epsi;
    else
        epsi_lower=epsi;
    end
end

L=bound;
epsilon=epsi;


end