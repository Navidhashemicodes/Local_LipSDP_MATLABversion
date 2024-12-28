function [bound,time,status] =  lipsdp_local_lip(net,alpha_param,beta_param,options,mode)

language = options.language;
solver = options.solver;
verbose = options.verbose;

% if(strcmp(language,'yalmip'))
%     warning('YALMIP not yet available for lipsdp_local_lip. CVX will be used instead');
% end


bb = cat(1,net.biases{1:end-1});

% input dimension
dim_in = size(net.weights{1},2);

% output dimension
%dim_out = size(net.weights{end},1);

dim_last_hidden = size(net.weights{end},2);

% total number of neurons
num_neurons = length(bb);


if(strcmp(language,'cvx'))
    if(verbose)
        cvx_begin sdp
    else
        cvx_begin sdp quiet
    end
    variable lambda(num_neurons,1)
    variable rho_sq nonnegative
elseif(strcmp(language,'yalmip'))
    lambda = sdpvar(num_neurons,1);
    rho_sq = sdpvar(1);
end


Q = [-2*diag(alpha_param.*beta_param)*diag(lambda) diag((alpha_param+beta_param))*diag(lambda);...
    diag((alpha_param+beta_param))*diag(lambda) -2*diag(lambda)];

A = sparse(blkdiag(net.weights{1:end-1}));
A = sparse([A zeros(size(A,1),dim_last_hidden)]);
B = sparse([zeros(num_neurons,dim_in) eye(num_neurons)]);

Mmid = [A;B]'*Q*[A;B];

%CM_out = sparse([eye(dim_in+num_neurons) -eye(dim_in+num_neurons) zeros(dim_in+num_neurons,1)]);
El = [zeros(dim_last_hidden,dim_in+num_neurons-dim_last_hidden) eye(dim_last_hidden)];
E0 = [eye(dim_in) zeros(dim_in,num_neurons)];

Mout = (rho_sq*(E0'*E0)-((net.weights{end}*El)'*net.weights{end}*El));


if(strcmp(mode,'lower'))
    
    if(strcmp(language,'cvx'))
        
        maximize(rho_sq)
        subject to
        Mmid+Mout<=0;
        
        if(strcmp(net.activation,'relu'))
            Ipn = intersect(find(alpha_param==0),find(beta_param==1));
            lambda(Ipn)>=0;
        else
            lambda>=0;
        end
        
        cvx_end
        
        bound = sqrt(rho_sq);
        time= cvx_cputime;
        status = cvx_status;
        
    elseif(strcmp(language,'yalmip'))
        
        options = sdpsettings('solver',solver,'verbose',verbose);
        
        if(strcmp(net.activation,'relu'))
            Ipn = intersect(find(alpha_param==0),find(beta_param==1));
            constraints = [lambda(Ipn)>=0];
        else
            constraints = [lambda>=0];
        end
        
        constraints = [constraints, Mmid+Mout<=0];
        out = optimize(constraints,-rho_sq,options);
        
        bound = sqrt(value(rho_sq));
        time= out.solvertime;
        status = out.info;
    end
    
    
elseif(strcmp(mode,'upper'))
    
    if(strcmp(language,'cvx'))
        
        minimize(rho_sq)
        subject to
        Mmid-Mout<=0;
        
        if(strcmp(net.activation,'relu'))
            Ipn = intersect(find(alpha_param==0),find(beta_param==1));
            lambda(Ipn)>=0;
        else
            lambda>=0;
        end
        
        cvx_end
        
        bound = sqrt(rho_sq);
        time= cvx_cputime;
        status = cvx_status;
        
    elseif(strcmp(language,'yalmip'))
        
        options = sdpsettings('solver',solver,'verbose',verbose);
        
        if(strcmp(net.activation,'relu'))
            Ipn = intersect(find(alpha_param==0),find(beta_param==1));
            constraints = [lambda(Ipn)>=0];
        else
            constraints = [lambda>=0];
        end
        
        constraints = [constraints, Mmid-Mout<=0];
        out = optimize(constraints,rho_sq,options);
        
        bound = sqrt(value(rho_sq));
        time= out.solvertime;
        status = out.info;
    end
    
end



message = ['method: local lipsdp', '| solver: ', solver, '| bound: ', num2str(bound), '| time: ', num2str(time), '| status: ', status];
disp(message);



end