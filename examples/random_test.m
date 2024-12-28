clc;
clear all;

rng('default');

dims = [2,10,10,2];
activation = 'tanh';

dir = pwd;
dir = fileparts(dir);
addpath(genpath([dir , '\matlab_engine']))
addpath(genpath([dir , '\yalmip']))

net = nnsequential(dims,activation);

num_neurons = sum(net.dims(2:end-1));

center = zeros(net.dims(1),1);
epsilon = 100*ones(net.dims(1),1);

[alpha_param , beta_param] = net.slope_bounds(center, epsilon);

options.language = 'yalmip';
options.verbose = 1;
options.solver = 'mosek';

mode = 'upper';


[bound,~,~] =  lipsdp_local_lip(net,alpha_param,beta_param,options,mode);
disp(bound);

%%%%%%%%%%%%%%%%%%%%%%

center = zeros(net.dims(1),1);
epsilon = 0.1*ones(net.dims(1),1);

[alpha_param , beta_param] = net.slope_bounds(center, epsilon);


[bound,~,~] =  lipsdp_local_lip(net,alpha_param,beta_param,options,mode);
disp(bound);

