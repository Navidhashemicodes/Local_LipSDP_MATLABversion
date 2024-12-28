clear
clc

dir = pwd;
dir = fileparts(dir);
dir = fileparts(dir);
addpath(genpath([dir , '\matlab_engine']))
addpath(genpath([dir , '\yalmip']))


load('Input_data.mat')
load('Trained_model.mat')
options.language='yalmip';
options.solver='mosek';
options.verbose=1;
mode='upper';
% Input = S;
[ alpha_param , beta_param ] = net.slope_bounds( Input, 100*ones(size(Input)) );
L_global = lipsdp_local_lip(net,alpha_param , beta_param, options, mode);

[L, epsilon, rho ]  =  Robustness_bisection(net, Input, L_global);
clc



epsilon_global=epsilon.*L/L_global;


save('Results.mat', 'L', 'epsilon', 'L_global', "epsilon_global", "rho", "net", "Input")