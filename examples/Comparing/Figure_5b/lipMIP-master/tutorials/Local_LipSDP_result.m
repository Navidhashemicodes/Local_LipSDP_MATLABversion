clear all
clc
close all

dir = pwd;
for i=1:5
    dir = fileparts(dir);
end
addpath(genpath([dir , '\matlab_engine']))
addpath(genpath([dir , '\yalmip']))


load('2by2comparison.mat')

Net = net;
net = nnsequential(Net.dims, Net.activation);
net.weights = Net.weights;
net.biases = Net.biases;

options.language='yalmip';
options.solver='mosek';
options.verbose=1;
mode='upper';
repeated=0;

epsilons= linspace(0.001, 1, 10);



dims = net.dims;
dimin= dims(1);
X=zeros(dimin,1);
dimout= dims(end);
len = length(epsilons);

bound = zeros(1, len);
time = zeros(1, len);
status = cell(1, len);

parfor i=1:len

    epsilon=epsilons(i)*ones(size(X));
    [alpha_param, beta_param] = net.slope_bounds(X , epsilon);
    [bound(i),time(i),status{i}] =   lipsdp_local_lip(net,alpha_param,beta_param,options,mode);

end

clc

disp(['The Lipstchitz certificates, for each epsilon in linspace(0.001,1,10), from Local_LipSDP are: ', num2str(bound), ' respectively.'] )


