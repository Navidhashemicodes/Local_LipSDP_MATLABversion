clear all
clc
close all

dir = pwd;
dir = fileparts(dir);
dir = fileparts(dir);
addpath(genpath([dir , '\matlab_engine']))
addpath(genpath([dir , '\yalmip']))


% Prepare tiled layout
t = tiledlayout(2, 4, 'TileSpacing', 'Compact', 'Padding', 'Compact');

cases = {
    [50,20,10],[50,20,20,10], [50,20,20,20,10],[50,20,20,20,20,10], ...
    [100,50,10], [100,50,50,10], [100,50,50,50,10], [100,50,50,50,50,10] };
case_results = cell(1,size(cases, 1));
case_info = cell(1,size(cases, 1));


options.language='yalmip';
options.solver='mosek';
options.verbose=1;
mode='upper';
repeated=0;

NNN=50;
epsilons= [linspace(0.01, 0.51, 26), 0.7, 0.9 , 1 ];

% Iterate over each case
for m = 1:size(cases, 2)


    dims = cases{m};
    dimin= dims(1);
    X=zeros(dimin,1);
    dimout= dims(end);
    
    if m>4
        NNN = 20;
    end

    bound = zeros(NNN , length(epsilons));
    time = zeros(NNN , length(epsilons));
    status = cell(NNN , length(epsilons));
    bound_global = zeros(NNN,1);

    
    for ii=1:NNN
        net = nnsequential(dims , 'tanh');

        parfor jj=1:length(epsilons)

            epsilon=epsilons(jj)*ones(size(X));
            [alpha_param, beta_param] = net.slope_bounds(X , epsilon);
            [bound(ii,jj),time(ii,jj),status{ii,jj}] =   lipsdp_local_lip(net,alpha_param,beta_param,options,mode);

        end
        epsilon=100*ones(size(X));
        [alpha_param, beta_param] = net.slope_bounds(X , epsilon);
        bound_global(ii,1) =   lipsdp_local_lip(net,alpha_param,beta_param,options,mode);
    end

    case_results{m} = bound ./ bound_global; 
    case_info{m} = { time, status };
    % Plot data
    nexttile;
    plot_areaerrorbar(case_results{m}, epsilons);
    
    % Format axes
    ax = gca;
    ax.FontSize = 16;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    xlim([0, 1]);
    ylim([0.6, 1.1]);
    ax.XTick = [0 0.2 0.4 0.6 0.8 1];
    ax.YTick = [0.6 0.7 0.8 0.9 1 1.1];
    ax.LineWidth = 4;
    ylabel('$L_{loc}/L_{global}$', 'Interpreter', 'latex', 'FontSize', 25);
    
    % Add labels and annotations
    text(0.5, 0.55, '$\epsilon$', 'FontSize', 40, 'Interpreter', 'latex');
    text(-0.2, 1.05, sprintf('(%c)', 'a' + m - 1), 'FontSize', 35, 'Interpreter', 'latex');
    box on;
end

print('-painters','-depsc','Casestudy1');

save("Experiment_1.mat", "cases", "case_results", "case_info");