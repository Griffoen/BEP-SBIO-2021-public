function x = PlotBellCurve_LvsR(T,Y,Options)
Y0 = Y;
t_end = Options(1);

options = odeset('RelTol',1e-9,'nonnegative',1);

% set loop ratios:
ratios = 10.^[-6:0.5:2]; %[10^-6, 10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 1, 10, 100]; %logarithmic x scale
%ratios = 10.^[-3.5:0.25:-1];
NR_ratios = size(ratios,2);
NR_odes = 22;

% custom time:
t_step        = 1;  % [sec]
t_end         = 3600*t_end; % [sec]
T = 0:t_step:t_end; % [sec]

timepoints = size(T,2);
% configure output store:
Yout = zeros(size(ratios,2) ,timepoints, NR_odes);


%% Target
target = 'SEAPex';
if strcmp(target,'RJ2L')
    targetID = 21;
elseif strcmp(target,'RJ2Lp')
    targetID = 22;
end


%% Mainloop
% note:
% the upper for loop is parfor capable, if we need it. just know that
% parfor and plotting does not go well together. however, feel free to
% analyse the results from parfor, or save them for later plotting 

% store the amout of receptor for ligand variation
RJ = Y0(19);

parfor i = 1:size(ratios,2) % vary one constant
        ii=[]; % Keep ii empty to disable constant variation
        
        if strcmp(target,'RJ2L')
            ii=[1003]; % special flag to disable KP1. This prevents conversion RJ2L to RJ2Lp
        end
        
        Y0i=Y0;
        Y0i(18) = RJ*ratios(i);

        % solve the model
        [~, Yout(i,:,:) ]=ode15s( @(t,y) ODEs(t,y,ii) ,T,Y0i,options);
end

Last = size(Yout,2)-1;


%% Next plots
c_min = 0.3;
% change color settings
newcolors = zeros(NR_ratios,3);
newcolors(:,1) = 0;
newcolors(:,2) = c_min:(1-c_min)/(NR_ratios-1):1;
%newcolors(:,3) = 1:-1/(NR_ratios-1):0;
set(gca, 'ColorOrder', newcolors, 'NextPlot', 'replacechildren');

Legend_data{NR_ratios} = [];
for j = 1:NR_ratios
    Legend_data{j} = ['RJ:L ratio; [1:', num2str(ratios(j)),']'];
end

% see if we have a steady state
figure(1)
set(gca, 'ColorOrder', newcolors, 'NextPlot', 'replacechildren');
% all RJ2L
% all_RJ2L = sum(Yout(:,:,[targetID]),3);
all_SEAP = sum(Yout(:,:,[17]),3);
% plot(T/3600, all_RJ2L,'-', 'LineWidth',1.5);
plot(T/3600, all_SEAP,'-', 'LineWidth',1.5);
legend(Legend_data)
xlabel('Time [hour]')
ylabel('Concentration [AU]')
% title('RJ2Lp')
title('SEAPex')



figure(2)
set(gca, 'ColorOrder', newcolors, 'NextPlot', 'replacechildren');
% Bell curve of RJ2Lp

Last = size(Yout,2)-1;
BellTime = Last; %1*60*60 ; %[Sec] time of concentrations in Bell curve 
% RJ2Lp_bell = all_RJ2L(:,BellTime);
SEAP_bell = all_SEAP(:,BellTime);
% plot(ratios, RJ2Lp_bell,'.', 'MarkerSize',15);
%plot(ratios, SEAP_bell,'.', 'MarkerSize',15);
semilogx(ratios, SEAP_bell,'.', 'MarkerSize',15);  %logarithmic x scale
xlabel('RJ:L ratio');
ylabel('Concentration [AU]');
%xlim([10^-4, 1]);
%title("Bell curve SEAPex at " + BellTime/3600 + " hours");
Title = "\fontsize{14}\color{black}\bfBell curve SEAPex at " + BellTime/3600 + " hours";
subtitle = '\fontsize{10}\color{gray}\rmK_D_1 = 0.1';
title({Title; subtitle});
end
