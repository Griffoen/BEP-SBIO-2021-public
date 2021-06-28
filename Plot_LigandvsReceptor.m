function x = Plot_LigandvsReceptor(T,Y,labels,Options)
Y0 = Y;

options = odeset('RelTol',1e-9,'nonnegative',1);

% set loop ratios:
ratios = [10^-4, 10^-3, 10^-2, 10^-1, 1]; %[0:0.5:5]/100;
NR_ratios = size(ratios,2);
NR_odes = 22;
timepoints = size(T,2);
% configure output store:
Yout = zeros(size(ratios,2) ,timepoints, NR_odes);

%% Mainloop
% note:
% the upper for loop is parfor capable, if we need it. just know that
% parfor and plotting does not go well together. however, feel free to
% analyse the results from parfor, or save them for later plotting 

% store the amout of receptor for ligand variation
RJ = Y0(19);

parfor i = 1:size(ratios,2) % vary one constant
        ii=[]; % Keep ii empty to disable constant variation
        
        %ii=[1001]; % special flag to disable KB2. This prevents STAT3p from binding to RJ2Lp
        
        Y0i=Y0;
        Y0i(18) = RJ*ratios(i);

        % solve the model
        [~, Yout(i,:,:) ]=ode15s( @(t,y) ODEs(t,y,ii) ,T,Y0i,options);
end


%% Plot the ligand vs receptor concentration

% first figure; plot receptor/ligand behavoir
X = [18:22,2,4];
Yindex = [0,.5,1]*(NR_ratios-1)+1;
figure(1)
axes = zeros(1,3);
for i = 1:3
    axes(i) = subplot(1,3,i);
    plot(T/3600, squeeze(Yout(Yindex(i),:,X)),'LineWidth',1.5);
    title(['Ratio [Receptor:Ligand] ; [1:', num2str(ratios(Yindex(i))),']'])
end
linkaxes(axes,'xy')
xlabel('Time [hour]')
ylabel('Concentration [nM]')
legend(labels(X))

%% Next plots
c_min =0.3;
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

figure(2)
set(gca, 'ColorOrder', newcolors, 'NextPlot', 'replacechildren');
% all RJ2Lp
all_RJ2Lp = sum(Yout(:,:,[22,2,4]),3);
plot(T/3600, all_RJ2Lp,'-', 'LineWidth',1.5);
legend(Legend_data)
xlabel('Time [hour]')
ylabel('Concentration [nM]')
title('Total RJ2Lp')
Ylim = ylim;


figure(3)
set(gca, 'ColorOrder', newcolors, 'NextPlot', 'replacechildren');
% RJ2Lp
RJ2Lp = sum(Yout(:,:,[22]),3);
plot(T/3600, RJ2Lp,'-', 'LineWidth',1.5);
legend(Legend_data)
xlabel('Time [hour]')
ylabel('Concentration [nM]')
title('RJ2Lp')

figure(4)
set(gca, 'ColorOrder', newcolors, 'NextPlot', 'replacechildren');
% all RJ2Lp_STAT3c
all_RJ2Lp_STAT3c = sum(Yout(:,:,[2]),3);
plot(T/3600, all_RJ2Lp_STAT3c,'-', 'LineWidth',1.5);
legend(Legend_data)
xlabel('Time [hour]')
ylabel('Concentration [nM]')
title('RJ2Lp_STAT3c')

figure(5)
set(gca, 'ColorOrder', newcolors, 'NextPlot', 'replacechildren');
% all STAT3cp
all_STAT3cp = sum(Yout(:,:,[3]),3);
plot(T/3600, all_STAT3cp,'-', 'LineWidth',1.5);
legend(Legend_data)
xlabel('Time [hour]')
ylabel('Concentration [nM]')
title('STAT3cp')

figure(6)
set(gca, 'ColorOrder', newcolors, 'NextPlot', 'replacechildren');
% all STAT in cytoplasm
all_STAT3c = sum(Yout(:,:,[1:4,6,7]),3);
plot(T/3600, all_STAT3c,'-', 'LineWidth',1.5);
legend(Legend_data)
xlabel('Time [hour]')
ylabel('Concentration [nM]')
title('Total STAT3c')

figure(7)
set(gca, 'ColorOrder', newcolors, 'NextPlot', 'replacechildren');
% all STAT3cpd
all_statcpd = sum(Yout(:,:,[7]),3);
plot(T/3600, all_statcpd,'-', 'LineWidth',1.5);
legend(Legend_data)
xlabel('Time [hour]')
ylabel('Concentration [nM]')
title('STAT3cpd')

figure(8)
set(gca, 'ColorOrder', newcolors, 'NextPlot', 'replacechildren');
% all STAT3npd
all_statnpd = sum(Yout(:,:,[8]),3);
plot(T/3600, all_statnpd,'-', 'LineWidth',1.5);
legend(Legend_data)
xlabel('Time [hour]')
ylabel('Concentration [nM]')
title('STAT3npd')

figure(9)
set(gca, 'ColorOrder', newcolors, 'NextPlot', 'replacechildren');
% all SEAPex
all_SEAP = sum(Yout(:,:,[17]),3);
plot(T/3600, all_SEAP,'-', 'LineWidth',1.5);
legend(Legend_data)
xlabel('Time [hour]')
ylabel('Concentration [nM]')
title('SEAPex')


figure(10);
set(gca, 'ColorOrder', newcolors, 'NextPlot', 'replacechildren');
% RJ2Lp_STAT3p
all_InHib = sum(Yout(:,:,[4]),3);
plot(T/3600, all_InHib,'-', 'LineWidth',1.5);
legend(Legend_data)
xlabel('Time [hour]')
ylabel('Concentration [nM]')
title('RJ2Lp\_STAT3p')
ylim(Ylim)

figure(11);
set(gca, 'ColorOrder', newcolors, 'NextPlot', 'replacechildren');
% RJ2Lp_STAT3p
all_InHib = sum(Yout(:,:,[1]),3);
plot(T/3600, all_InHib,'-', 'LineWidth',1.5);
legend(Legend_data)
xlabel('Time [hour]')
ylabel('Concentration [nM]')
title('STAT3c')

figure(12);
set(gca, 'ColorOrder', newcolors, 'NextPlot', 'replacechildren');
% mRNAc_SEAP
all_InHib = sum(Yout(:,:,[14]),3);
plot(T/3600, all_InHib,'-', 'LineWidth',1.5);
legend(Legend_data)
xlabel('Time [hour]')
ylabel('Concentration [nM]')
title('mRNAc_SEAP')

figure(13);
set(gca, 'ColorOrder', newcolors, 'NextPlot', 'replacechildren');
% nuclear STAT3
all_InHib = sum(Yout(:,:,[8,8,9,11,12]),3);
plot(T/3600, all_InHib,'-', 'LineWidth',1.5);
legend(Legend_data)
xlabel('Time [hour]')
ylabel('Concentration [nM]')
title('nuclear STAT3')



figure() % maker shure the last figure gets its title.
% reset the color sceme
reset(gca)
end