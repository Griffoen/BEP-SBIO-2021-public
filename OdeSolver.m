Defaults

options = odeset('RelTol',1e-9,'nonnegative',1);

%% Mainloop
% note:
% the upper for loop is parfor capable, if we need it. just know that
% parfor and plotting does not go well together. However, feel free to
% analyse the results from parfor, or save them for later plotting 

for i = 1:1 % vary one constant
    for j = 1:1 % vary another
        
        %ii = [i,j];
        ii=[]; % Keep ii empty to disable constant variation
        
        % solve the model
        [T,Y]=ode15s( @(t,y) ODEs(t,y,ii) ,T,Y0,options);
        
        % any analasis here:
        % or define a function later to analyze the results with.
        % out = DoAnalasis(Y,); % or something like it.
        colors = [0.4 0.4 0.4; 0.8 0.5 0.2; 0.9 0.75 0; 0.5 0.2 0.5; 0.25 0.80 0.54; 1 0 0; 0.2 0.2 0.2]; %STAT visualisation
        lines = {'-', '--', ':'};
        Title = "\fontsize{14}\color{black}\bfSEAP secretion";
        PlotYx(T,Y, [14:17],labels, colors, lines, Title);
        xlim([0, 8]);
        ylim([0, 500]);
        
        colors = [0.8 0.5 0.2; 0.5 0.2 0.5; 0.2 0.4 0.9; 1 1 0; 0 1 1; 0.5 0.5 0.5; 0 1 0]; %receptor and ligand
        Title = "\fontsize{14}\color{black}\bfSTAT circle";
        PlotYx(T,Y, [1,7:9],labels, colors, lines, Title); %,3,12,11,6,2,4
        xlim([0, 5]);
        ylim([0, 500]);
        
        colors = [0.8 0.5 0.2; 0.5 0.2 0.5; 0.2 0.4 0.9; 1 1 0; 0 1 1; 0.5 0.5 0.5; 0 1 0]; %receptor and ligand
        Title = "\fontsize{14}\color{black}\bfReceptor species";
        PlotYx(T,Y, [19:22,2],labels, colors, lines, Title);
        ylim([0 12]);
    end
end

%% Plot
function PlotYx(t,Y,x,labels,colors,lines,Title)
figure; hold on;
 
set(gca, 'ColorOrder', colors, 'LineStyleOrder', lines');
reset(gca);
plot(t/3600, Y(:,x), 'Linewidth', 1.5);
%ylim([0,200]);
xlabel('Time [hour]')
ylabel('Concentration [nM]')
legend(labels(x))

% Create a title and a subtitle with the initial conditions
subtitle = '\fontsize{10}\color{gray}\rmInitial [nM]: ';
Y0 = Y(1,:);
for n = 1:length(Y0)
    if Y0(n) ~= 0
        subtitle = strcat(subtitle, labels(n), '=', num2str(Y0(n)), ',');
    end
end
subtitle = string(subtitle);
title({Title; subtitle});
end