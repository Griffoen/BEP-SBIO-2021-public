function x = Plot_SEAP_Variations(T,Y,Options)
Y0 = Y;

options = odeset('RelTol',1e-9,'nonnegative',1);

%% select the variated parameter

item = Options; % 1: vary ligand concentrations, 2: vary receptor concentrations, 3: vary KB7, 4: vary Kd1

% set loop variations
switch item
    case 1
        % variations in ligand concentration
        range = [1e-3 1e-2, 1e-1, 1, 1e1, 1e2 1e3]*(Y0(18)); % variations ligand concentration
        
        NR_parameter = 18;
        
        NR_variations = size(range,2);
        Legend_data{NR_variations} = [];
        vary = 'Ligand concentration [nM]';
        for j = 1:NR_variations
            Legend_data{j} = ['Ligand concentration: ', num2str(range(j)),'[nM]'];
        end
        
    case 2
        % variations in receptor concentration
        range = [1e-3 1e-2, 1e-1, 1, 1e1, 1e2 1e3]*Y0(19); % variations receptor concentration
        
        NR_parameter = 19;
        
        NR_variations = size(range,2);
        Legend_data{NR_variations} = [];
        vary = 'Receptor concentration [nM]';
        for j = 1:NR_variations
            Legend_data{j} = ['Receptor concentration: ', num2str(range(j)),'[nM]'];
        end
        
    case 3
        % variations in receptor ligand binding constant (KB7)
        range = [1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3]*0.43; %KB7 is 0.43
        
        NR_parameter = 100; % in order to not adjust the Yout values
        
        NR_variations = size(range,2);
        Legend_data{NR_variations} = [];
        vary = 'KB7 [1/(nM*s)]';
        for j = 1:NR_variations
            Legend_data{j} = ['KB7: ', num2str(range(j)), '[1/(nM*s)]'];
        end
        
    case 4
        % variations in Kd1
        range = [1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3]* 0.1; %Kd1 is 0.1
        
        NR_parameter = 101; % in order to not adjust the Yout values
        
        NR_variations = size(range,2);
        Legend_data{NR_variations} = [];
        vary = 'Kd1 [nM]'; 
        for j = 1:NR_variations
            Legend_data{j} = ['Kd: ', num2str(range(j)),'[nM]'];
        end
end
        

NR_odes = 22;
timepoints = size(T,2);
% configure output store:
Yout = zeros(size(range,2) ,timepoints, NR_odes); 

%% Mainloop
% note:
% the upper for loop is parfor capable, if we need it. just know that
% parfor and plotting does not go well together. however, feel free to
% analyse the results from parfor, or save them for later plotting 

% vary the ligand concentration
parfor i = 1:NR_variations % vary one constant
        
        if NR_parameter == 100;
            ii = [1004,i]; % execute constant variation for KB7
        elseif NR_parameter == 101;
            ii = [1005,i]; % execute constant variation for Kd
        else 
            ii=[]; % Keep ii empty to disable constant variation
        end
            
        
        %ii=[1001]; % special flag to disable KB2. This prevents STAT3p from binding to RJ2Lp
        
        Y0i=Y0;
        
        % set the beginconcentration
        if NR_parameter < 100;
            Y0i(NR_parameter) = range(i); % vary the concentration of the parameter
        end

        % solve the model
        [~, Yout(i,:,:) ]=ode15s( @(t,y) ODEs(t,y,ii) ,T,Y0i,options);
end

%% Plot

newcolors = zeros(NR_variations,3);
newcolors(:,1) = 0;
newcolors(:,2) = 0:1/(NR_variations-1):1;
%newcolors(:,3) = 1:-1/(NR_variations-1):0;
set(gca, 'ColorOrder', newcolors, 'NextPlot', 'replacechildren');

figure(1)
set(gca, 'ColorOrder', newcolors, 'NextPlot', 'replacechildren');
% all SEAPex
all_SEAP = sum(Yout(:,:,[17]),3);
plot(T/3600, all_SEAP,'-', 'LineWidth',1.5);
legend(Legend_data)
xlabel('Time [hour]')
ylabel('Concentration [nM]')
title('SEAPex')

%% Plot 2, at certain timepoint
figure(2)
MeasurementTime = 3600; %[sec]
SEAP_values = all_SEAP(:,MeasurementTime);
%Make sure only relevant values get plotted
SEAP_mask = SEAP_values >= 1e-6;
SEAP_values = SEAP_values(SEAP_mask);
range_values = range(SEAP_mask);
%Plot everything on logarithmic scale
semilogx(range_values, SEAP_values,'.', 'MarkerSize',15);
xlabel(vary);
ylabel('Concentration SEAPex [nM]');
title("SEAPex concentration at " + MeasurementTime/3600 + " hours dependent on " + vary);
end
