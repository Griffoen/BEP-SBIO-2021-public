%script to create all the plots needed for a report/validation
PlotOptions = char('all','Tests', 'STAT/RJL', 'LigandvsReceptor', 'SEAP variations', 'BellCurve');

nOpt = size(PlotOptions,1);
    
disp(['scenarios available:'])    
disp(['------------------------------------'])
    for iplot=1:nOpt
        disp([ num2str( iplot,'%10.0f' ) '. ' PlotOptions(iplot,:) ])
    end
disp(['------------------------------------'])
disp([' '])
    
myChoice=input(['choose plots (1-' num2str( nOpt,'%10.0f' ) '): ']);



if myChoice == 2 || myChoice == 1
%First run TEST_ODE to see if the ODEs are correct and function as expected
TEST_ODE;

input('Press any key to continue');
close all; clearvars -except myChoice;
end

if myChoice == 3 || myChoice == 1
%plot STAT3 and receptor over time
OdeSolver

input('Press any key to continue');
close all; clearvars -except myChoice;
end

if myChoice == 4 || myChoice == 1
Defaults
Plot_LigandvsReceptor(T,Y0,labels,[]);
input('Press any key to continue');
close all; clearvars -except myChoice;
end

if myChoice == 5 || myChoice == 1
Defaults;
item = 1; % 1: vary ligand concentrations, 2: vary receptor concentrations, 3: vary KB7, 4: vary Kd
Plot_SEAP_Variations(T,Y0,item)
input('Press any key to continue');
close all; clearvars -except myChoice;
end

if myChoice == 6 || myChoice == 1
Defaults;
time = 96;
Plot_BellCurve_LvsR(T,Y0,[time]);
input('Press any key to continue');
close all; clearvars -except myChoice;
end

if myChoice == 7 || myChoice == 1
Vary_Dimerization_KD;
input('Press any key to continue');
close all; clearvars -except myChoice;
end


