Defaults

options = odeset('RelTol',1e-9,'nonnegative',1);

resultsY = [];

%%
% test = [1,2,3,5,8,10]; %STAT3npd
% testthreshold = 6;
% testindex = find(test>testthreshold, 1, 'first');
% time = [1,2,3,4,5,6]

threshold = 400;
Kd2=1090;
KDarray = 10.^[-3:0.25:1.5]*Kd2;

for i = 1:length(KDarray)
    
ii = [1007, i];
[T,Y] = ode15s( @(t,y) ODEs(t,y,ii) ,T,Y0,options);
 

index = find(Y(:,8)>threshold, 1, 'first');
timevalue = T(index)/60;

if isempty(timevalue)
    disp('Error: STAT3npd does not reach threshold within modelled time')
else
    resultsY(i) = timevalue;
end

end

semilogx(KDarray, resultsY, '.', 'MarkerSize',15);
xlabel('K_D_2 [nM]');
ylabel('Time [minutes]');
title(['Time for STAT3npd to reach ', num2str(threshold), ' nM dependent on K_D_2']);
hold on
standard = find(KDarray == 1090);
plot(KDarray(standard), resultsY(standard), '.r', 'MarkerSize', 15);

