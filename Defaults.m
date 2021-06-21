%% Model defaults
% Call this matlabfile at the strating of a ODE solver to set the defaults
close all;

%% Configure time
t_step        = 1;               % [sec]
t_end         = 3600*24;          % [sec] +/- 4 hours
T             = 0:t_step:t_end;  % [sec]

%% Define the labels so that the legend entry for ODE 3,5 are given as labels([3,5])
labels = {  'STAT3c',       'RJ2Lp-STAT3c',  'STAT3cp',      'RJ2Lp-STAT3cp', 'PPX', ...
            'PPX-STAT3cp',  'STAT3cpd',     'STAT3npd',     'STAT3np',      'PPN', ...
            'PPN-STAT3np',  'STAT3n',       'mRNAn-SEAP',   'mRNAc-SEAP',   'SEAPer', ...
            'SEAPg',        'SEAPex',       'Ligand',       'RJ',           'RJL', ...
            'RJ2L',         'RJ2Lp'};

%% Defining initial conditions
% The values bellow are in [nM]
Y0_STAT3c   = 1000;
Y0_3zero    = zeros(1,3);
Y0_PPX      = 50;
Y0_4zero    = zeros(1,4);
Y0_PPN      = 60;
Y0_7zero    = zeros(1,7);
Y0_L        = 16.9/100;       % The optimum value for the ligand is 1/100 the receptor, found via a Bell curve.
Y0_RJ       = 16.9;
Y0_3zero;
%     1         2-4        5       6-9       10     11-17      18    19     20-22 
Y0 = [Y0_STAT3c, Y0_3zero, Y0_PPX, Y0_4zero, Y0_PPN, Y0_7zero, Y0_L, Y0_RJ, Y0_3zero];