%% Configurations
% get a reusable default for Y0
Defaults;

options = odeset('RelTol',1e-9,'nonnegative',1);
Time = 60*[0:0.1:60*2]; % [sec]. it is now up to 2 hours

% due to floating point Precision we cant use a == b, 
% so instead we use abs(a - mean(b)) < some small value.
Precision = 1e-5;
% below negligible is a substance considert nonexistant
negligible = 1e-3;

%% The tests
% these should stay constant:       %nL is the name for L
nR =        'Receptors';    R =     [2,2,4,4,19:22, 21, 22]; %2&4 & 21&22 contain 2 receptors
nSTAT3 =    'STAT3';        STAT3 = [1:4, 6:9, 11,12, 7,8 ]; % 7 & 8 are listed twice, as they are a dimer so contain 2x STAT3
nPPX =      'PPX';          PPX =   [5,6];
nPPN =      'PPN';          PPN =   [10,11];
nSEAP =     'SEAP';         SEAP =  [15:17];

Test_isConstant(Y0,Time,options, Precision, R,      nR,     2)
Test_isConstant(Y0,Time,options, Precision, STAT3,  nSTAT3, 3)
Test_isConstant(Y0,Time,options, Precision, PPX,    nPPX,   4)
Test_isConstant(Y0,Time,options, Precision, PPN,    nPPN,   5)
%Total amount of SEAP should remain constant
Test_isSEAPConstant(Y0,Time,options, Precision, SEAP,   nSEAP,   17)

Test_6_PPX_suppreses_STAT3cpd_morePPX(Y0,Time,options,negligible);
Test_7_PPX_suppreses_STAT3cpd_PPXbindingConstant(Y0,Time,options,negligible);
Test_9_mRNA_gives_SEAP(Y0,Time,options,negligible);
Test_10_NoDeDimerisation_STAT3nd_more_mRNA(Y0,Time,options,negligible);
Test_11_NoPPX_more_mRNA(Y0,Time,options,negligible);
Test_12_weakPPX_more_mRNA(Y0,Time,options,negligible);
Test_13_NoPPN_more_mRNA(Y0,Time,options,negligible);
Test_14_weakPPN_more_mRNA(Y0,Time,options,negligible);

% plots
PLOT_noRJ2Lp_check_R(Y0,Time,options,negligible);
PLOT_STAT(Y0,Time,options,negligible);


%% The test functions
function Test_isSEAPConstant(Y0,Time,options, P, ODE, name, testnr)
    %Test SEAP case with a starting SEAP concentration
    ii = int64(16);
    Y0(15) = 5;
    
    [~,Y] = ode15s( @(t,y) ODEs(t,y,ii), Time, Y0, options);
    
    all_ODE = sum(Y(:,ODE),2);
    % test the condition:
    if all( abs(all_ODE - mean(all_ODE) < P) )
        disp(['[ Info ] Test ',num2str(testnr),' passed, the sum of ',name,' is constant.'])
    else
        % the test failed:
        figure(testnr)
        plot(Time,all_ODE)
        title(['The sum of all ',name,' should stay constant'])
        legend(['The sum of all ',name])
        warning('[ TEST FAILED ] Test %d: The sum of all %s should be constant!',testnr, name)
    end
end

function Test_isConstant(Y0,Time,options, P, ODE, name, testnr)
    [~,Y] = ode15s( @(t,y) ODEs(t,y,[]), Time, Y0, options);
    
    all_ODE = sum(Y(:,ODE),2);
    % test the condition:
    if all( abs(all_ODE - mean(all_ODE) < P) )
        disp(['[ Info ] Test ',num2str(testnr),' passed, the sum of ',name,' is constant.'])
    else
        % the test failed:
        figure(testnr)
        plot(Time,all_ODE)
        title(['The sum of all ',name,' should stay constant'])
        legend(['The sum of all ',name])
        warning('[ TEST FAILED ] Test %d: The sum of all %s should be constant!',testnr, name)
    end
end

function Test_6_PPX_suppreses_STAT3cpd_morePPX(Y0,Time,options, N)
    % increaseing the amount of PPX should prevent STAT3cpd formation
    Y0(5) = 1e20; % add lots of PPX

    [~,Y]=ode15s( @(t,y) ODEs(t,y,[]) ,Time,Y0,options);
    
    STAT3cpd = Y(:,7);

    % test the condition:
    if all( STAT3cpd < N)
        disp(['[ Info ] Test ',num2str(6),' passed.'])
    else
        % the test failed:
        figure(6)
        plot(Time,STAT3cpd,'-b',Time,0*STAT3cpd+N,'--r')
        title('exess PPX should supress STAT3cpd formation')
        legend('STAT3cpd','MAX allowed')
        warning('[ TEST FAILED ] Test 6: exess PPX does not supress STAT3cpd formation!')
    end
end

function Test_7_PPX_suppreses_STAT3cpd_PPXbindingConstant(Y0,Time,options, N)
    % increaseing the binding constant of PPX should prevent STAT3cpd formation
    
    ii = int64(7); % overide the binding constant
    
    [~,Y]=ode15s( @(t,y) ODEs(t,y,ii) ,Time,Y0,options);
    
    STAT3cpd = Y(:,7);

    % test the condition:
    if all( STAT3cpd < N )
        disp(['[ Info ] Test ',num2str(7),' passed.'])
    else
        % the test failed:
        figure(7)
        plot(Time,STAT3cpd,'-b',Time,0*STAT3cpd+N,'--r')
        title('stronger PPX binding should supress STAT3cpd formation')
        legend('STAT3cpd','MAX allowed')
        warning('[ TEST FAILED ] Test 7: stronger PPX binding does not supress STAT3cpd formation!')
    end
end

function Test_9_mRNA_gives_SEAP(Y0,Time,options,N)
    % mRNA should give SEAP secretion
    Y0(:)  =0;  % everything else is 0
    Y0(13) = 1; % mRNA is 1 nM
    
    [~,Y]=ode15s( @(t,y) ODEs(t,y,[]) ,Time,Y0,options);
    
    % 13 to 17 are all components of mRNA to the secretion pathway
    mRNAtoSEAP = Y(:,13:17);
    
    % test the condition:
    if any( mRNAtoSEAP(17) > N)
        disp(['[ Info ] Test ',num2str(9),' passed.'])
    else
        %else the test failed:
        figure(9)
        plot(Time,mRNAtoSEAP)
        title('mRNA should give SEAP')
        legend('mRNAn-SEAP','mRNAc-SEAP','SEAPer','SEAPg','SEAPex')
        warning('[ TEST FAILED ] Test 9: mRNA should result in SEAP excretion!')
    end
end

function Test_10_NoDeDimerisation_STAT3nd_more_mRNA(Y0,Time,options,N)
    % no de-dimersiation of STAT3npd should give more mRNA
    
    [~,Y]=ode15s( @(t,y) ODEs(t,y,[]) ,Time,Y0,options);
    ii = int64(10);
    [~,Y2]=ode15s( @(t,y) ODEs(t,y,ii) ,Time,Y0,options);
    
    
    Ignore_first_x_pc = 0.05;
    mRNAn_SEAP =   Y(floor(Ignore_first_x_pc*end):end ,13);
    mRNAn_SEAP2 = Y2(floor(Ignore_first_x_pc*end):end ,13);
    % adjust time aswell
    Time = Time(floor(Ignore_first_x_pc*end):end);
    
    % test the condition:
    if all(mRNAn_SEAP2 > mRNAn_SEAP)
        disp(['[ Info ] Test ',num2str(10),' passed.'])
    else
        % the test failed:
        figure(10)
        plot(Time,mRNAn_SEAP,'-r',Time,mRNAn_SEAP2,'-b')
        title('no de-dimersiation of STAT3npd should give more mRNA')
        legend('mRNAn-SEAP original','mRNAn-SEAP no-dedimerisation')
        warning('[ TEST FAILED ] Test 10: no de-dimersiation of STAT3npd should give more mRNA!')
    end
end

function Test_11_NoPPX_more_mRNA(Y0,Time,options,N)
    % no PPX should give more mRNA
    
    [~,Y]=ode15s( @(t,y) ODEs(t,y,[]) ,Time,Y0,options);
    
    Y0(5)=0;
    [~,Y2]=ode15s( @(t,y) ODEs(t,y,[]) ,Time,Y0,options);
    
    Ignore_first_x_pc = 0.05;
    mRNAn_SEAP =   Y(floor(Ignore_first_x_pc*end):end ,13);
    mRNAn_SEAP2 = Y2(floor(Ignore_first_x_pc*end):end ,13);
    % adjust time aswell
    Time = Time(floor(Ignore_first_x_pc*end):end);
    
    % test the condition:
    if all(mRNAn_SEAP2 > mRNAn_SEAP)
        disp(['[ Info ] Test ',num2str(11),' passed.'])
    else
        % the test failed:
        figure(11)
        plot(Time,mRNAn_SEAP,'-r',Time,mRNAn_SEAP2,'-b')
        title('no PPX should give more mRNA')
        legend('mRNAn-SEAP original','mRNAn-SEAP without PPX')
        warning('[ TEST FAILED ] Test 11: no PPX should give more mRNA!')
    end
end

function Test_12_weakPPX_more_mRNA(Y0,Time,options,N)
    %  weaker PPX binding should give more mRNA
    
    [~,Y]=ode15s( @(t,y) ODEs(t,y,[]) ,Time,Y0,options);
    
    ii = int64(12);
    [~,Y2]=ode15s( @(t,y) ODEs(t,y,ii) ,Time,Y0,options);
    
    Ignore_first_x_pc = 0.05;
    mRNAn_SEAP =   Y(floor(Ignore_first_x_pc*end):end ,13);
    mRNAn_SEAP2 = Y2(floor(Ignore_first_x_pc*end):end ,13);
    % adjust time aswell
    Time = Time(floor(Ignore_first_x_pc*end):end);
    
    % test the condition:
    if all(mRNAn_SEAP2 > mRNAn_SEAP)
        disp(['[ Info ] Test ',num2str(12),' passed.'])
    else
        % the test failed:
        figure(12)
        plot(Time,mRNAn_SEAP,'-r',Time,mRNAn_SEAP2,'-b')
        title('weaker PPX binding should give more mRNA')
        legend('mRNAn-SEAP original','mRNAn-SEAP weaker PPX')
        warning('[ TEST FAILED ] Test 12: weaker PPX binding should give more mRNA!')
    end
end

function Test_13_NoPPN_more_mRNA(Y0,Time,options,N)
    % no PPN should give more mRNA
    
    [~,Y]=ode15s( @(t,y) ODEs(t,y,[]) ,Time,Y0,options);
    
    Y0(10)=0;
    [~,Y2]=ode15s( @(t,y) ODEs(t,y,[]) ,Time,Y0,options);
    
    Ignore_first_x_pc = 0.05;
    mRNAn_SEAP =   Y(floor(Ignore_first_x_pc*end):end ,13);
    mRNAn_SEAP2 = Y2(floor(Ignore_first_x_pc*end):end ,13);
    % adjust time aswell
    Time = Time(floor(Ignore_first_x_pc*end):end);
    
    % test the condition:
    if all(mRNAn_SEAP2 > mRNAn_SEAP)
        disp(['[ Info ] Test ',num2str(13),' passed.'])
    else
        % the test failed:
        figure(13)
        plot(Time,mRNAn_SEAP,'-r',Time,mRNAn_SEAP2,'-b')
        title('no PPN should give more mRNA')
        legend('mRNAn-SEAP original','mRNAn-SEAP without PPN')
        warning('[ TEST FAILED ] Test 13: no PPN should give more mRNA!')
    end
end

function Test_14_weakPPN_more_mRNA(Y0,Time,options,N)
    %  weaker PPN binding should give more mRNA
    
    [~,Y]=ode15s( @(t,y) ODEs(t,y,[]) ,Time,Y0,options);
    
    ii = int64(14);
    [~,Y2]=ode15s( @(t,y) ODEs(t,y,ii) ,Time,Y0,options);
    
    Ignore_first_x_pc = 0.05;
    mRNAn_SEAP =   Y(floor(Ignore_first_x_pc*end):end ,13);
    mRNAn_SEAP2 = Y2(floor(Ignore_first_x_pc*end):end ,13);
    % adjust time aswell
    Time = Time(floor(Ignore_first_x_pc*end):end);
    
    % test the condition:
    if all(mRNAn_SEAP2 > mRNAn_SEAP)
        disp(['[ Info ] Test ',num2str(14),' passed.'])
    else
        % the test failed:
        figure(14)
        plot(Time,mRNAn_SEAP,'-r',Time,mRNAn_SEAP2,'-b')
        title('weaker PPN binding should give more mRNA')
        legend('mRNAn-SEAP original','mRNAn-SEAP weaker PPN')
        warning('[ TEST FAILED ] Test 14: weaker PPN binding should give more mRNA!')
    end
end

% the plot functions
function PLOT_noRJ2Lp_check_R(Y0,Time,options,N)
    %  noRJ2Lp_check_R
    ii = int64(15);
    [~,Y]=ode15s( @(t,y) ODEs(t,y,ii) ,Time,Y0,options);
    
    % no test
    all_R = Y(:,[18:21]);
    %

    figure(15)
    plot(Time,all_R)
    title('noRJ2Lp check R')
    legend('L','RJ','RJL','RJ2L')
end

function PLOT_STAT(Y0,Time,options,N)
    [~,Y]=ode15s( @(t,y) ODEs(t,y,[]) ,Time,Y0,options);
    
    % no test
    all_STAT = Y(:,[1:4, 6:9, 11,12 ]);
    %
    figure(16)
    plot(Time,all_STAT)
    title('STAT3')
    legend('STAT3c','RJ2Lp-STAT3','STAT3cp', 'RJ2Lp-STAT3p','PPX-STAT3cp','STAT3cpd','STAT3npd','STAT3np', 'PPN-STAT3np','STAT3n')
end