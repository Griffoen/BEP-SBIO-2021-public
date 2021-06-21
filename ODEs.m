%% The main model function
function [dy] = ODEs(t,y,i)
%{
This is the naming scheme for the constants:
K{}  - some constant
KTxy - a transport constant, where y gets transported to x. ex: KTcSTAT3
KB   - a binding constant.
KBu  - an unbinding constant
KM   - a constant related to Michaelis-Menten
KMH  - a hill kinetics exponent
KD   - degradation constant for mRNA
Kt   - translation rate
KP   - phosphorylation
%}

%% Variables currently used with their ODE number given by Y(x)
%--- Signal pathway
STAT3c          = y(1);
RJ2Lp_STAT3c    = y(2);
STAT3cp         = y(3);
RJ2Lp_STAT3cp   = y(4);
PPX             = y(5);
PPX_STAT3cp     = y(6);
STAT3cpd        = y(7);
STAT3npd        = y(8);
STAT3np         = y(9);
PPN             = y(10);
PPN_STAT3np     = y(11);
STAT3n          = y(12);
%--- SEAP production and secretion
mRNAn_SEAP      = y(13);
mRNAc_SEAP      = y(14);
SEAPer          = y(15);
SEAPg           = y(16);
SEAPex          = y(17);
%--- Receptor
L               = y(18);
RJ              = y(19);
RJL             = y(20);
RJ2L            = y(21);
RJ2Lp           = y(22);
%---

%% Constants
% Michaelis-Menten
KM1  = 0.01;  % [nM/s]
KM2  = 400;   % [nM]
KM3  = 6e-3;  % [1/s]
KM4  = 10;    % [nM]
Rtot = 1100;  % [nM] concentration ribosomes

% Hill exponent 
KMH1 = 1;     % [-] Used for mRNA transcription

% Transport
KTcSTAT3     = 0.05;    % [1/s]
KTnSTAT3cpd  = 0.005;   % [1/s]
KTcmRNA_SEAP = 0.001;   % [1/s]
KTerSEAP     = 0.0005;  % [1/s]
KTgSEAP      = 0.0005;  % [1/s]

% Binding
Kd1 = 0.1;      % [nM]       (range 1000-0.1 nM) % Kd of ligand/receptor binding
Kd2 = 1090;     % [nM]                           % Kd of the dimerization of the receptor units
KB1 = 0.008;    % [1/(nM*s)]
KB2 = 0.005;    % [1/(nM*s)]
KB3 = 0.02;     % [1/(nM*s)]
KB4 = 0.001;    % [1/(nM*s)]
KB5 = KB3;      % [1/(nM*s)]
KB6 = 0.001;    % [1/(nM*s)]
KB7 = 0.43;     % [1/(nM*s)] (diffusion controlled)
KB8 = 0.000312; % [1/(nM*s)]

% Unbinding
KBu1  = 0.8;     % [1/s]
KBu2  = 0.4;     % [1/s]
% KBu3  = nan;   % [1/s]  % KBu3 is unused
KBu4  = 0.5;     % [1/s]
KBu5  = 0.1;     % [1/s]
KBu6  = 0.2;     % [1/s]
KBu7  = 0.003;   % [1/s]
KBu8  = KBu5;    % [1/s]
KBu9  = 0.2;     % [1/s]
KBu10 = 0.005;   % [1/s]
KBu11 = KB7*Kd1; % [1/s]
KBu12 = KB8*Kd2; % [1/s]

% Degradation 
KDn1 = 0;        % [1/s]
KDc1 = 0.0005;   % [1/s]

% Phosphorylation  
KP1 = 0.005;     % [1/s]


%% Constant vatiation
% The code below varies constants for the purpose of testing the model
% response under various conditions. The first half is reserved for model
% validation and should not be modified or extended for other reasons. The 
% second half is for model variation and predictions and can be modified freely.

% If there is intention to change variables, I is non-empty
if not(isempty(i)) 
    if isinteger(i)
        %% RESERVED CODE: validation testing
        % This if-statement checks the type of i.
        % Matlab stores REAL numbers as a "double" by default.
        % this does not qualify as an integer for this test, as only 
        % int8, int16, int32, int64, uint8, uint16, uint32, and uint64 do.
        % the functionality tests in TEST_ODE.m use this to selectively
        % override constants for the purpose of model validation.
        switch i
            case 7
                % increase the binding constant from PPX to STAT3cp
                KB4 = 1e20;
            case 8
                 % decrease the unbinding constant from PPX-STAT3cp
                KBu6 = 0;
            case 10
                % zero the de-dimerisation constant of STATnpd
                KBu8 = 0;
            case 12
                % zero the binding of PPX to STAT3cp 
                KB4 = 0;
            case 14
                % zero the binding of PPN to STAT3np
                KB6 = 0;
            case 15
                % prevent conversion from RJ2L to RJ2Lp
                KP1 = 0;
            case 16
                % Prevent export of SEAP mRNA
                KTcmRNA_SEAP = 0;
        end
    else
        %% NON-RESERVED CODE: model variation
        % The code below can be modified or extended to include extra model
        % variations. For extension, create a new 'case' block with an unused
        % number. In calls to ODEs.m make sure the 'i' parameter is assigned a
        % list starting with that number.
        switch i(1)
            case 1001 % used for Plot_LigandvsReceptor.m
                % prevent STAT3p from binding to RJ2Lp
                KB2 = 0; 
            case 1002

                VAR = [1/10, 1, 10]*KBu4;
                KBu4 = VAR(i(2));
                
            case 1003
                
                KP1 = 0; % prevent conversion of RJ2L to RJ2Lp
                
            case 1004
                
                VAR = [1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3]*KB7;
                KB7 = VAR(i(2));
                
            case 1005
                
                VAR = [1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3]*Kd1;
                Kd1 = VAR(i(2));
                KBu11 = KB7 * Kd1;
                
            case 1006
                
                VAR = [0, 50, 100];
                PPX = VAR(i(2));
            case 1007
                
                VAR = 10.^[-3:0.25:1.5]*Kd2;
                Kd2 = VAR(i(2));
                KBu12 = KB8*Kd2;
        end
    end
end


%% The ODE's

%% [ 1 ] STAT3c
% 
% [+] Transportation of STAT3 from the nucleus to the cytoplasm
% [STAT3n] -> [STAT3c]
% 
% [+/-] Binding of STAT3 to the receptor complex
% [RJ2Lp]+[STAT3c] <-> [RJ2Lp-STAT3c]
%
% [+] Dissociation of the PPX-STAT3cp complex, releasing dephosphorylated STAT3
% [PPX-STAT3cp] -> [PPX]+[STAT3c]

D_STAT3c = KTcSTAT3*STAT3n ...
    - KB1*STAT3c*RJ2Lp + KBu1*RJ2Lp_STAT3c ...
    + KBu7*PPX_STAT3cp;

%% [ 2 ] RJ2Lp_STAT3c
%
% [+/-] Binding of STAT3 to the receptor complex
% [RJL2p]+[STAT3c] <-> [RJ2Lp-STAT3c]
%
% [-] Release of phosphorylated STAT3 from receptor
% [RJ2Lp-STAT3c] -> [RJ2Lp]+[STAT3cp]

D_RJ2Lp_STAT3c = KB1*STAT3c*RJ2Lp - KBu1*RJ2Lp_STAT3c ...
    - KBu2*RJ2Lp_STAT3c;

%% [ 3 ] STAT3cp
%
% [+] Release of phosphorylated STAT3 from receptor
% [RJ2Lp-STAT3c] -> [RJ2Lp]+[STAT3cp]
%
% [+/-] Binding of STAT3cp to the receptor complex
% [RJ2Lp]+[STAT3cp] <-> [RJ2Lp-STAT3cp]
%
% [+/-] Dimerisation of STAT3 in the cytoplasm
% 2[STAT3cp] <-> [STAT3cp-STAT3cp]
%
% [+/-] Binding of cellular phosphatase to STAT3cp
% [PPX]+[STAT3cp] <-> [PPX-STAT3cp]

D_STAT3cp = KBu2*RJ2Lp_STAT3c ...
    - KB2*STAT3cp*RJ2Lp + KBu4*RJ2Lp_STAT3cp ...
    - 2*KB3*STAT3cp*STAT3cp + 2*KBu5*STAT3cpd ...
    - KB4*PPX*STAT3cp + KBu6*PPX_STAT3cp;

%% [ 4 ] RJ2Lp_STAT3cp
%
% [+/-] Binding of STAT3cp to the receptor complex
% [RJ2Lp]+[STAT3cp] <-> [RJ2Lp-STAT3cp]

D_RJ2Lp_STAT3cp = KB2*STAT3cp*RJ2Lp - KBu4*RJ2Lp_STAT3cp;

%% [ 5 ] PPX
%
% [+] Dissociation of the PPX-STAT3cp complex, releasing dephosphorylated STAT3
% [PPX-STAT3cp] -> [PPX]+[STAT3c]
%
% [+/-] Binding of cellular phosphatase to STAT3cp
% [PPX]+[STAT3cp] <-> [PPX-STAT3cp]

D_PPX = KBu7*PPX_STAT3cp ...
    - KB4*PPX*STAT3cp + KBu6*PPX_STAT3cp;

%% [ 6 ] PPX_STAT3cp
%
% [+/-] Binding of cellular phosphatase to STAT3cp
% [PPX]+[STAT3cp] <-> [PPX-STAT3cp]
%
% [-] Dissociation of the PPX-STAT3cp complex, releasing dephosphorylated STAT3
% [PPX-STAT3cp] -> [PPX]+[STAT3c]

D_PPX_STAT3cp = KB4*PPX*STAT3cp - KBu6*PPX_STAT3cp ...
    - KBu7*PPX_STAT3cp;

%% [ 7 ] STAT3cpd
%
% [+/-] Dimerisation of STAT3 in the cytoplasm
% 2[STAT3cp] <-> [STAT3cpd]
% 
% [-] Transportation of STAT3pd from the cytoplasm to the nucleus
% [STAT3cpd] -> [STAT3npd]

D_STAT3cpd = KB3*STAT3cp*STAT3cp - KBu5*STAT3cpd ...
    - KTnSTAT3cpd*STAT3cpd;

%% [ 8 ] STAT3npd
%
% [+] Transportation of STAT3pd from the cytoplasm to the nucleus
% [STAT3cpd] -> [STAT3npd]
% 
% [+/-] Dimerisation of STAT3np
% 2[STAT3np] <-> [STAT3npd]
%

D_STAT3npd = KTnSTAT3cpd*STAT3cpd ...
    - KBu8*STAT3npd + KB5*STAT3np*STAT3np;

%% [ 9 ] STAT3np
%
% [+/-] Dimerisation of STAT3np
% 2[STAT3np] <-> [STAT3npd]
%
% [+/-] Binding of nuclear phosphatase to STAT3p
% [PPN]+[STAT3np] <-> [PPN-STAT3np]

D_STAT3np = 2*KBu8*STAT3npd - 2*KB5*STAT3np*STAT3np ...
    - KB6*PPN*STAT3np + KBu9*PPN_STAT3np;

%% [ 10 ] PPN
%
% [+] Disassociation of the PPN-STAT3p complex, releasing dephosphorylated STAT3
% [PPN-STAT3np] -> [PPN]+[STAT3n]
%
% [+/-] Binding of nuclear phosphatase to STAT3p
% [PPN]+[STAT3np] <-> [PPN-STAT3np]

D_PPN = KBu10*PPN_STAT3np ...
    - KB6*PPN*STAT3np + KBu9*PPN_STAT3np;

%% [ 11 ] PPN_STAT3np
% 
% [+/-] Binding of nuclear phosphatase to STAT3p
% [PPN]+[STAT3np] <-> [PPN-STAT3np]
%  
% [-] Disassociation of the PPN-STAT3np complex, releasing dephosphorylated STAT3
% [PPN-STAT3np] -> [PPN]+[STAT3n]

D_PPN_STAT3np = KB6*PPN*STAT3np - KBu9*PPN_STAT3np ...
    - KBu10*PPN_STAT3np;

%% [ 12 ] STAT3n
% 
% [+] Disassociation of the PPN-STAT3np complex, releasing dephosphorylated STAT3
% [PPN-STAT3np] -> [PPN]+[STAT3n]
% 
% [-] Transportation of STAT3 from the nucleus to the cytoplasm 
% [STAT3n] -> [STAT3c]

D_STAT3n = KBu10*PPN_STAT3np ...
    - KTcSTAT3*STAT3n;

%% [ 13 ] mRNAn_SEAP
%
% [+] Transcription of SEAPmRNA 
% d[mRNAn]/dt= KM1*[STAT3npd]/(KM2+[STAT3npd])
%
% [-] Transportation of SEAPmRNA from the nucleus to the cytoplasm
% [mRNAn] -> [mRNAc]
%
% [-] degradation
% [mRNAn] -> [-]

D_mRNAn_SEAP = KM1*STAT3npd^KMH1 /( KM2^KMH1 + STAT3npd^KMH1 ) ...
    - KTcmRNA_SEAP*mRNAn_SEAP ...
    - KDn1*mRNAn_SEAP;

%% [ 14 ] mRNAc_SEAP
% 
% [+] Transportation of SEAPmRNA from the nucleus to the cytoplasm
% [mRNAn] -> [mRNAc]
% 
% [-] Degradation of mRNA in the cytoplasm
% [mRNAc] -> [-]

D_mRNAc_SEAP = KTcmRNA_SEAP*mRNAn_SEAP ...
    - KDc1*mRNAc_SEAP;

%% [ 15 ] SEAPer
% 
% [+] Translation of mRNA to SEAP, via Michaelis Menten
% d[SEAPer]/dt = KM3*[mRNAc_SEAP]*Rtot/(KM4 + Rtot)
% 
% [-] transport of SEAP from ER to golgi
% [SEAPer] -> [SEAPg]

D_SEAPer = KM3*mRNAc_SEAP*Rtot/(KM4 + Rtot) ... 
    - KTerSEAP*SEAPer;

%% [ 16 ] SEAPg
% 
% [+] Transport of SEAP from ER to golgi
% [SEAPer] -> [SEAPg]
%
% [-] Transport of SEAP from golgi to extracellular
% [SEAPg] -> [SEAPex]

D_SEAPg = KTerSEAP*SEAPer ...
    - KTgSEAP*SEAPg;

%% [ 17 ] SEAPex
%
% [+] Transport of SEAP from golgi to extracellular
% [SEAPg] -> [SEAPex]

D_SEAPex = KTgSEAP*SEAPg;

%% [ 18 ] L
% 
% [+/-] Binding of ligand to Receptor-Jak
% [RJ]+[L] <-> [LRJ]

D_L = 0; %- KB7*RJ*L + KBu11*RJL;

%% [ 19 ] RJ
%
% [+/-] Binding of ligand to Receptor-Jak
% [RJ]+[L] <-> [LRJ]
% 
% [+/-] Binding of second Receptor to Ligand-Receptor-Jak
% [LRJ]+[RJ] <-> [RJ2L]

D_RJ = - KB7*RJ*L + KBu11*RJL ...
    - KB8*RJL*RJ + KBu12*RJ2L;

%% [ 20 ] RJL
%
% [+/-] Binding of ligand to Receptor-Jak
% [RJ]+[L] <-> [LRJ]
% 
% [+/-] Binding of second Receptor to Ligand-Receptor-Jak
% [LRJ]+[RJ] <-> [RJ2L]

D_RJL = KB7*RJ*L - KBu11*RJL ...
    - KB8*RJL*RJ + KBu12*RJ2L;

%% [ 21 ] RJ2L
% 
% [+/-] Binding of second Receptor to Ligand-Receptor-Jak
% [LRJ]+[RJ] <-> [RJ2L]
% 
% [-] Phosphorylation of the receptor complex
% [RJ2L] -> [RJ2Lp]

D_RJ2L = KB8*RJL*RJ - KBu12*RJ2L ...
    - KP1*RJ2L;

%% [ 22 ] RJ2Lp
%
% [+] Phosphorylation of the receptor complex
% [RJ2L] -> [RJ2Lp]
% 
% [+] Release of phosphorylated STAT3 from receptor
% [RJ2Lp-STAT3c] -> [RJ2Lp]+[STAT3cp]
% 
% [+/-] Binding of STAT3 to the receptor complex
% [RJ2Lp]+[STAT3c] <-> [RJ2Lp-STAT3c]
% 
% [+/-] Binding of STAT3cp to the receptor complex
% - competitieve inhibitie nog even checken
% [RJ2Lp]+[STAT3cp] <-> [RJ2Lp-STAT3cp]

D_RJ2Lp = KP1*RJ2L ...
    + KBu2*RJ2Lp_STAT3c ...
    - KB1*STAT3c*RJ2Lp + KBu1*RJ2Lp_STAT3c...
    - KB2*STAT3cp*RJ2Lp + KBu4*RJ2Lp_STAT3cp;

%% Collecting the results with their ODE number
dy = [
    D_STAT3c        % = 1;
    D_RJ2Lp_STAT3c  % = 2;
    D_STAT3cp       % = 3;
    D_RJ2Lp_STAT3cp % = 4;
    D_PPX           % = 5;
    D_PPX_STAT3cp   % = 6;
    D_STAT3cpd      % = 7;
    D_STAT3npd      % = 8;
    D_STAT3np       % = 9;
    D_PPN           % = 10;
    D_PPN_STAT3np   % = 11;
    D_STAT3n        % = 12;
    %--- SEAP production and secretion ---%
    D_mRNAn_SEAP    % = 13;
    D_mRNAc_SEAP    % = 14;
    D_SEAPer        % = 15;
    D_SEAPg         % = 16;
    D_SEAPex        % = 17;
    %--- Receptor ---%
    D_L             % = 18;
    D_RJ            % = 19;
    D_RJL           % = 20;
    D_RJ2L          % = 21;
    D_RJ2Lp         % = 22;
    ];
end