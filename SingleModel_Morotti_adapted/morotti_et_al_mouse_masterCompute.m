%% Main file: morotti_et_al_mouse_masterCompute
% 
% Lyon et al. 2021. ADAPTED FROM MOROTTI MODEL AVAILABLE WITH PUBLICATION BELOW
% Please cite the following paper when using this model:
% Morotti S, Edwards AG, McCulloch AD, Bers DM & Grandi E. (2014). A novel
% computational model of mouse myocyte electrophysiology to assess the 
% synergy between Na(+) loading and CaMKII. Journal of Physiology, doi:
% 10.1113/jphysiol.2013.266676.
% 
% This model was built upon the code of the Soltis and Saucerman model
% of rabbit ventricular EC coupling.
% Reference: Soltis AR & Saucerman JJ. (2010). Synergy between CaMKII 
% substrates and beta-adrenergic signaling in regulation of cardiac
% myocyte Ca(2+) handling. Biophysical Journal 99, 2038-2047.
% 
% This file loads initial conditions, calls the ode solver and plots
% simulation results.

% close all;
% clear all;
% clearvars -except cleftca RyRaff cleftcaval clamp nb_beat srca
clearvars
clc;
tic
%% Parameters for external modules
PKP2 = 1; BARS = 0; %flag for presence of PKP2 and beta adrenergic stimulation
nb_beat = 2; clamp = 0; %number of beat to simulate and flag for clamp recording
fname = 'yfinal_PKP2_1Hz'; %initial conditions

% ECC and CaM modules
% freq = 1;                   % [Hz] - Stimulation frequency
% cycleLength = 1e3/freq;     % [ms]
cycleLength = 1000; %350
CaMtotDyad = 418;           % [uM]
BtotDyad = 1.54/8.293e-4;   % [uM]
CaMKIItotDyad = 120;        % [uM]
CaNtotDyad = 3e-3/8.293e-4; % [uM]
PP1totDyad = 96.5;          % [uM]
CaMtotSL = 5.65;            % [uM]
BtotSL = 24.2;              % [uM]
CaMKIItotSL = 120*8.293e-4; % [uM]
CaNtotSL = 3e-3;            % [uM]
PP1totSL = 0.57;            % [uM]
CaMtotCyt = 5.65;           % [uM]
BtotCyt = 24.2;             % [uM]
CaMKIItotCyt = 120*8.293e-4;% [uM]
CaNtotCyt = 3e-3;           % [uM] 
PP1totCyt = 0.57;           % [uM]

% ADJUST CaMKII ACTIVITY LEVELS (expression = 'WT', 'OE', or 'KO')
expression = 'WT'%;

CKIIOE = 0; % Should be zero during 'WT' and 'KO' runs
if strcmp(expression,'OE') % OE
    CKIIOE = 1; % Flag for CaMKKII-OE (0=WT, 1=OE)
    n_OE=6;
    CaMKIItotDyad = 120*n_OE;         % [uM] 
    CaMKIItotSL = 120*8.293e-4*n_OE;  % [uM]
    CaMKIItotCyt = 120*8.293e-4*n_OE; % [uM]
elseif strcmp(expression,'KO')
    CaMKIItotDyad = 0;          % [uM] 
    CaMKIItotSL = 0;            % [uM]
    CaMKIItotCyt = 0;           % [uM]
end

%plb_val=38; % RABBIT
plb_val=106; % MOUSE

% Parameters for CaMKII module
LCCtotDyad = 31.4*.9;       % [uM] - Total Dyadic [LCC] - (umol/l dyad)
LCCtotSL = 0.0846;          % [uM] - Total Subsarcolemmal [LCC] (umol/l sl)
if PKP2, LCCtotDyad = 31.4*.7; end
RyRtot = 382.6;             % [uM] - Total RyR (in Dyad)
PP1_dyad = 95.7;            % [uM] - Total dyadic [PP1]
PP1_SL = 0.57;              % [uM] - Total Subsarcolemmal [PP1]
PP2A_dyad = 95.76;          % [uM] - Total dyadic PP2A
OA = 0;                     % [uM] - PP1/PP2A inhibitor Okadaic Acid
PLBtot = plb_val;           % [uM] - Total [PLB] in cytosolic units

% Parameters for BAR module
if BARS == 0, Ligtot = 0; else Ligtot = 0.1; end %.1               % [uM] - SET LIGAND CONCENTRATION (0 or 0.1)
LCCtotBA = 0.025;           % [uM] - [umol/L cytosol]
RyRtotBA = 0.135;           % [uM] - [umol/L cytosol]
PLBtotBA = plb_val;         % [uM] - [umol/L cytosol]
TnItotBA = 70;              % [uM] - [umol/L cytosol]
IKstotBA = 0.025;           % [uM] - [umol/L cytosol]
ICFTRtotBA = 0.025;         % [uM] - [umol/L cytosol]
PP1_PLBtot = 0.89;          % [uM] - [umol/L cytosol]
IKurtotBA = 0.025;          % [uM] - [umol/L cytosol] MOUSE
PLMtotBA = 48;              % [uM] - [umol/L cytosol] MOUSE

% For Recovery from inactivation of LTCC
recoveryTime = 10; % initialize to smallest value

% Parameter varied in protocol simulation
variablePar = 20; % initilization
%% Collect all parameters and define mass matrix for BAR module

p = [cycleLength,recoveryTime,variablePar,CaMtotDyad,BtotDyad,CaMKIItotDyad,...
    CaNtotDyad,PP1totDyad,CaMtotSL,BtotSL,CaMKIItotSL,CaNtotSL,PP1totSL,...
    CaMtotCyt,BtotCyt,CaMKIItotCyt,CaNtotCyt,PP1totCyt,...
    LCCtotDyad,RyRtot,PP1_dyad,PP2A_dyad,OA,PLBtot,LCCtotSL,PP1_SL,...
    Ligtot,LCCtotBA,RyRtotBA,PLBtotBA,TnItotBA,IKstotBA,ICFTRtotBA,...
    PP1_PLBtot,IKurtotBA,PLMtotBA,CKIIOE];

%% Establish and define global variables

global tStep tArray I_Ca_store I_to_store I_Na_store I_K1_store ibar_store 
global gates Jserca IKs_store Jleak ICFTR Incx Jrel
global I_ss_store dVm_store Ipca_store I_NaK_store I_Nabk_store I_kr_store
global I_kur1_store I_kur2_store
tStep = 1;
tArray = zeros(1,1e6);
I_Ca_store=zeros(1,1e6);
I_to_store=zeros(3,1e6);
I_Na_store = zeros(1,1e6);
I_K1_store = zeros(1,1e6);
ibar_store=zeros(1,1e6);
gates = zeros(2,1e6);
Jserca = zeros(1,1e6);
IKs_store = zeros(1,1e6);
Jleak = zeros(1e6,2);
Jrel  = zeros(1,1e6);
ICFTR = zeros(1,1e6);
Incx = zeros(1,1e6);
I_kur1_store = zeros(1,1e6);
I_kur2_store = zeros(1,1e6);
I_ss_store = zeros(1,1e6);
dVm_store = zeros(1,1e6);
Ipca_store = zeros(1,1e6);
I_NaK_store = zeros(1,1e6);
I_Nabk_store = zeros(1,1e6);
I_kr_store = zeros(1,1e6);

%% Assign initial conditions
% WT %%%%%%%%%%%%%%%%%%%%%%%
load(fname)
% load yfin_WT_1Hz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y0n=yfinal;
if clamp, y0n(36) = cleftca; y0n(31)=srca; end %steadystate

%% Run simulation
% tic
% tend = 1e3;
tend = nb_beat*cycleLength;
tspan = [0 tend];
load('testJacobian.mat'); %precomputed Jacobian to speed up code 
options = odeset('RelTol',1e-5,'MaxStep',2, 'JPattern',S);
[t,y] = ode15s(@morotti_et_al_mouse_masterODEfile,tspan,y0n,options,p,clamp); %%steadystate 
yfinal = y(end,:)';
% toc
%% Output variables

tArray = tArray(1:tStep);
Ica = I_Ca_store(1:tStep);
Ito = I_to_store(1,1:tStep);
Itof = I_to_store(2,1:tStep);
Itos = I_to_store(3,1:tStep);
INa = I_Na_store(1:tStep);
IK1 = I_K1_store(1:tStep);
s1 = gates(1,1:tStep);
k1 = gates(2,1:tStep);
Jserca = Jserca(1:tStep);
Iks = IKs_store(1:tStep);
Jleak = Jleak(1:tStep,:);
Jrel = Jrel(1:tStep);
ICFTR = ICFTR(1:tStep);
Incx = Incx(1:tStep);
Ikur1 = I_kur1_store(1:tStep);
Ikur2 = I_kur2_store(1:tStep);
Iss = I_ss_store(1:tStep);
dVm = dVm_store(1:tStep);
Ipca = Ipca_store(1:tStep);
INaK = I_NaK_store(1:tStep);
INabk = I_Nabk_store(1:tStep);
Ikr = I_kr_store(1:tStep);

toc

%% Plot traces
figure
subplot(2,2,1), plot(t,y(:,38),'LineWidth',2); ylabel('[Ca]i (mM)'), hold on
subplot(2,2,2), plot(t,y(:,37),'LineWidth',2); ylabel('[Ca]sarcolemmal (mM)'), hold on
subplot(2,2,3), plot(t,y(:,36),'LineWidth',2); ylabel('[Ca]cleft (mM)'), hold on
