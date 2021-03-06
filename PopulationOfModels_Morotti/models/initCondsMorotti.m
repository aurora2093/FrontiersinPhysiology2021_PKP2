function [y0n,p] = initCondsMorotti(PKP2flag,BARSflag)

PKP2 = PKP2flag;
% nb_beat = 1; 
% ECC and CaM modules
freq = 1;                   % [Hz] - CHANGE DEPENDING ON FREQUENCY
cycleLength = 1e3/freq;     % [ms]
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
Ligtot = 0;
if BARSflag == 1,Ligtot = 0.1, end%.1               % [uM] - SET LIGAND CONCENTRATION (0 or 0.1)
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
load F:\Documents\BME\BME_CircAdapt\PKP2_project_cell\SIMULATIONS\morotti_code\yfin_WT_1Hz
% load yfinal_normal_3000ms_pace
% load yfinal_normal_3000ms_pace_BAR
% load yfinal_normal_3000ms_pace_ICabx1.5.mat
% load 'yfinal_normal_3000ms_pace_RyR30%.mat'
% load yfinal_normal_350ms_pace
% load yfinal_normal_350ms_pace_BAR
% load yfinal_PKP2_3000ms_pace
% load yfinal_PKP2_3000ms_pace_ICabx1.5.mat
% load yfinal_PKP2_3000ms_pace_BAR
% load yfinal_PKP2_3000ms_pace_ICabx1.5_BAR.mat
% load F:\Documents\BME\BME_CircAdapt\PKP2_project_cell\SIMULATIONS\morotti_code\yfinal_normal_3000ms_pace.mat
% load yfinal_PKP2_3000ms_pace_ICabx1.5_BAR_verticilide.mat
% load yfinal_PKP2_3000ms_pace_ICabx1.5_BAR_verticilide_passiveleakonly.mat
% load 'yfinal_PKP2_3000ms_pace_RyR30%_v2_2603.mat'
% load yfinal_PKP2_350ms_pace
% load yfinal_PKP2_350ms_pace_BAR
% load yfinal_PKP2_350ms_pace_ICabx2.mat
% load 'yfinal_PKP2_350ms_pace_RyR30%.mat'
% load yfin_PKP2_350mspaced
% load yfin_WT_1Hz_120sISO
%load yfin_WT_NaGain_1Hz
% CaMKII-OE %%%%%%%%%%%%%%%%
%load yfin_OE_1Hz
%load yfin_OE_loop_1Hz
%load yfin_OE_NoNaGain_1Hz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y0n=yfinal;