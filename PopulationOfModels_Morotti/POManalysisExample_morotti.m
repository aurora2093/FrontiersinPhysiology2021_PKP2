% Lyon et al. 2021. Adapted from below: 
% This script is an example of how to create a POM, calibrate it to some
% experimental data and obtain results from the calibration process. The
% script uses the TP06 model to create the POM, only a small population
% (1000 models) is created at a single cycle length (600 ms). This POM is
% then calibrated to three different experimental data:
% 'ARIdata/time_point_1.mat', 'ARIdata/time_point_2.mat',
% 'ARIdata/time_point_3.mat'. The population is saved in 'POM.mat' and the
% ePOMs are saved in 'ePOMs/ePOM1.mat', 'ePOMs/ePOM2.mat',
% 'ePOMs/ePOM3.mat'.
%
% Several outputs are saved in 'Results/':
%   'Boxplots/chanX.png' contains the boxplots of the statistical
%   distributions of the ePOMs' parameters. Each image contains the results
%   of calibrating the POM to channel X at of the three time points.
%
% This script was developed in the Multiscale Cardiovascular Engineering 
% Group at University College London (MUSE-UCL) by Carlos Ledezma. It's 
% protected by a Creative Commons Attribution-ShareAlike 4.0 International 
% license (https://creativecommons.org/licenses/by-sa/4.0/).
%
% For further details please refer to the published work. This paper must
% be cited if using any of the routines or data provided with this script:
% C. Ledezma et al. "Bridging organ and cellular-level behavior in ex-vivo
% experimental platforms using populations of models of cardiac
% electrophysiology" (2018). ASME Journal of Engineering and Science in
% Medical Diagnostics and Therapy. doi: 10.1115/1.4040589.

addpath('func','models');
addpath FrontiersPaper_code/PopulationOfModels_Morotti/func
addpath FrontiersPaper_code/SingleModel_Morotti_adapted 

paramNames = {'GNa', 'GbarNal', 'ICa_scale', 'GtoFast',...
    'Gkr', 'Gk1', 'IbarNCX', 'IbarNaK', 'fSerca','fRyR'};

fnametype = ['WT2'];

% Significance level to use in the statistical analysis
pMax = 0.001;

% Name for the population of models
POMpath = ['POM_' fnametype '.mat'];

% Paths to the ARI signals
ARIpath = {'ARIdata/time_point_1.mat',....
           'ARIdata/time_point_2.mat',...
           'ARIdata/time_point_3.mat' };
       
% Create folder to save the ePOMs
ePOMroot = 'ePOMs';
if isempty(dir(ePOMroot))
    mkdir(ePOMroot);
end
% Paths to the ePOMs
ePOMpath = {'ePOMs/ePOM1.mat', 'ePOMs/ePOM2.mat', 'ePOMs/ePOM3.mat'};

% Create folders to save the results
if isempty(dir('Results'))
    mkdir('Results')
end
if isempty(dir('Results/Boxplots'))
    mkdir('Results/Boxplots')
end
if isempty(dir(['Results/MWUtestResults' num2str(pMax*100)]))
    mkdir(['Results/MWUtestResults' num2str(pMax*100)]);
end

% Create a population of TP06 models and save it
% Note that this is just an example, for the parameters required to have a
% sufficiently larged, converged, population refer to the paper

POM = generatePoM(3,1000,200,500,0.5,1.5); % TO USE FOR THE MOROTTI MODEL 
save(POMpath,'POM');
movefile('MyPoM',['MyPoM_' fnametype]) 


%% Get all variables for all model populations

nameend = 'PKP2_CABx1.5_from240421_500beats'; %Example name of previously generated population of models
savecurrentsfname = ['myPoM_' nameend];
load(['POM_' nameend '.mat'])
PKP2flag = 1; experimflag = 1; BARSflag = 0; 
for i=1:length(POM)
save_currents(i,POM,savecurrentsfname,PKP2flag,experimflag,BARSflag)   
end


