addpath('func','models');
addpath FrontiersPaper_code/PopulationOfModels_Morotti/func
addpath FrontiersPaper_code/SingleModel_Morotti_adapted 


paramNames = {'GNa', 'GbarNal', 'ICa_scale', 'GtoFast',...
    'Gkr', 'Gk1', 'IbarNCX', 'IbarNaK', 'fSerca','fRyR'};
% Name for the population of models
POMpath = ['POM_' fnametype '.mat'];

       
%% Create other population from initial 
nameend = 'PKP2_CABx2.5_from240421_500beats'; %name of final population
pop = 'POM_PKP2_CABx1.5_from240421.mat'; % inital population to start from
POMpath2 = ['POM_' nameend];

PKP2flag = 1; experimflag = 0; BARSflag = 0; CabFac = 2.5; 

POM = generatePoM_others(3,1000,500, pop,PKP2flag,experimflag,BARSflag,CabFac); %
save(POMpath2,'POM');

savecurrentsfname = ['myPoM_' nameend];
load(['POM_' nameend '.mat'])
for i=1:length(POM)
save_currents(i,POM,savecurrentsfname,PKP2flag,experimflag,BARSflag,CabFac)   
end
