% params = generateSamples(popSize, meanParams, lowMult*meanParams, highMult*meanParams ,...
%     popSize, 2, []);
function PoM = generatePoM_others(model,CL,numBeats,initpop,PKP2flag,experimflag,BARSflag,CabFac)

% !!!!!!!!! MAKE SURE YOU ADJUTS THE PKP2/EXPERIM FLAGS  !!!!!!!!!!! 
addpath('models','func')

% PKP2flag = 1; experimflag = 1; BARSflag = 1;
mod = @myMorotti;
% load('POM_WT.mat')
load(initpop);
savePath = 'myPoM_other/'; 

numIt = length(POM);

for i = 1 : numIt
    result(i) = struct('V',{{0}},'time',{{0}},'APD',[],'yfin',{{0}});
end

if isempty(gcp('nocreate'))
    parpool();
end

% Temporary directory to save PoM pieces
if isempty(dir('myPoM_other'))
    mkdir('myPoM_other')
end

POMSparams = zeros(length(POM),length(POM(1).params)); POMSyfin = zeros(length(POM),length(POM(1).yfin));
for i = 1:length(POM)
    POMSparams(i,:) = POM(i).params; 
    POMSyfin(i,:) = POM(i).yfin; 
end 
% params = zeros(1,length(POM(1).params)); X0 = zeros(length(POM(1).yfin),1);
parfor_progress(numIt);
% Solve model for every combination of parameters
parfor i = 1 : numIt
    for j = 1 : length(CL)
        if model == 1
            X0 = initCondsORd;
        elseif model == 2
            X0 = initCondsTP06();
        elseif model ==3 
            [X0,p] = initCondsMorotti(PKP2flag, BARSflag);
        else
            error('Model has to be 1 (ORd) or 2 (TP06)');
        end
%         currentPOM = POM(i); 
%         params = currentPOM.params;
%         X0 = currentPOM.yfin'; 
        params = POMSparams(i,:);
        X0 = POMSyfin(i,:)'; 
        S = getJacobian;
%         options = [];
        options = odeset('RelTol',1e-5,'MaxStep',2, 'JPattern',S);
%         odeset('MaxStep',1,'RelTol',1e-7,'AbsTol',1e-9);
        for n=1:numBeats
%             [result(i).time{j} , X]=ode15s(mod,[0 CL(j)],X0,options,1,params(i,:),celltype);
%               [result(i).time{j} , X]=ode15s(mod,[0 CL(j)],X0,options,p,params(i,:));
              [result(i).time{j} , X]=ode15s(mod,[0 CL(j)],X0,options,p,params,PKP2flag,experimflag,CabFac);
            X0=X(size(X,1),:);
        end
        
%         result(i).V{j} = X(:,1);
        result(i).V{j} = X(:,39);
        result(i).APD(j) = findAPD(result(i).V{j},result(i).time{j},90);
        result(i).yfin = X0;     
    end
    
    % Save each model
    V = result(i).V;
    APD = result(i).APD;
    time = result(i).time;
    yfin = result(i).yfin;
    parsave([savePath 'PoM' num2str(i) '.mat'],{V , time , params , APD, CL, yfin});
    parfor_progress;
end
parfor_progress(0);

clear result

disp('Saving population to a single variable...')
PoM = struct('V',[],'t',[],'params',[],'APD',[],'CL',[],'yfin',[]);

k = 1;
for i = 1 : numIt
    allPhysio = true;
    load(['myPoM_other/PoM' num2str(i) '.mat']); % Load the record    
    if allPhysio % If everything is in order, save the action potentials
        PoM(k).V = var{1};
        PoM(k).t = var{2};
        PoM(k).params = var{3};
        PoM(k).APD = var{4};
        PoM(k).CL = var{5};
        PoM(k).yfin = var{6};
        k = k + 1;
    end
%     delete(['myPoM/PoM' num2str(i) '.mat']);
end

delete(gcp);
end 