
load('POM_WT_calibrated_240421.mat')


threshold = 10*10^-6; threshold2 = 0.004; ampthreshold = 1*10^-4;

%% BAR PLOTS

% pop1 = 'myPoM_WT_from240421_500beats';
% pop2 = 'myPoM_WT_BARS_from240421_500beats';
% pop3 = 'myPoM_PKP2_from240421_500beats';
% pop4 = 'myPoM_PKP2_BARS_from240421_500beats';
% pop5 = 'myPoM_PKP2_CABx1.5_from240421_500beats';
% pop6 = 'myPoM_PKP2_BARS_CABx1.5_from240421_500beats';
% pop7 = 'myPoM_PKP2_CABx2.5_from240421_500beats';
% pop8 = 'myPoM_PKP2_BARS_CABx2.5_from240421_500beats';

pop1 = 'myPoM_PKP2_BARS_from240421_500beats';
pop2 = 'myPoM_PKP2_BARS_CABx1.5_from240421_500beats';
pop3 = 'myPoM_PKP2_BARS_CABx1.5_RyR20_from240421_500beats';
pop4 = 'myPoM_PKP2_BARS_CABx2.5_from240421_500beats';
pop5 = 'myPoM_PKP2_BARS_CABx2.5_RyR20_from240421_500beats';
pop6 = 'myPoM_PKP2_BARS_CABx2.5_RyR20_NCX2_from240421_500beats';
% pop1 = 'myPoM_PKP2_from240421_500beats';
% pop2 = 'myPoM_PKP2_BARS_from240421_500beats';
% pop3 = 'myPoM_PKP2_BARS_CABx2.5_from240421_500beats';
POP = {pop1 pop2 pop3 pop4 pop5 pop6};

Ca_amp = zeros(length(POP),length(POM)); CaRel_occ = zeros(length(POP),length(POM)); 
DAD_occ = zeros(length(POP),length(POM)); cleftCa = zeros(length(POP),length(POM));

for P = 1:length(POP)
    pop = POP{P};
    for i =1:length(POM)
        fname = [pop '/PoM' num2str(i) '_currents.mat'];
        load(fname)
        Ca_amp(P,i) = max(y(:,38));  cleftCa(P,i) = min(y(:,36));
        dCa = diff(y(:,38))./diff(time); dV = diff(y(:,39))./diff(time);
        x = find(time>400,1,'first');
        % if ~isempty(find(dCa(x,:)>threshold)), CaRel_occ(1,i) = 1; end
        Ca = y(:,38); if ~isempty(find(Ca(x:end,:)>min(Ca)+ampthreshold)), CaRel_occ(P,i) = 1; end
        if ~isempty(find(dV(x,:)>threshold2)), DAD_occ(P,i) = 1; end
    end
end

DADmodel_RyR = []; noDADmodel_RyR = [];
for i =1:length(POM)
    fname = ['myPoM_PKP2_BARS_CABx2.5_RyR20_from240421_500beats/PoM' num2str(i) '_currents.mat'];
    load(fname)
    Ca_amp(length(POP)+1,i) = max(y(:,38)); cleftCa(length(POP)+1,i) = min(y(:,36));
    dCa = diff(y(:,38))./diff(time); dV = diff(y(:,39))./diff(time);
    x = find(time>400,1,'first');
    Ca = y(:,38);
    if ~isempty(find(Ca(x:end,:)>min(Ca)+ampthreshold))
        CaRel_occ(length(POP)+1,i) = 1; DADmodel_RyR = [DADmodel_RyR i];
    else
        noDADmodel_RyR = [noDADmodel_RyR i];
    end
    if ~isempty(find(dV(x,:)>threshold2)), DAD_occ(length(POP)+1,i) = 1; end
end

for i =1:length(POM)
    if ismember(i,DADmodel_RyR)
        fname = ['myPoM_PKP2_BARS_CABx2.5_RyR20_NCX2_from240421_500beats/PoM' num2str(i) '_currents.mat'];
        load(fname)
    else
        fname = ['myPoM_PKP2_BARS_CABx2.5_RyR20_from240421_500beats/PoM' num2str(i) '_currents.mat'];
        load(fname)
    end
    load(fname)
    Ca_amp(length(POP)+2,i) = max(y(:,38)); cleftCa(length(POP)+2,i) = min(y(:,36));
    dCa = diff(y(:,38))./diff(time); dV = diff(y(:,39))./diff(time);
    x = find(time>400,1,'first');
    if ~isempty(find(dCa(x:end,:)>threshold)), CaRel_occ(length(POP)+2,i) = 1; end
    Ca = y(:,38); if ~isempty(find(Ca(x:end,:)>min(Ca)+ampthreshold)), CaRel_occ(length(POP)+2,i) = 1; end
    if ~isempty(find(dV(x,:)>threshold2)), DAD_occ(length(POP)+2,i) = 1; end
end


S = sum(CaRel_occ,2)./length(POM).*100; 
figure, bar(S([1 2 4 5 6]))

figure, bar(sum(DAD_occ,2))

%% DAD MODEL ANALYSIS

DADmodel = [1,2,3,4,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,26,28,29,30,31,32,33,34,35,36,38,40,41,42,44,45,47,48,49,50,51,52,53,55,58,59,60];
DADmodel_RyR = [1,2,4,6,7,9,10,11,12,13,15,16,19,20,21,25,28,30,32,34,38,40,42,43,44,47,48,51,53,55,58,59,60];
noDADmodel_RyR = [3,5,8,14,17,18,22,23,24,26,27,29,31,33,35,36,37,39,41,45,46,49,50,52,54,56,57,61];
[sharedvals,idx] = intersect(DADmodel,DADmodel_RyR,'stable');

figure
[y0n,p] = initCondsMorotti(1,1); RyRtot = p(20); LCCtotDyad = p(19); diasCa = [];
% for i = sharedvals(1:10)
for i = DADmodel_RyR(11)
    fname = ['myPoM_PKP2_BARS_CABx2.5_RyR20_NCX2_from240421/PoM' num2str(i) '_currents.mat'];
    load(fname)
    t = time;
    ax1 = subplot(4,4,1),hold on, plot(t,y(:,39),'-','LineWidth',2); ylabel('Em (mV)');
    ax2 = subplot(4,4,2),hold on,plot(t,100*y(:,87+45+4)./RyRtot,'-','LineWidth',2); ylabel('RyRp (%)'), hold on
    ax3 = subplot(4,4,3),hold on,plot(t,100*y(:,87+45+2)./LCCtotDyad,'-','LineWidth',2); ylabel('LTCCp (%)');
    ax4 = subplot(4,4,4),hold on,plot(t,y(:,38),'-','LineWidth',2); ylabel('[Ca]i (mM)');
    ax5 = subplot(4,4,5),hold on,plot(t,y(:,30)+y(:,31),'-','LineWidth',2); ylabel('tot [Ca]SR (mM)'); xlabel('Time (s)');
    ax6 = subplot(4,4,6),hold on, plot(t, y(:,31),'-','LineWidth',2); ylabel('Free SR Ca(mM)'); xlabel('Time (s)');
    xx = find(tArray>60,1,'first'); xxx=find(tArray>1000,1,'first');
    ax7 = subplot(4,4,7),hold on, plot(tArray(1:xx), Jleak(1:xx,1),'-','LineWidth',2); ylabel('RyR Jleak Systolic (pA/pF)'); xlabel('Time (s)');
    %     ax8 = subplot(4,4,8),hold on, plot(tArray(xxx:end), Jleak(xxx:end,1),'LineWidth',2); ylabel('RyR Jleak Diastolic (pA/pF)'); xlabel('Time (s)');
    ax8 = subplot(4,4,8),hold on, plot(t, y(:,37),'-','LineWidth',2); ylabel('Ca sl'); xlabel('Time (s)');
    ax9 = subplot(4,4,9),hold on, plot(tArray, Jserca(:,:),'-','LineWidth',2); ylabel('JSerca'); xlabel('Time (s)');
    ax10 = subplot(4,4,10),hold on, plot(t, y(:,46),'-','LineWidth',2); ylabel('Background Ca(mM)'); xlabel('Time (s)');
    ax11 = subplot(4,4,11),hold on, plot(tArray, Incx(:,:),'-','LineWidth',2); ylabel('Incx'); xlabel('Time (s)');
    ax12 = subplot(4,4,12),hold on, plot(t, y(:,36),'-','LineWidth',2); ylabel('Cleft Ca(mM)'); xlabel('Time (s)');
    ax13 = subplot(4,4,13), plot(t,y(:,14),'-','LineWidth',2); ylabel('RyR R'), hold on
    ax14 = subplot(4,4,14), plot(t,y(:,15),'-','LineWidth',2); ylabel('RyR O'), hold on
    ax15 = subplot(4,4,15), plot(t,y(:,16),'-','LineWidth',2); ylabel('RyR I'), hold on
    %     ax16 = subplot(4,4,16), plot(tArray,Jrel,'-','LineWidth',2); ylabel('Jrel'), hold on
    ax16 = subplot(4,4,16), plot(t, y(:,45),'-','LineWidth',2); ylabel('NCXvar'), hold on
    linkaxes([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 ax10 ax11 ax12 ax13 ax14 ax15 ax16],'x')
    
    %     subplot(2,2,1), plot(t,y(:,4)), hold on
    %     subplot(2,2,2), plot(t,y(:,2)), hold on
    %     diasCa = [diasCa min(y(:,38))];
end

boxplot_conductivity;

Ca_decay;

%% Diastolic cleft Ca // DAD occurence

%1:normal 2PKP2 3PKP2BARS 4PKP2BARSRYR 5PKP2BARSCAB 6PKP2BARSCABRYR
%7PKP2BARSCABRYRNCX
figure


for i =[1 2 3 4 5 6]
    D=[]; C = []; erro= [];
    C = [C mean(cleftCa(i,:),2)]; erro = [erro std(cleftCa(i,:))];
    D = [D length(find(CaRel_occ(i,:)>0))/length(POM)*100];
    errorbar(C,D,erro,'o','LineWidth',3,'MarkerSize',7), hold on
end

A = [C' D']; corrcoef(A); 
%% PLOT THE RYR STATES FIGURE

figure
% i = 39;
% i = 1;
% i = 44
i = 15;
fname = ['myPoM_WT_from240421/PoM' num2str(i) '_currents.mat'];
fname = ['myPoM_PKP2_BARS_CABx1.5_from240421_500beats/PoM' num2str(i) '_currents.mat'];

load(fname)
t = time;

ax1 = subplot(3,6,1),hold on,plot(t,y(:,38),'-','LineWidth',2); ylabel('[Ca]i (mM)');

ax2 = subplot(3,6,2), plot(t,y(:,14),'-','LineWidth',2); ylabel('RyR R'), hold on
ax3 = subplot(3,6,3), plot(t,y(:,15),'-','LineWidth',2); ylabel('RyR O'), hold on
ax4 = subplot(3,6,4), plot(t,y(:,16),'-','LineWidth',2); ylabel('RyR I'), hold on
ax5 = subplot(3,6,5),hold on, plot(t, y(:,36),'-','LineWidth',2); ylabel('Cleft Ca(mM)'); xlabel('Time (s)');
ax6 = subplot(3,6,6), plot(tArray,Jrel,'-','LineWidth',2); ylabel('Jrel'), hold on

linkaxes([ax1 ax2 ax3 ax4 ax5 ax6],'x')

