PoM = struct('V',[],'t',[],'params',[],'APD',[],'CL',[],'yfin',[]);
CL = 1000;

k = 1;
for i = 1 : 500
    allPhysio = true;
    load(['myPoM_WT2/PoM' num2str(i) '.mat']); % Load the record
    for j = 1 : length(CL) % Loop through the cycle lengths to check that all are physiological
        if ~isPhysio(var{1}{j},var{2}{j},var{4}(j))
            allPhysio = false;
        end
    end
    
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

figure
for i=1:length(PoM)
    plot(cell2mat(PoM(i).t),cell2mat(PoM(i).V)), hold on
end