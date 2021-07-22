
rate1 = []; diastCa = [];
% POM = DADmodel_RyR;
for i =1:length(POM)
    

fname = ['myPoM_PKP2_BARS_CABx2.5_RyR20_from240421_500beats/PoM' num2str(i) '_currents.mat'];
load(fname)
t = time; 
    
Ca = y(:,38); resampCa = resample(Ca,t,'spline'); cleft = y(:,36);
gradient = diff(resampCa); gradient = gradient(10:400);
rate1 = [rate1 min(gradient(find(gradient<0)))];
   
diastCa = [diastCa min(cleft)];
end

for i =1:length(POM)
    
fname = ['myPoM_PKP2_BARS_CABx2.5_RyR20_NCX2_from240421_500beats/PoM' num2str(i) '_currents.mat'];
load(fname)
t = time; 
    
Ca = y(:,38); resampCa = resample(Ca,t,'spline'); cleft = y(:,36);
gradient = diff(resampCa); gradient = gradient(10:400);
rate1 = [rate1 min(gradient(find(gradient<0)))];
diastCa = [diastCa min(cleft)];
    
end

g = [zeros(1,length(POM)) ones(1,length(POM))];
boxplot([rate1(1:length(POM)) rate1(length(POM)+1:end)],g)
[h,p] = ttest2(rate1(1:length(POM)), rate1(length(POM)+1:end)); p

g = [zeros(1,length(POM)) ones(1,length(POM))];
boxplot([diastCa(1:length(POM)) diastCa(length(POM)+1:end)],g)
[h,p] = ttest2(diastCa(1:length(POM)), diastCa(length(POM)+1:end)); p


rate2 = [];
for i =1:length(POM)
    
fname = ['myPoM_PKP2_BARS_CABx2.5_RyR20_from240421/PoM' num2str(i) '_currents.mat'];
load(fname)
t = time; 
    
Ca = y(:,38); resampCa = resample(Ca,t,'spline');
[a,b] = max(resampCa(1:400)); Capart = resampCa(b:end);
l = find(Capart < min(resampCa)+min(resampCa)*5/100,1,'first');
    
rate2 = [rate2 l];
end
for i =1:length(POM)
    
fname = ['myPoM_PKP2_BARS_CABx2.5_RyR20_NCX2_from240421/PoM' num2str(i) '_currents.mat'];
load(fname)
t = time; 
    
Ca = y(:,38); resampCa = resample(Ca,t,'spline');
[a,b] = max(resampCa(1:400)); Capart = resampCa(b:end);
l = find(Capart < min(resampCa)+min(resampCa)*5/100,1,'first');
    
rate2 = [rate2 l];
end

g = [zeros(1,length(POM)) ones(1,length(POM))];
boxplot([rate2(1:length(POM)) rate2(length(POM)+1:end)],g)
[h,p] = ttest2(rate2(1:length(POM)), rate2(length(POM)+1:end)); p





ncx = [];
for i =1:length(POM)
    
fname = ['myPoM_PKP2_BARS_CABx1.5_from240421/PoM' num2str(i) '_currents.mat'];
load(fname)
t = time; 
    
Ca = y(:,38); resampCa = resample(Ca,t,'spline');
[a,b] = max(resampCa(1:400)); Capart = resampCa(b:end);
l = find(Capart < min(resampCa)+min(resampCa)*5/100,1,'first');
    
ncx = [ncx min(Incx(:,:))];
end
for i =1:length(POM)
    
fname = ['myPoM_PKP2_BARS_CABx2.5_from240421/PoM' num2str(i) '_currents.mat'];
load(fname)
t = time; 
    
Ca = y(:,38); resampCa = resample(Ca,t,'spline');
[a,b] = max(resampCa(1:400)); Capart = resampCa(b:end);
l = find(Capart < min(resampCa)+min(resampCa)*5/100,1,'first');
    
ncx = [ncx min(Incx(:,:))];
end

g = [zeros(1,length(POM)) ones(1,length(POM))];
boxplot([ncx(1:length(POM)) ncx(length(POM)+1:end)],g)
[h,p] = ttest2(ncx(1:length(POM)), ncx(length(POM)+1:end)); p