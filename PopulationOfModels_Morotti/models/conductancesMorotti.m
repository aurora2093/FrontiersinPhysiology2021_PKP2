function cond = conductancesMorotti()
%endo = 0, epi = 1, M = 2
GNa = 10;

GbarNal = .0065*(1+0)*2; % deltGbarNal_CKII = 0

ICa_scale = 5.25;

IbarNCX = 1;

IbarNaK = 5;

% ks = 25;
% 
% kleak = 2*5.348e-6; 
% 
% Kmf = 0.3e-3; % [mM] changed from rabbit (0.246e-3) % from Yang-Saucerman
% Kmr = 2.1;

fSerca = 1; 
fRyR = 1; 

GtoFast=0.44;

Ko = 5.4;
Gkr = 0.03*sqrt(Ko/5.4);
Gk1 = 0.3*sqrt(Ko/5.4);

cond = [GNa, GbarNal, ICa_scale, GtoFast, Gkr, Gk1, IbarNCX, IbarNaK, fSerca, fRyR];