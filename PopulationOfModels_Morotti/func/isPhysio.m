function isPhysio = isPhysio(V,t,APD)
% Find out if an AP is physiological

% Derivative of V
dV = diff(V)./diff(t);

Vpart = V(find(t>80,1,'first'):end,:);

% Upstroke duration
UPD = t(find(dV < 0, 1));

% if APD > 100 && APD < 400 && max(V) > 0 && min(V) < -64 && UPD < 10
% if APD > 20 && APD < 80 && max(V) > 25 && min(V) < -75 && UPD < 10 && all(Vpart <min(V)+1)
if APD > 20 && APD < 80 && max(V) > 25 && min(V) < -75 && UPD < 10 && all(Vpart <min(V)+1.9)
    isPhysio = true;        
else
    isPhysio = false;
end
