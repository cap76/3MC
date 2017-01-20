function [X,Y,Z,Energy_WLC] = FixedEnd(Parameters,X0x,Y0x,Z0x);
%Generates an arbtrary length Markov chain for chromosomes with ends fixed.
%A single update is defined: a crankshaft rotation about an axis connecting
%two randomly chosen beads. Example of use: FixedEndsScript.

addpath('./Functions','./netlab')
%__________________________________________________________________________
%Physical parameters
len    = Parameters(1,2); % Length of chromosome (No. of Khun segments).
k      = Parameters(1,3); % Rigidity constant.
Radius = Parameters(1,4); % Radius of NP.
a      = Parameters(1,6); % Khun segment length.
U      = Parameters(1,5); % Confinement potential.
beta   = 1;               % Inverse temperature.
%__________________________________________________________________________
%Unphysical parameters
No          = Parameters(1,1);    % No. steps in Markov chain.
%__________________________________________________________________________
%Empty arrays for storing data points during MC.
thet        = rand(No,1)*2*pi;
Bond3       =    randint(No,1,2)*(len-1) +1;
Bond4       =    randint(No,1,len-2)+2;
lenB        = randint(No,1,len-2)+2;
stepc       = rand(1,No);
X           = zeros(No,len);
Y           = zeros(No,len);
Z           = zeros(No,len);
En          = zeros(1,1);
Energy      = zeros(No,1);
Energy_WLC  = zeros(No,1);
%__________________________________________________________________________
%Set first state in the Markov chain equal to initial conditions.
X(1,:)  = X0x;
Y(1,:)  = Y0x;
Z(1,:)  = Z0x;
%__________________________________________________________________________
%Calculate bending energies.
Delta       = sqrt(diff(X0x).^2 + diff(Y0x).^2 + diff(Z0x).^2);
E1          = [(diff(X0x)./Delta)',(diff(Y0x)./Delta)',(diff(Z0x)./Delta)'];
En(1,1)     = sum(diag(E1*E1',1));
%__________________________________________________________________________
%Calculation of energy due to confinement.
ll   = length(find(X(1,:).^2 + Y(1,:).^2 + Z(1,:).^2 <= (Radius+1)^2));
E_c  = ll*0 + (len-ll)*U;
%__________________________________________________________________________
%Total starting energy.
Energy(1,1)     = -(k/a)*En + E_c + 3*sum(diff(X0x).^2 + diff(Y0x).^2 + diff(Z0x).^2)/(2*a^2);
Energy_WLC(1,1) = -(k/a)*En;
%__________________________________________________________________________
%Begin Metropolis-Hastings step.
for i = 2:No

maxi  = len - lenB(i,1);
Bond  = randint2(1,1,maxi)+1;
Bond2 = Bond + lenB(i,1);
%__________________________________________________________________________
%Change coordinates system.
dx      = X(i-1,Bond);
dy      = Y(i-1,Bond);
dz      = Z(i-1,Bond);
Xsub    = [X(i-1,Bond:Bond2)'-dx,Y(i-1,Bond:Bond2)'-dy,Z(i-1,Bond:Bond2)'-dz];
usub    = [X(i-1,Bond2) - dx,Y(i-1,Bond2) - dy,Z(i-1,Bond2) - dz];
udiv    = sqrt((X(i-1,Bond2) - dx).^2 + (Y(i-1,Bond2) - dy).^2 + (Z(i-1,Bond2) - dz).^2);
[Pr]    = Quaternion3(thet(i,1),usub/udiv,Xsub);
%__________________________________________________________________________
%Set new coordinates for the Markov chain step i.
X(i,:) = X(i-1,:);
Y(i,:) = Y(i-1,:);
Z(i,:) = Z(i-1,:);
X(i,Bond:Bond2) = Pr(:,1)+dx;
Y(i,Bond:Bond2) = Pr(:,2)+dy;
Z(i,Bond:Bond2) = Pr(:,3)+dz;
%__________________________________________________________________________
rx      = diff(X(i,:));
ry      = diff(Y(i,:));
rz      = diff(Z(i,:));
rx1      = diff(X(i-1,:));
ry1      = diff(Y(i-1,:));
rz1      = diff(Z(i-1,:));
Delta   = sqrt(rx.^2 + ry.^2 + rz.^2);
Delta1   = sqrt(rx1.^2 + ry1.^2 + rz1.^2);
%__________________________________________________________________________
%Calculate bond energies.
   E1 = [[(rx(1,max(Bond-1,1):min(Bond,len-1))./Delta(1,max(Bond-1,1):min(Bond,len-1))), ...
        (rx(1,Bond2-1:min(Bond2,len-1))./Delta(1,Bond2-1:min(Bond2,len-1)))]', ...
        [(ry(1,max(Bond-1,1):min(Bond,len-1))./Delta(1,max(Bond-1,1):min(Bond,len-1))), ...
        (ry(1,Bond2-1:min(Bond2,len-1))./Delta(1,Bond2-1:min(Bond2,len-1)))]', ...
        [(rz(1,max(Bond-1,1):min(Bond,len-1))./Delta(1,max(Bond-1,1):min(Bond,len-1))), ...
        (rz(1,Bond2-1:min(Bond2,len-1))./Delta(1,Bond2-1:min(Bond2,len-1)))]'];
   E2 = [[(rx1(1,max(Bond-1,1):min(Bond,len-1))./Delta1(1,max(Bond-1,1):min(Bond,len-1))), ...
        (rx1(1,Bond2-1:min(Bond2,len-1))./Delta1(1,Bond2-1:min(Bond2,len-1)))]', ...
        [(ry1(1,max(Bond-1,1):min(Bond,len-1))./Delta1(1,max(Bond-1,1):min(Bond,len-1))), ...
        (ry1(1,Bond2-1:min(Bond2,len-1))./Delta1(1,Bond2-1:min(Bond2,len-1)))]', ...
        [(rz1(1,max(Bond-1,1):min(Bond,len-1))./Delta1(1,max(Bond-1,1):min(Bond,len-1))), ...
        (rz1(1,Bond2-1:min(Bond2,len-1))./Delta1(1,Bond2-1:min(Bond2,len-1)))]'];

    
En1(1,1) = sum(diag(E1*E1',1));
En2(1,1) = sum(diag(E2*E2',1));
%__________________________________________________________________________
%Calculate total energy of the new state.
ll   = length(find((X(i,:).^2 + Y(i,:).^2 + Z(i,:).^2 <= (Radius+1)^2)));
E_c  = (ll)*0 + ((len) -ll)*U;
Energy_WLC(i,1) = Energy_WLC(i-1,1) + (k/a)*En2 - (k/a)*En1;
Energy(i,1) = Energy_WLC(i,1) + E_c + 3*sum(rx.^2 + ry.^2 +rz.^2)/(2*a^2);
%Set conditions for acceptance.
condition = rand(1,1);
dE        = abs(Energy(i,1) - Energy(i-1,1));
%Accept or reject new state.
if ((Energy(i,1)<Energy(i-1,1)) || (condition<exp(-beta*dE)))
else
X(i,:)  = X(i-1,:);
Y(i,:)  = Y(i-1,:);
Z(i,:)  = Z(i-1,:);
Energy(i,1)     = Energy(i-1,1);
Energy_WLC(i,1) = Energy_WLC(i-1,1);
end
end

%subsample       = 1:50:100000;
%X = X(subsample,:); 
%Y = Y(subsample,:); 
%Z = Z(subsample,:); 
%Energy_WLC = Energy_WLC(subsample,1);
