function [X,Y,Z,Energy_WLC] = BouquetTV(varargin);
%
%Description
%
%[X,Y,Z,Energy] = BOUQUET(Parameters,x,y,z,VMF) generates an arbtrary 
%length Markov Chain for chromosomes in a Bouquet formation with boundary 
%conditions drawn from a von Mises-Fisher distribution. Three types of 
%update are defined: (i) a crankshaft rotation about an axis connecting 
%two randomly determined bonds; (ii) a crankshaft rotation of one arm; 
%and (iii) rotation of chromosomes via random axis. Example of use:
%BouquetScript.

addpath('./Functions','./netlab')

switch nargin
case 0
	warning('No inputs 1')
case 1
	if length(varargin{1})==6
		disp('No initial conformation input. Generating Conformation now.')
		Parameters = varargin{1};
		[x0,y0,z0,theta,phi]    = tMultivMF(0,0,1,1e-9,Parameters(1,4),1);
		[x1,y1,z1]              = RandFlightCon([x0,y0,z0],[x0,y0,z0],Parameters(1,6),Parameters(1,2));
		X0x = x1; Y0x = y1; Z0x = z1;
        VMF = [1e-9,0,0,1];
        disp('Initial Configuration Generated.')
	else
		warning('Not enough parameters specified in input 2')
	end

case 2
	if length(varargin{1})==6 && length(varargin{2})==4
		Parameters = varargin{1};
		VMF = varargin{2}
		[x0,y0,z0,theta,phi]    = tMultivMF(0,0,1,1e-9,Parameters(1,4),1);
		[x1,y1,z1]              = RandFlightCon([x0,y0,z0],[x0,y0,z0],Parameters(1,6),Parameters(1,2));
		X0x = x1; Y0x = y1; Z0x = z1;
		disp('No initial conformation input. Generating Conformation now.')
	else
		warning('Not enough parameters specified in one or more inputs')
	end

case 3
		warning('Wrong number of inputs specified 3')
case 4
		warning('Wrong number of inputs specified 4')
case 5
				Parameters = varargin{1};
				X0x = varargin{2}; 
				Y0x = varargin{3}; 
				Z0x = varargin{4};
				VMF = varargin{5};
                
end

%__________________________________________________________________________
%Physical parameters
len    = Parameters(1,2); % Length of chromosome (No. of Khun segments).
Radius = Parameters(1,3); % Radius of NP.
a      = Parameters(1,5); % Khun segment length.
U      = Parameters(1,4); % Confinement potential.
k      = Parameters(1,6:end)'; % Rigidity constant.

beta   = 1;               % Inverse temperature.
%__________________________________________________________________________
%Unphysical parameters
No          = Parameters(1,1);    % No. steps in Markov chain.
mu          = VMF(1,2:4);         % Mean cluster direciton.
Xmu1        = [X0x(1,1),Y0x(1,1),Z0x(1,1)];
Xmu2        = [X0x(1,len),Y0x(1,len),Z0x(1,len)];
kappa       = VMF(1,1);           % Cluster strength.
%__________________________________________________________________________
%Empty arrays for storing data points during MCMC.
thet            = rand(No,1)*2*pi;
[xxx,yxx,zxx]   = Points_on_Sphere(1,No);
uu              = [xxx',yxx',zxx'];
Bond3           = randint2(No,1,2)*(len-1) +1;
Bond4           = randint2(No,1,len-2)+2;
lenB            = randint2(No,1,len-2)+2;
X               = zeros(No,len);
stepc           = rand(1,No);
Y               = zeros(No,len);
Z               = zeros(No,len);
En              = zeros(1,1);
Energy          = zeros(No,1);
Energy_WLC      = zeros(No,1);
PvMF1           = zeros(No,1);
PvMF2           = zeros(No,1);
%__________________________________________________________________________
%Set first state in the Markov chain equal to initial conditions.
X(1,:)  = X0x;
Y(1,:)  = Y0x;
Z(1,:)  = Z0x;
%__________________________________________________________________________
%Calculate bending energies.
Delta       = sqrt(diff(X0x).^2 + diff(Y0x).^2 + diff(Z0x).^2);
E1          = [(diff(X0x)./Delta)',(diff(Y0x)./Delta)',(diff(Z0x)./Delta)'];
En(1,1)     = sum(k.*diag(E1*E1',1));
%__________________________________________________________________________
%Calculation of energy due to confinement.
ll   = length(find(X(1,:).^2 + Y(1,:).^2 + Z(1,:).^2 <= (Radius+1)^2));
E_c  = ll*0 + (len-ll)*U;
%__________________________________________________________________________
%Total starting energy.
Energy(1,1)     = -(1/a)*En + E_c + 3*sum(diff(X0x).^2 + diff(Y0x).^2 + diff(Z0x).^2)/(2*a^2);
Energy_WLC(1,1) = -(1/a)*En;
PvMF1(1,1)       = (kappa*(mu*Xmu1'./sqrt(sum(Xmu1.^2))));
PvMF2(1,1)       = (kappa*(mu*Xmu2'./sqrt(sum(Xmu2.^2))));
%__________________________________________________________________________
%Begin Metropolis-Hastings step.
for i = 2:No
stepchoice = stepc(1,i);  
if stepchoice < (1/3)
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
E1          = [(rx./Delta)',(ry./Delta)',(rz./Delta)'];
En(1,1)     = sum(k.*diag(E1*E1',1));
  %  keyboard
%    En1(1,1) = sum(k(max(Bond-1,1):min(Bond,len-1)).*diag(E1*E1',1));
%En2(1,1) = sum(k(max(Bond-1,1):min(Bond,len-1)).*diag(E2*E2',1));
%__________________________________________________________________________

%Calculation of energy due to confinement.
ll   = length(find(X(1,:).^2 + Y(1,:).^2 + Z(1,:).^2 <= (Radius+1)^2));
E_c  = ll*0 + (len-ll)*U;
%__________________________________________________________________________
%Total starting energy.
Energy(i,1)     = -(1/a)*En + E_c + 3*sum(diff(X0x).^2 + diff(Y0x).^2 + diff(Z0x).^2)/(2*a^2);
Energy_WLC(i,1) = -(1/a)*En;
PvMF1(i,1)       = (kappa*(mu*Xmu1'./sqrt(sum(Xmu1.^2))));
PvMF2(i,1)       = (kappa*(mu*Xmu2'./sqrt(sum(Xmu2.^2))));
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

elseif stepchoice > (2/3)
Bond = Bond3(i,1); %randint(1,1,2)*(len-1) +1;
Bond2 = Bond4(i,1);%randint(1,1,len-2)+2;
%__________________________________________________________________________
dx = X(i-1,Bond2);
dy = Y(i-1,Bond2);
dz = Z(i-1,Bond2);
Xsub = [X(i-1,min(Bond,Bond2):max(Bond,Bond2))'-dx, ...
        Y(i-1,min(Bond,Bond2):max(Bond,Bond2))'-dy, ...
        Z(i-1,min(Bond,Bond2):max(Bond,Bond2))'-dz];

usub = [X(i-1,Bond2),Y(i-1,Bond2),...
        Z(i-1,Bond2)]/sqrt(X(i-1,Bond2).^2+Y(i-1,Bond2).^2+ ...
        Z(i-1,Bond2).^2);
[Pr] = Quaternion3(thet(i,1),usub,Xsub);
%__________________________________________________________________________
% Set new coordinates for the Markov chain step i.
X(i,:) = X(i-1,:);
Y(i,:) = Y(i-1,:);
Z(i,:) = Z(i-1,:);
X(i,min(Bond,Bond2):max(Bond,Bond2)) = Pr(:,1) + dx;
Y(i,min(Bond,Bond2):max(Bond,Bond2)) = Pr(:,2) + dy;
Z(i,min(Bond,Bond2):max(Bond,Bond2)) = Pr(:,3) + dz;
%__________________________________________________________________________
rx = diff(X(i,:));
ry = diff(Y(i,:));
rz = diff(Z(i,:));
rx1 = diff(X(i-1,:));
ry1 = diff(Y(i-1,:));
rz1 = diff(Z(i-1,:));
Delta = sqrt(rx.^2 + ry.^2 + rz.^2);
Delta1 = sqrt(rx1.^2 + ry1.^2 + rz1.^2);
%__________________________________________________________________________
E1          = [(rx./Delta)',(ry./Delta)',(rz./Delta)'];
En(1,1)     = sum(k.*diag(E1*E1',1));
  %  keyboard
%    En1(1,1) = sum(k(max(Bond-1,1):min(Bond,len-1)).*diag(E1*E1',1));
%En2(1,1) = sum(k(max(Bond-1,1):min(Bond,len-1)).*diag(E2*E2',1));
%__________________________________________________________________________

%Calculation of energy due to confinement.
ll   = length(find(X(1,:).^2 + Y(1,:).^2 + Z(1,:).^2 <= (Radius+1)^2));
E_c  = ll*0 + (len-ll)*U;
%__________________________________________________________________________
%Total starting energy.
Energy(i,1)     = -(1/a)*En + E_c + 3*sum(diff(X0x).^2 + diff(Y0x).^2 + diff(Z0x).^2)/(2*a^2);
Energy_WLC(i,1) = -(1/a)*En;
%PvMF1(i,1)       = (kappa*(mu*Xmu1'./sqrt(sum(Xmu1.^2))));
%PvMF2(i,1)       = (kappa*(mu*Xmu2'./sqrt(sum(Xmu2.^2))));
%__________________________________________________________________________
% Calculate total energy of the new state.
%ll   = length(find((X(i,:).^2 + Y(i,:).^2 + Z(i,:).^2 <= (Radius+1)^2)));
%E_c  = (ll)*0 + ((len) -ll)*U;
%Energy_WLC(i,1) = Energy_WLC(i-1,1) + (1/a)*En2 - (1/a)*En1;
%Energy(i,1) = Energy_WLC(i,1) + E_c + 3*sum(rx.^2 + ry.^2 +rz.^2)/(2*a^2);
Xmu1=[X(i,1),Y(i,1),Z(i,1)];
Xmu2=[X(i,len),Y(i,len),Z(i,len)];
PvMF1(i,1) = (kappa *(mu*Xmu1'./sqrt(sum(Xmu1.^2))));
PvMF2(i,1) = (kappa *(mu*Xmu2'./sqrt(sum(Xmu2.^2))));
%Set conditions for acceptance.
condition1 = rand(1,1);
condition2 = rand(1,1);
condition3 = rand(1,1);
dE        = abs(Energy(i,1) - Energy(i-1,1));
% Accept or reject new state.
if ((((PvMF1(i,1)>PvMF1(i-1,1)) || (condition1<exp(-abs(PvMF1(i,1)-PvMF1(i-1,1))))) && ((PvMF2(i,1)>PvMF2(i-1,1)) || (condition2<exp(-abs(PvMF2(i,1)-PvMF2(i-1,1)))))) && ((Energy(i,1)<Energy(i-1,1)) || (condition3<exp(-beta*dE))))
else
X(i,:) = X(i-1,:);
Y(i,:) = Y(i-1,:);
Z(i,:) = Z(i-1,:);
Energy(i,1)     = Energy(i-1,1);
Energy_WLC(i,1) = Energy_WLC(i-1,1);
PvMF1(i,1)       = PvMF1(i-1,1);
PvMF2(i,1)       = PvMF2(i-1,1);
end

else
Xsub                        = [X(i-1,:)',Y(i-1,:)',Z(i-1,:)'];
[Pr] = Quaternion3(thet(i,1),uu(i,:),Xsub);
X(i,:) = Pr(:,1)';
Y(i,:) = Pr(:,2)';
Z(i,:) = Pr(:,3)';
Energy_WLC(i,1) = Energy_WLC(i-1,1);
Energy(i,1)     = Energy(i-1,1);
Xmu1            = [X(i,1),Y(i,1),Z(i,1)];
Xmu2            = [X(i,len),Y(i,len),Z(i,len)];
PvMF1(i,1)       = (kappa *(mu*Xmu1'./sqrt(sum(Xmu1.^2))));
PvMF2(i,1)       = (kappa *(mu*Xmu2'./sqrt(sum(Xmu2.^2))));
%Set conditions for acceptance.
condition1 = rand(1,1);
condition2 = rand(1,1);
% Accept or reject new state.
if (((PvMF1(i,1)>PvMF1(i-1,1)) || (condition1<exp(-abs(PvMF1(i,1)-PvMF1(i-1,1))))) && ((PvMF2(i,1)>PvMF2(i-1,1)) || (condition2<exp(-abs(PvMF2(i,1)-PvMF2(i-1,1))))))
else
X(i,:) = X(i-1,:);
Y(i,:) = Y(i-1,:);
Z(i,:) = Z(i-1,:);
PvMF1(i,1)       = PvMF1(i-1,1);
PvMF2(i,1)       = PvMF2(i-1,1);
end
end
end
subsample       = 1:100:No;%linspace(50000,100000,5001);
X = X(subsample,:); 
Y = Y(subsample,:); 
Z = Z(subsample,:); 
Energy_WLC = Energy_WLC(subsample,1);

%Addressx1 = (['C:\WLCData\X' num2str(Chromval) '.mat']);
%Addressy1 = (['C:\WLCData\Y' num2str(Chromval) '.mat']);
%Addressz1 = (['C:\WLCData\Z' num2str(Chromval) '.mat']);
%AddressE1 = (['C:\WLCData\E' num2str(Chromval) '.mat']);
%save(Addressx1,'X')
%save(Addressy1,'Y')
%save(Addressz1,'Z')