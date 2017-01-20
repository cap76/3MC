function [X,Y,Z,Energy] = FreeWLCTV(varargin);
% Generates an arbtrary length Markov chain for chromosomes in a Rabl
% formation with boundary conditions drawn from a von Mises-Fisher
% distribution. Two types of update are defined: (i) a crankshaft rotation
% about an axis connecting two randomly determined bonds; 
% (ii) rotation of chromosomes via random axis.

addpath('./Functions','./netlab')

switch nargin
case 0
	warning('No inputs')
case 1
	if length(varargin{1})==6
		disp('No initial conformation input. Generating Conformation now.')
		Parameters = varargin{1};
		[x1,y1,z1]              = RandFlightCon([0,0,0],[1,1,1],Parameters(1,6),Parameters(1,2));
		X0x = x1; Y0x = y1; Z0x = z1;
        disp('Initial Configuration Generated.')
	else
		warning('Not enough parameters specified in input')
	end

case 2
        warning('Wrong number of inputs specified')
case 3
		warning('Wrong number of inputs specified')
case 4
				Parameters = varargin{1};
				X0x = varargin{2}; 
				Y0x = varargin{3}; 
				Z0x = varargin{4};
end

%__________________________________________________________________________
%Physical parameters
len    = Parameters(1,2); %Length of chromosome (No. of Khun segments).
Radius = Parameters(1,3); %Radius of NP.
a      = Parameters(1,5); %Khun segment length.
U      = Parameters(1,4); %Confinement potential.
k      = Parameters(1,6:end)'; %Rigidity constant.
beta   = 1;               %Inverse temperature.
%__________________________________________________________________________
%Unphysical parameters
No          = Parameters(1,1);                      %No. steps in Markov chain.
%__________________________________________________________________________
%Empty arrays for storing data points during MC.
thet        = rand(No,1)*2*pi;          %Angle of rotation ~ U(0,2pi)
[xxx,yxx,zxx]= Points_on_Sphere(1,No);  %Axis of rotation.
uu          = [xxx',yxx',zxx'];         %Axis of rotation.
Bond1       =randint(No,1,len-1)+1;     %Bond of rotation for telomere arm 1.
step        = randn(3,No)*Radius/10;    %Random walk size for the polymer.
stepc       = rand(1,No);               %Choice of update.
X           = zeros(No,len);            %Array for X coordinates.
Y           = zeros(No,len);            %Array for Y coordinates.
Z           = zeros(No,len);            %Array for Z coordinate.
En          = zeros(1,1);               %Array for total free energy.
Energy      = zeros(No,1);              %Array for total free energy.
Energy_WLC  = zeros(No,1);              %Array for bending energy.
PvMF1       = zeros(No,1);              %Array for von Mises Fisher likelihood.
condition1  = rand(No,1);               %Acceptance criterior for energy.
condition2  = rand(No,1);               %Acceptance criterior for vMF.
%__________________________________________________________________________
%Set first state in the Markov chain equal to initial conditions.
X(1,:)  = X0x - X0x(1,1);
Y(1,:)  = Y0x - Y0x(1,1);
Z(1,:)  = Z0x - Z0x(1,1);
%__________________________________________________________________________
%Calculate bending energies.
Delta   = sqrt(diff(X0x).^2 + diff(Y0x).^2 + diff(Z0x).^2);
E1      = [(diff(X0x)./Delta)',(diff(Y0x)./Delta)',(diff(Z0x)./Delta)'];
%keyboard
En(1,1) = sum(k.*diag(E1*E1',1));
%__________________________________________________________________________
%Calculation of energy due to confinement.
ll  = length(find(X(1,:).^2 + Y(1,:).^2 + Z(1,:).^2 <= (Radius+1)^2));
E_c = ll*0 + (len-ll)*U;
%__________________________________________________________________________
%Total starting energy.
Energy(1,1)     = -(1/a)*En + E_c + 3*sum(diff(X0x).^2 + diff(Y0x).^2 + diff(Z0x).^2)/(2*a^2);
Energy_WLC(1,1) = -(1/a)*En;
%__________________________________________________________________________
%Begin Metropolis-Hastings step.
for i = 2:No
stepchoice = stepc(1,i);  
if stepchoice < (1/3)
dx      = X(i-1,Bond1(i,1));
dy      = Y(i-1,Bond1(i,1));
dz      = Z(i-1,Bond1(i,1));

X(i,:) = X(i-1,:);
Y(i,:) = Y(i-1,:);
Z(i,:) = Z(i-1,:);

Xsub    = [X(i-1,Bond1(i,1):len)'-dx,Y(i-1,Bond1(i,1):len)'-dy,Z(i-1,Bond1(i,1):len)'-dz];
[Pr]    = Quaternion3(thet(i,1),uu(i,:),Xsub);
X(i,Bond1(i,1):len) = Pr(:,1)+dx;
Y(i,Bond1(i,1):len) = Pr(:,2)+dy;
Z(i,Bond1(i,1):len) = Pr(:,3)+dz;
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
Delta   = sqrt((rx).^2 + (ry).^2 + (rz).^2);
E1      = [((rx)./Delta)',((ry)./Delta)',((rz)./Delta)'];
%keyboard
En(1,1) = sum(k.*diag(E1*E1',1));
%__________________________________________________________________________
%Calculate total energy of the new state.
ll   = length(find((X(i,:).^2 + Y(i,:).^2 + Z(i,:).^2 <= (Radius+1)^2)));
E_c  = (ll)*0 + ((len) -ll)*U;

Energy(i,1)     = -(1/a)*En + E_c + 3*sum((rx).^2 + (ry).^2 + (rz).^2)/(2*a^2);
Energy_WLC(i,1) = -(1/a)*En;
PvMF1(i,1) = PvMF1(i-1,1);
%Set conditions for acceptance.
dE        = abs(Energy(i,1) - Energy(i-1,1));
%Accept or reject new state.
if ((Energy(i,1)<Energy(i-1,1)) || (condition1(i,1)<exp(-beta*dE)))
else
X(i,:)  = X(i-1,:);
Y(i,:)  = Y(i-1,:);
Z(i,:)  = Z(i-1,:);
Energy(i,1)     = Energy(i-1,1);
Energy_WLC(i,1) = Energy_WLC(i-1,1);
end

elseif stepchoice > (2/3)

Xsub = [X(i-1,:)',Y(i-1,:)',Z(i-1,:)'];
[Pr] = Quaternion3(thet(i,1),uu(i,:),Xsub);
X(i,:) = Pr(:,1)';
Y(i,:) = Pr(:,2)';
Z(i,:) = Pr(:,3)';

rx      = diff(X(i,:));
ry      = diff(Y(i,:));
rz      = diff(Z(i,:));
Energy_WLC(i,1) = Energy_WLC(i-1,1);
Energy(i,1) = Energy(i-1,1);
%Set conditions for acceptance.
%Accept or reject new state.

else

X(i,:) = X(i-1,:) + (1/3)*step(1,i);
Y(i,:) = Y(i-1,:) + (1/3)*step(2,i);
Z(i,:) = Z(i-1,:) + (1/3)*step(3,i);
rx = diff(X(i,:));
ry = diff(Y(i,:));
rz = diff(Z(i,:));
   
%__________________________________________________________________________
%Calculate total energy of the new state.
ll   = length(find((X(i,:).^2 + Y(i,:).^2 + Z(i,:).^2 <= (Radius+1)^2)));
E_c  = (ll)*0 + ((len) -ll)*U;
Energy_WLC(i,1) = Energy_WLC(i-1,1);
Energy(i,1) = Energy_WLC(i,1) + E_c + 3*sum(rx.^2 + ry.^2 +rz.^2)/(2*a^2);
PvMF1(i,1) = PvMF1(i-1,1);
%Set conditions for acceptance.
dE        = abs(Energy(i,1) - Energy(i-1,1));
%Accept or reject new state.
if ((Energy(i,1)<Energy(i-1,1)) || (condition1(i,1)<exp(-beta*dE)))
else
X(i,:)  = X(i-1,:);
Y(i,:)  = Y(i-1,:);
Z(i,:)  = Z(i-1,:);
Energy(i,1)     = Energy(i-1,1);
Energy_WLC(i,1) = Energy_WLC(i-1,1);
end
end
end
subsample       = 1:100:No;%linspace(50000,100000,5001);
X = X(subsample,:); 
Y = Y(subsample,:); 
Z = Z(subsample,:); 
Energy = Energy(subsample,1);
%Save to hard drive:
%Addressx1 = (['C:\WLCData\X' num2str(Chromval) '.mat']);
%Addressy1 = (['C:\WLCData\Y' num2str(Chromval) '.mat']);
%Addressz1 = (['C:\WLCData\Z' num2str(Chromval) '.mat']);
%AddressE1 = (['C:\WLCData\E' num2str(Chromval) '.mat']);
%save(Addressx1,'X')
%save(Addressy1,'Y')
%save(Addressz1,'Z')