function [X,Y,Z,Energy] = RablPrior(varargin);
%
%Description (version 3)
%
%[X,Y,Z,Energy] = Rabl3(Parameters,x,y,z,dx,dy,dz) generates an 
%arbtrary length Markov chain for chromosomes in a Rabl formation 
%with boundary conditions drawn from a von Mises-Fisher
%distribution. Two types of update are defined: (i) a crankshaft 
%rotation about an axis connecting two randomly determined bonds; 
%(ii) rotation of chromosomes via random axis. Example of use: RablScript3.

addpath('./Functions','./netlab')

switch nargin
case 0
	warning('No inputs')
case 1
	if length(varargin{1})==7
		disp('No initial conformation input. Generating Conformation now.')
		Parameters = varargin{1};
        [x0,y0,z0,theta,phi] = tMultivMF(0,0,1,1e-9,Parameters(1,4),1);
        [x1,y1,z1] = RandFlightCon([x0,y0,z0],[x0,y0,z0],Parameters(1,6),Parameters(1,2));
        dx = x1(1,Parameters(1,7)); dy = y1(1,Parameters(1,7)); dz = z1(1,Parameters(1,7));
        x1 = x1 - dx; y1 = y1 - dy; z1 = z1 - dz + Parameters(1,4);
		X0x = x1; Y0x = y1; Z0x = z1;
        VMF = [1e-9,0,0,1];
        disp('Initial Configuration Generated.')
	else
		warning('Not enough parameters specified in input')
	end

case 2
	if length(varargin{1})==6 && length(varargin{2})==4
		Parameters = varargin{1};
		VMF = varargin{2}
        [x0,y0,z0,theta,phi] = tMultivMF(0,0,1,1e-9,Parameters(1,4),1);
        [x1,y1,z1] = RandFlightCon([x0,y0,z0],[x0,y0,z0],Parameters(1,6),Parameters(1,2));
        dx = x1(1,Parameters(1,7)); dy = y1(1,Parameters(1,7)); dz = z1(1,Parameters(1,7));
        x1 = x1 - dx; y1 = y1 - dy; z1 = z1 - dz + Parameters(1,4);
		X0x = x1; Y0x = y1; Z0x = z1;
        %disp('No initial conformation input. Generating Conformation now.')
	else
		warning('Not enough parameters specified in one or more inputs')
	end

case 3
		warning('Wrong number of inputs specified')
case 4
		warning('Wrong number of inputs specified')
case 5
				Parameters = varargin{1};
				X0x = varargin{2}; 
				Y0x = varargin{3}; 
				Z0x = varargin{4};
				VMF = varargin{5};
end

%__________________________________________________________________________
%Physical parameters
len    = Parameters(1,2); %Length of chromosome (No. of Khun segments).
k      = Parameters(1,3); %Rigidity constant.
Radius = Parameters(1,4); %Radius of NP.
a      = Parameters(1,6); %Khun segment length.
U      = Parameters(1,5); %Confinement potential.
Cen    = Parameters(1,7); %Centromere location.
beta   = 1;               %Inverse temperature.
%__________________________________________________________________________
%Unphysical parameters
No          = Parameters(1,1);                      %No. steps in Markov chain.
mu          = VMF(1,2:4);                           %Mean vector.
Xmu1        = [X0x(1,Cen),Y0x(1,Cen),Z0x(1,Cen)];   %Position of walker on NP.
kappa       = VMF(1,1);                             %Spread parameter.
%__________________________________________________________________________
%Empty arrays for storing data points during MC.
thet        = rand(No,1)*2*pi;          %Angle of rotation ~ U(0,2pi)
thet1        = rand(No,1)*2*pi;          %Angle of rotation ~ U(0,2pi)
[xxx,yxx,zxx]= Points_on_Sphere(1,No);  %Axis of rotation.
uu          = [xxx',yxx',zxx'];         %Axis of rotation.
[xxx,yxx,zxx]= Points_on_Sphere(1,No);  %Axis of rotation.
uu1          = [xxx',yxx',zxx'];         %Axis of rotation.
Bond1       =randint2(No,1,Cen-1)+2;      %Bond of rotation for telomere arm 1.
Bond2       =randint2(No,1,len-Cen)+Cen;%Bond of rotation for telomere arm 2.
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
X(1,:)  = X0x;
Y(1,:)  = Y0x;
Z(1,:)  = Z0x;
%__________________________________________________________________________
%Calculate bending energies.
Delta   = sqrt(diff(X0x).^2 + diff(Y0x).^2 + diff(Z0x).^2);
E1      = [(diff(X0x)./Delta)',(diff(Y0x)./Delta)',(diff(Z0x)./Delta)'];
En(1,1) = sum(diag(E1*E1',1));
%__________________________________________________________________________
%Calculation of energy due to confinement.
ll  = length(find(X(1,:).^2 + Y(1,:).^2 + Z(1,:).^2 <= (Radius+1)^2));
E_c = ll*0 + (len-ll)*U;
%__________________________________________________________________________
%Total starting energy.
PvMF1(1,1)      = (kappa*(mu*Xmu1'./sqrt(sum(Xmu1.^2))));
Energy(1,1)     = -(k/a)*En + E_c + 3*sum(diff(X0x).^2 + diff(Y0x).^2 + diff(Z0x).^2)/(2*a^2) - PvMF1;
Energy_WLC(1,1) = -(k/a)*En;
%__________________________________________________________________________
%Begin Metropolis-Hastings step.
for i = 2:No
stepchoice = stepc(1,i);  
if stepchoice < (0.5)
dx      = X(i-1,Bond1(i,1));
dy      = Y(i-1,Bond1(i,1));
dz      = Z(i-1,Bond1(i,1));

X(i,:) = X(i-1,:);
Y(i,:) = Y(i-1,:);
Z(i,:) = Z(i-1,:);

Xsub    = [X(i-1,1:Bond1(i,1))'-dx,Y(i-1,1:Bond1(i,1))'-dy,Z(i-1,1:Bond1(i,1))'-dz];
[Pr]    = Quaternion3(thet(i,1),uu(i,:),Xsub);
X(i,1:Bond1(i,1)) = Pr(:,1)+dx;
Y(i,1:Bond1(i,1)) = Pr(:,2)+dy;
Z(i,1:Bond1(i,1)) = Pr(:,3)+dz;

dx      = X(i,Bond2(i,1));
dy      = Y(i,Bond2(i,1));
dz      = Z(i,Bond2(i,1));

Xsub    = [X(i-1,Bond2(i,1):len)'-dx,Y(i-1,Bond2(i,1):len)'-dy,Z(i-1,Bond2(i,1):len)'-dz];
[Pr]    = Quaternion3(thet1(i,1),uu1(i,:),Xsub);
X(i,Bond2(i,1):len) = Pr(:,1)+dx;
Y(i,Bond2(i,1):len) = Pr(:,2)+dy;
Z(i,Bond2(i,1):len) = Pr(:,3)+dz;
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
   E1 = [[rx(1,:)./Delta(1,:)]', [ry(1,:)./Delta(1,:)]', [rz(1,:)./Delta(1,:)]'];
   E2 = [[rx1(1,:)./Delta1(1,:)]', [ry1(1,:)./Delta1(1,:)]', [rz1(1,:)./Delta1(1,:)]'];
  
En1(1,1) = sum(diag(E1*E1',1));
En2(1,1) = sum(diag(E2*E2',1));
%__________________________________________________________________________
%Calculate total energy of the new state.
ll   = length(find((X(i,:).^2 + Y(i,:).^2 + Z(i,:).^2 <= (Radius+1)^2)));
PvMF1(i,1) = PvMF1(i-1,1);
E_c  = (ll)*0 + ((len) -ll)*U;
Energy_WLC(i,1) = Energy_WLC(i-1,1) + (k/a)*En2 - (k/a)*En1;
Energy(i,1) = Energy_WLC(i,1) + E_c + 3*sum(rx.^2 + ry.^2 +rz.^2)/(2*a^2) - PvMF1(i,1) + PvMF1(i-1,1);
%Set conditions for acceptance.
condition = rand(1,1);
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

else

Xsub = [X(i-1,:)',Y(i-1,:)',Z(i-1,:)'];
[Pr] = Quaternion3(thet(i,1),uu(i,:),Xsub);
X(i,:) = Pr(:,1)';
Y(i,:) = Pr(:,2)';
Z(i,:) = Pr(:,3)';

%deltaxyz = sqrt(X(i,Cen)^2 + Y(i,Cen)^2 + Z(i,Cen)^2);
%X(i,:) = X(i,:) - X(i,Cen) + Radius*(X(i,Cen)/deltaxyz);
%Y(i,:) = Y(i,:) - Y(i,Cen) + Radius*(Y(i,Cen)/deltaxyz);
%Z(i,:) = Z(i,:) - Z(i,Cen) + Radius*(Z(i,Cen)/deltaxyz);

rx      = diff(X(i,:));
ry      = diff(Y(i,:));
rz      = diff(Z(i,:));
%__________________________________________________________________________
ll   = length(find((X(i,:).^2 + Y(i,:).^2 + Z(i,:).^2 <= (Radius+1)^2)));
Xmu1            = [X(i,Cen),Y(i,Cen),Z(i,Cen)];
PvMF1(i,1)       = (kappa *(mu*Xmu1'./sqrt(sum(Xmu1.^2))));
E_c  = (ll)*0 + ((len) -ll)*U;
Energy_WLC(i,1) = Energy_WLC(i-1,1);
Energy(i,1) = Energy_WLC(i,1) + E_c + 3*sum(rx.^2 + ry.^2 +rz.^2)/(2*a^2)  - PvMF1(i,1) + PvMF1(i-1,1);
condition = rand(1,1);
%Set conditions for acceptance.
dE        = abs(Energy(i,1) - Energy(i-1,1));
%Accept or reject new state.
if ((Energy(i,1)<Energy(i-1,1)) || (condition1(i,1)<exp(-beta*dE)))
else
X(i,:) = X(i-1,:);
Y(i,:) = Y(i-1,:);
Z(i,:) = Z(i-1,:);
Energy_WLC(i,1) = Energy_WLC(i-1,1);
Energy(i,1) = Energy(i-1,1);
PvMF1(i,1) = PvMF1(i-1,1);
end
end
end
subsample       = 1:5:No;%linspace(50000,100000,5001);
X = X(subsample,:); 
Y = Y(subsample,:); 
Z = Z(subsample,:); 
Energy_WLC = Energy_WLC(subsample,1);
%Save to hard drive:
%Addressx1 = (['C:\WLCData\X' num2str(Chromval) '.mat']);
%Addressy1 = (['C:\WLCData\Y' num2str(Chromval) '.mat']);
%Addressz1 = (['C:\WLCData\Z' num2str(Chromval) '.mat']);
%AddressE1 = (['C:\WLCData\E' num2str(Chromval) '.mat']);
%save(Addressx1,'X')
%save(Addressy1,'Y')
%save(Addressz1,'Z')