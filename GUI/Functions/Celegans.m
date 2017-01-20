function [net] = Celegans(varargin);

%C.elegans - Markov chain Monte Carlo algorithm for simulating chromosome trajectories when one end is fixed to the NP but the other is unfixed as observed in meiotic cells in C.elegans 
%
% usage: [X,Y,Z,Energy] = Celegans(Parameters,x,y,z,VMF) 
%
% where
%
% Parameters        is a row vector containing all internal parameters [No,len,k,R,U,a,Cen], where
%	No              is the length of the MCMC chain
%	len             is the number of segments that make the chromosome
%   k               represents the persistence length of the chromosome
%   R               the radius of a spherical nuclear periphery (NP)
%   U               the confining potential of the NP
%   a               the length of each segment of the chromosome
%   Cen             the position (bead no.) of the centromere 
%   x               start and end positions of chromosome (x values)
%   y               start and end positions of chromosome (y values)
%   z               start and end positions of chromosome (z values)
% VMF               Paremeters of the von Mises-Fisher distribution: a
%                   vector [k1,xhat1,yhat1,zhat1] where,
%   k1              represents the spread parameter of telomere 1 
%   [xhat,yhat,zhat]represents the mean vector
% X,Y,Z             represents an [No x len] matrix detailing the No sample configurations of the chromosome
% E                 is a column vector with the corresponding energies
%                   associated with each chromosome confuguration (may be
%                   used to numerically estimate partition funciton Z)
%
% Example use can be found by running the script demCelegans.m
%
%    Copyright (C) 2010  Christopher Penfold
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

%PEB added this in case user cancels
net = [];



switch nargin
case 0
warning('No inputs')
case 1
if length(varargin{1})==7
disp('No initial conformation input. Generating Conformation now.')
Parameters = varargin{1};
[x1,y1,z1]              = RandFlightCon([0,0,0],[1,1,1],Parameters(1,6),Parameters(1,2));
X0x = x1; Y0x = y1; Z0x = z1;
dx = x1(1,1);
dy = y1(1,1);
dz = z1(1,1);
x1 = x1 - dx;
y1 = y1 - dy;
z1 = z1 - dz + Parameters(1,4);
VMF = [1e-10,0,0,1];
disp('Initial Configuration Generated.')
else
warning('Not enough parameters specified in input')
end
case 2
warning('Wrong number of inputs specified')
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
mu          = VMF(1,2:4);         % Mean cluster direciton.
kappa       = VMF(1,1);           % Cluster strength.

%__________________________________________________________________________
%Empty arrays for storing data points during MC.
thet        = rand(No,1)*2*pi;          %Angle of rotation ~ U(0,2pi)
[xxx,yxx,zxx]= Points_on_Sphere(1,No);  %Axis of rotation.
Xmu1        = [X0x(1,1),Y0x(1,1),Z0x(1,1)];   %Position of walker on NP.
uu          = [xxx',yxx',zxx'];         %Axis of rotation.
Bond1       = randint(No,1,len-1)+1;    %Bond of rotation for telomere arm 1.
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
X(1,:)  = X0x;
Y(1,:)  = Y0x;
Z(1,:)  = Z0x;
%Contribution from microtubules
PvMF1(1,1)      = (kappa*(mu*Xmu1'./sqrt(sum(Xmu1.^2))));			   
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
Energy(1,1)     = -(k/a)*En + E_c + 3*sum(diff(X0x).^2 + diff(Y0x).^2 + diff(Z0x).^2)/(2*a^2) - PvMF1(1,1);
Energy_WLC(1,1) = -(k/a)*En;
%__________________________________________________________________________
%Begin Metropolis-Hastings step.
for i = 2:No

    stepchoice = stepc(1,i);
				   
if stepchoice<0.5				   
				   
    dx      = X(i-1,Bond1(i,1));
    dy      = Y(i-1,Bond1(i,1));
    dz      = Z(i-1,Bond1(i,1));
    X(i,:)  = X(i-1,:);
    Y(i,:)  = Y(i-1,:);
    Z(i,:)  = Z(i-1,:);
    Xsub    = [X(i-1,Bond1(i,1):len)'-dx,Y(i-1,Bond1(i,1):len)'-dy,Z(i-1,Bond1(i,1):len)'-dz];
    [Pr]    = Quaternion3(thet(i,1),uu(i,:),Xsub);
    X(i,Bond1(i,1):len) = Pr(:,1)+dx;
    Y(i,Bond1(i,1):len) = Pr(:,2)+dy;
    Z(i,Bond1(i,1):len) = Pr(:,3)+dz;
%__________________________________________________________________________
    rx      = diff(X(i,:));
    ry      = diff(Y(i,:));
    rz      = diff(Z(i,:));
    rx1     = diff(X(i-1,:));
    ry1     = diff(Y(i-1,:));
    rz1     = diff(Z(i-1,:));
    Delta   = sqrt(rx.^2 + ry.^2 + rz.^2);
    Delta1  = sqrt(rx1.^2 + ry1.^2 + rz1.^2);
%__________________________________________________________________________
%Calculate bond energies.
    E1 = [[rx(1,max(1,Bond1(i,1)-1):min(Bond1(i,1)+1,len-1))./Delta(1,max(1,Bond1(i,1)-1):min(Bond1(i,1)+1,len-1))]',...
      [ry(1,max(1,Bond1(i,1)-1):min(Bond1(i,1)+1,len-1))./Delta(1,max(1,Bond1(i,1)-1):min(Bond1(i,1)+1,len-1))]',...
      [rz(1,max(1,Bond1(i,1)-1):min(Bond1(i,1)+1,len-1))./Delta(1,max(1,Bond1(i,1)-1):min(Bond1(i,1)+1,len-1))]'];
    E2 = [[rx1(1,max(1,Bond1(i,1)-1):min(Bond1(i,1)+1,len-1))./Delta1(1,max(1,Bond1(i,1)-1):min(Bond1(i,1)+1,len-1))]',...
      [ry1(1,max(1,Bond1(i,1)-1):min(Bond1(i,1)+1,len-1))./Delta1(1,max(1,Bond1(i,1)-1):min(Bond1(i,1)+1,len-1))]',...
      [rz1(1,max(1,Bond1(i,1)-1):min(Bond1(i,1)+1,len-1))./Delta1(1,max(1,Bond1(i,1)-1):min(Bond1(i,1)+1,len-1))]'];
    En1(1,1) = sum(diag(E1*E1',1));
    En2(1,1) = sum(diag(E2*E2',1));
%__________________________________________________________________________
%Calculate total energy of the new state.
    PvMF1(i,1) = PvMF1(i-1,1);
    ll   = length(find((X(i,:).^2 + Y(i,:).^2 + Z(i,:).^2 <= (Radius+1)^2)));
    E_c  = (ll)*0 + ((len) -ll)*U;
    Energy_WLC(i,1) = Energy_WLC(i-1,1) + (k/a)*En2 - (k/a)*En1;
    Energy(i,1)     = Energy_WLC(i,1) + E_c + 3*sum(rx.^2 + ry.^2 +rz.^2)/(2*a^2) - PvMF1(i,1);
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

    else
				   
				   %Rotation of entire chromosome				   
				   				   
    Xsub = [X(i-1,:)',Y(i-1,:)',Z(i-1,:)'];
	[Pr] = Quaternion3(thet(i,1),uu(i,:),Xsub);				   
    X(i,:) = Pr(:,1)';
	Y(i,:) = Pr(:,2)';
	Z(i,:) = Pr(:,3)';						   
	rx     = diff(X(i,:));
	ry     = diff(Y(i,:));
	rz     = diff(Z(i,:));
	Energy_WLC(i,1) = Energy_WLC(i-1,1);
	Xmu1            = [X(i,1),Y(i,1),Z(i,1)];
	ll   = length(find((X(i,:).^2 + Y(i,:).^2 + Z(i,:).^2 <= (Radius+1)^2)));
	E_c  = (ll)*0 + ((len) -ll)*U;
	PvMF1(i,1)      = (kappa*(mu*Xmu1'./sqrt(sum(Xmu1.^2))));				   
	Energy(i,1) = Energy_WLC(i,1) + E_c + 3*sum(rx.^2 + ry.^2 +rz.^2)/(2*a^2) - PvMF1(i,1);		
    dE        = abs(Energy(i,1) - Energy(i-1,1));
		if ((Energy(i,1)<Energy(i-1,1)) || (condition1(i,1)<exp(-beta*dE)))
		else
		X(i,:)  = X(i-1,:);
		Y(i,:)  = Y(i-1,:);
		Z(i,:)  = Z(i-1,:);
		Energy(i,1)     = Energy(i-1,1);
		Energy_WLC(i,1) = Energy_WLC(i-1,1);
        		PvMF(i,1) = PvMF(i-1,1);
        end			 												 
end				  				   
    if ~mod(i, floor(No/20))
        %update ~20 times in total
        if UpdateProgress(['Running Analysis    ' num2str(i) '/' num2str(No)])
            return;
        end
    %
    %end

end %end loop
%printStr('Done');

subsample       = double(int64(No/2)):10:No;
X = X(subsample,:); 
Y = Y(subsample,:); 
Z = Z(subsample,:); 
Energy = Energy(:,1);

net.X = X;
net.Y = Y;
net.Z = Z;
net.E = Energy;
net.Parameters = Parameters;

end
%-----------------------------------------------------------------------