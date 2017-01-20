function [net] = SingleTether(Parameters,X0x,Y0x,Z0x);

%SingleTether - Markov chain Monte Carlo algorithm for simulating chromosome trajectories when one end is fixedA single tupe of update is defined in the MCMC search: (i) An off 
%lattice rotation of one arm about a randomly selected bead. 
%
% usage: [X,Y,Z,Energy] = SingleTether(Parameters,x,y,z,VMF) 
%
% where
%
% Parameters        is a row vector containing all internal parameters [No,len,k,R,U,a], where
%	No              is the length of the MCMC chain
%	len             is the number of segments that make the chromosome
%   k               represents the persistence length of the chromosome
%   R               the radius of a spherical nuclear periphery (NP)
%   U               the confining potential of the NP
%   a               the length of each segment of the chromosome
%   x               start and end positions of chromosome (x values)
%   y               start and end positions of chromosome (y values)
%   z               start and end positions of chromosome (z values)
% VMF               Paremeters of the von Mises-Fisher distribution: a
%                   matrix [k1,xhat1,yhat1,zhat1;k1,xhat1,yhat1,zhat1] where
%   k1              represents the spread parameter of telomere 1 
%   [xhat,yhat,zhat]represents the mean vector
% X,Y,Z             represents an [No x len] matrix detailing the No sample configurations of the chromosome
% E                 is a column vector with the corresponding energies associated with each chromosome confuguration
%
% Example use can be found by running the script demBouquet.m
%
%    Copyright (C) 2010  Christopher Penfold

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

addpath('./Functions','./netlab')

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
En(1,1) = sum(diag(E1*E1',1));
%__________________________________________________________________________
%Calculation of energy due to confinement.
ll  = length(find(X(1,:).^2 + Y(1,:).^2 + Z(1,:).^2 <= (Radius+1)^2));
E_c = ll*0 + (len-ll)*U;
%__________________________________________________________________________
%Total starting energy.
Energy(1,1)     = -(k/a)*En + E_c + 3*sum(diff(X0x).^2 + diff(Y0x).^2 + diff(Z0x).^2)/(2*a^2);
Energy_WLC(1,1) = -(k/a)*En;
%__________________________________________________________________________
%Begin Metropolis-Hastings step.
for i = 2:No
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

E1 = [[rx(1,max(1,Bond1(i,1)-1):min(Bond1(i,1)+1,len-1))./Delta(1,max(1,Bond1(i,1)-1):min(Bond1(i,1)+1,len-1))]',...
      [ry(1,max(1,Bond1(i,1)-1):min(Bond1(i,1)+1,len-1))./Delta(1,max(1,Bond1(i,1)-1):min(Bond1(i,1)+1,len-1))]',...
      [rz(1,max(1,Bond1(i,1)-1):min(Bond1(i,1)+1,len-1))./Delta(1,max(1,Bond1(i,1)-1):min(Bond1(i,1)+1,len-1))]'];
%keyboard
E2 = [[rx1(1,max(1,Bond1(i,1)-1):min(Bond1(i,1)+1,len-1))./Delta1(1,max(1,Bond1(i,1)-1):min(Bond1(i,1)+1,len-1))]',...
      [ry1(1,max(1,Bond1(i,1)-1):min(Bond1(i,1)+1,len-1))./Delta1(1,max(1,Bond1(i,1)-1):min(Bond1(i,1)+1,len-1))]',...
      [rz1(1,max(1,Bond1(i,1)-1):min(Bond1(i,1)+1,len-1))./Delta1(1,max(1,Bond1(i,1)-1):min(Bond1(i,1)+1,len-1))]'];
En1(1,1) = sum(diag(E1*E1',1));
En2(1,1) = sum(diag(E2*E2',1));
%__________________________________________________________________________
%Calculate total energy of the new state.
ll   = length(find((X(i,:).^2 + Y(i,:).^2 + Z(i,:).^2 <= (Radius+1)^2)));
E_c  = (ll)*0 + ((len) -ll)*U;
Energy_WLC(i,1) = Energy_WLC(i-1,1) + (k/a)*En2 - (k/a)*En1;
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

net.X = X;
net.Y = Y;
net.Z = Z;
net.E = Energy;
net.Parameters = Parameters;

%subsample       = 1:50:100000;
%X = X(subsample,:); 
%Y = Y(subsample,:); 
%Z = Z(subsample,:); 
%Energy = Energy(subsample,1);
%Save to hard drive:
%Addressx1 = (['C:\WLCData\X' num2str(Chromval) '.mat']);
%Addressy1 = (['C:\WLCData\Y' num2str(Chromval) '.mat']);
%Addressz1 = (['C:\WLCData\Z' num2str(Chromval) '.mat']);
%AddressE1 = (['C:\WLCData\E' num2str(Chromval) '.mat']);
%save(Addressx1,'X')
%save(Addressy1,'Y')
%save(Addressz1,'Z')