function [net] = Rabl2Bouquet(varargin)%Parameters,X0x,Y0x,Z0x,VMF);

%Rabl2BouquetPrior - Markov chain Monte Carlo algorithm for simulating chromosome trajectories when telomeres and centromeres are attached to the nuclear periphery and conditionally distributed according to a von Mises-Fisher distribution. 
%Three types of update are defined in the MCMC search: (i) A crankshaft
%rotation about an axis connecting two randomly determined bonds; (ii) An arm rotation of either arm; and (iii) Rotation of chromosomes via random axis. 
%
% usage: [X,Y,Z,Energy] = Rabl2BouquetPrior(Parameters,x,y,z,VMF) 
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
%   x (optional)    starting x cooridnates of chromosome
%   y (optional)    starting y cooridnates of chromosome
%   z (optional)    starting z cooridnates of chromosome
% VMF               Paremeters of the von Mises-Fisher distribution: a matrix [k1,xhat1,yhat1,zhat1;k1,xhat1,yhat1,zhat1;k3,xhat3,yhat3,zhat3] where
%   k1               represents the spread parameter of telomere 1, 2 and centromere respecitively 
%   [xhat,yhat,zhat]represents the mean vector
% X,Y,Z             represents an [No x len] matrix detailing the No sample configurations of the chromosome
% E                 is a column vector with the corresponding energies associated with each chromosome confuguration
%
% Example use can be found by running the script demRabl2Bouquet.m
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
net = [];


switch nargin
case 0
	warning('No inputs 1')
case 1
	if length(varargin{1})==7
		disp('No initial conformation input. Generating Conformation now.')
		Parameters = varargin{1};        
        [x0,y0,z0,theta,phi]    = tMultivMF(0,0,1,Parameters(1,3),Parameters(1,4),1);
        [x01,y01,z01,theta,phi] = tMultivMF(0,0,1,Parameters(1,3),Parameters(1,4),1);
        [x02,y02,z02,theta,phi] = tMultivMF(0,0,1,Parameters(1,3),Parameters(1,4),1);
        [x1,y1,z1] = RandFlightCon([x0,y0,z0],[x01,y01,z01],Parameters(1,6),Parameters(1,2)-Parameters(1,7));
        [x2,y2,z2] = RandFlightCon([x01,y01,z01],[x02,y02,z02],Parameters(1,6),Parameters(1,7)+1);
        x = [x1,x2(1,2:length(x2))];
        y = [y1,y2(1,2:length(x2))];
        z = [z1,z2(1,2:length(x2))];        
		[x0,y0,z0,theta,phi]    = tMultivMF(0,0,1,1e-9,Parameters(1,4),1);
		[x1,y1,z1]              = RandFlightCon([x0,y0,z0],[x0,y0,z0],Parameters(1,6),Parameters(1,2));
		X0x = x1; Y0x = y1; Z0x = z1;
        VMF = [1e-9,0,0,1;1e-9,0,0,1;1e-9,0,0,1];
        disp('Initial Configuration Generated.')
	else
		warning('Not enough parameters specified in input 2')
	end

case 2
	if length(varargin{1})==7 && size(varargin{2})==[3,2]
		Parameters = varargin{1};
		VMF = varargin{2};
        [x0,y0,z0,theta,phi]    = tMultivMF(0,0,1,Parameters(1,3),Parameters(1,4),1);
        [x01,y01,z01,theta,phi] = tMultivMF(0,0,1,Parameters(1,3),Parameters(1,4),1);
        [x02,y02,z02,theta,phi] = tMultivMF(0,0,1,Parameters(1,3),Parameters(1,4),1);
        [x1,y1,z1] = RandFlightCon([x0,y0,z0],[x01,y01,z01],Parameters(1,6),Parameters(1,2)-Parameters(1,7));
        [x2,y2,z2] = RandFlightCon([x01,y01,z01],[x02,y02,z02],Parameters(1,6),Parameters(1,7)+1);
        x = [x1,x2(1,2:length(x2))];
        y = [y1,y2(1,2:length(x2))];
        z = [z1,z2(1,2:length(x2))];        
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
k      = Parameters(1,3); % Rigidity constant.
Radius = Parameters(1,4); % Radius of NP.
a      = Parameters(1,6); % Khun segment length.
U      = Parameters(1,5); % Confinement potential.
Cen    = Parameters(1,7); % Centromere location
beta   = 1;               % Inverse temperature.
%__________________________________________________________________________
%Unphysical parameters
No           = Parameters(1,1);    % No. steps in Markov chain.
mu1          = VMF(1,2:4);
mu2          = VMF(2,2:4);
mu3          = VMF(3,2:4);
Xmu1         = [X0x(1,1),Y0x(1,1),Z0x(1,1)];
Xmu2         = [X0x(1,len),Y0x(1,len),Z0x(1,len)];
Xmu3         = [X0x(1,Cen),Y0x(1,Cen),Z0x(1,Cen)];
kappa1       = VMF(1,1);
kappa2       = VMF(2,1);
kappa3       = VMF(3,1);
%__________________________________________________________________________
%Empty arrays for storing data points during MC.
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
PvMF3           = zeros(No,1);
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
PvMF1(1,1)       = (kappa1*(mu1*Xmu1'./sqrt(sum(Xmu1.^2))));
PvMF2(1,1)       = (kappa2*(mu2*Xmu2'./sqrt(sum(Xmu2.^2))));
PvMF3(1,1)       = (kappa3*(mu3*Xmu3'./sqrt(sum(Xmu3.^2))));
Energy(1,1)     = -(k/a)*En + E_c + 3*sum(diff(X0x).^2 + diff(Y0x).^2 + diff(Z0x).^2)/(2*a^2) - (PvMF1(1,1)+PvMF2(1,1)+PvMF3(1,1));
Energy_WLC(1,1) = -(k/a)*En;
%__________________________________________________________________________
%Begin Metropolis-Hastings step.
for i = 2:No
stepchoice = stepc(1,i);  



if stepchoice < (1/3)

Xsub = flipdim(X(i-1,1:Cen),2);
Ysub = flipdim(Y(i-1,1:Cen),2);
Zsub = flipdim(Z(i-1,1:Cen),2);
%
[Xs,Ys,Zs] = ArmRotate(Xsub,Ysub,Zsub);
%
X(i,:) = X(i-1,:);
Y(i,:) = Y(i-1,:);
Z(i,:) = Z(i-1,:);
%
X(i,1:Cen) = flipdim(Xs,2);
Y(i,1:Cen) = flipdim(Ys,2);
Z(i,1:Cen) = flipdim(Zs,2);
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
   E1 = [[rx./Delta]', ...
         [ry./Delta]', ...
         [rz./Delta]'];
     
   E2 = [[rx1./Delta1]', ...
         [ry1./Delta1]', ...
         [rz1./Delta1]'];
%
En1(1,1) = sum(diag(E1*E1',1));
En2(1,1) = sum(diag(E2*E2',1));
%__________________________________________________________________________
% Calculate total energy of the new state.
Xmu1       =[X(i,1),Y(i,1),Z(i,1)];
Xmu2       =[X(i,len),Y(i,len),Z(i,len)];
Xmu3       =[X(i,Cen),Y(i,Cen),Z(i,Cen)];
PvMF1(i,1) = (kappa1 *(mu1*Xmu1'./sqrt(sum(Xmu1.^2))));
PvMF2(i,1) = (kappa2 *(mu2*Xmu2'./sqrt(sum(Xmu2.^2))));
PvMF3(i,1) = (kappa3 *(mu3*Xmu3'./sqrt(sum(Xmu3.^2))));
ll         = length(find((X(i,:).^2 + Y(i,:).^2 + Z(i,:).^2 <= (Radius+1)^2)));
E_c        = (ll)*0 + ((len) -ll)*U;
Energy_WLC(i,1) = Energy_WLC(i-1,1) + (k/a)*En2 - (k/a)*En1;
Energy(i,1)= Energy_WLC(i,1) + E_c + 3*sum(rx.^2 + ry.^2 +rz.^2)/(2*a^2) - (PvMF1(i,1)+PvMF2(i,1)+PvMF3(i,1));

%Set conditions for acceptance.
condition1 = rand(1,1);
condition2 = rand(1,1);
condition3 = rand(1,1);
condition4 = rand(1,1);
dE         = abs(Energy(i,1) - Energy(i-1,1));
%__________________________________________________________________________
%Accept or reject new state.
if ((Energy(i,1)<Energy(i-1,1)) || (condition1<exp(-beta*dE)))
else
X(i,:)  = X(i-1,:);
Y(i,:)  = Y(i-1,:);
Z(i,:)  = Z(i-1,:);
Energy(i,1)     = Energy(i-1,1);
Energy_WLC(i,1) = Energy_WLC(i-1,1);
PvMF1(i,1) = PvMF1(i-1,1);
PvMF2(i,1) = PvMF2(i-1,1);
PvMF3(i,1) = PvMF3(i-1,1);
end
%X(i,:) = X(i-1,:);
%Y(i,:) = Y(i-1,:);
%Z(i,:) = Z(i-1,:);
%Energy(i,:) = Energy(i-1,:);
%Energy_WLC(i,:) = Energy_WLC(i-1,:);

Xsub = X(i,1:Cen);
Ysub = Y(i,1:Cen);
Zsub = Z(i,1:Cen);
[Xs,Ys,Zs] = CrankshaftRotate(Xsub,Ysub,Zsub);


Xsub2 = X(i,:);
Ysub2 = Y(i,:);
Zsub2 = Z(i,:);
Xsub2(1:Cen) = Xs;
Ysub2(1:Cen) = Ys;
Zsub2(1:Cen) = Zs;

rx = diff(Xsub2);
ry = diff(Ysub2);
rz = diff(Zsub2);
rx1 = diff(X(i,:));
ry1 = diff(Y(i,:));
rz1 = diff(Z(i,:));
Delta = sqrt(rx.^2 + ry.^2 + rz.^2);
Delta1 = sqrt(rx1.^2 + ry1.^2 + rz1.^2);
%__________________________________________________________________________
   E1 = [[rx./Delta]', ...
         [ry./Delta]', ...
         [rz./Delta]'];
     
   E2 = [[rx1./Delta1]', ...
         [ry1./Delta1]', ...
         [rz1./Delta1]'];

En1(1,1) = sum(diag(E1*E1',1));
En2(1,1) = sum(diag(E2*E2',1));
%__________________________________________________________________________
% Calculate total energy of the new state.
Xmu1=[X(i,1),Y(i,1),Z(i,1)];
Xmu2=[X(i,len),Y(i,len),Z(i,len)];
Xmu3=[X(i,Cen),Y(i,Cen),Z(i,Cen)];
PvMF1(i,1) = (kappa1 *(mu1*Xmu1'./sqrt(sum(Xmu1.^2))));
PvMF2(i,1) = (kappa2 *(mu2*Xmu2'./sqrt(sum(Xmu2.^2))));
PvMF3(i,1) = (kappa3 *(mu3*Xmu3'./sqrt(sum(Xmu3.^2))));

ll   = length(find((Xsub2.^2 + Ysub2.^2 + Zsub2.^2 <= (Radius+1)^2)));
E_c  = (ll)*0 + ((len) -ll)*U;
Energy_WLC(i,1) = Energy_WLC(i-1,1) + (k/a)*En2 - (k/a)*En1;
Energy(i,1) = Energy_WLC(i,1) + E_c + 3*sum(rx.^2 + ry.^2 +rz.^2)/(2*a^2)  - (PvMF1(i,1)+PvMF2(i,1)+PvMF3(i,1));

%Set conditions for acceptance.
condition1 = rand(1,1);
dE         = abs(Energy(i,1) - Energy(i-1,1));
if ((Energy(i,1)<Energy(i-1,1)) || (condition1<exp(-beta*dE)))
X(i,:)  = Xsub2;
Y(i,:)  = Ysub2;
Z(i,:)  = Zsub2;
%Energy(i,1)     = Energy(i-1,1);
%Energy_WLC(i,1) = Energy_WLC(i-1,1);
else
X(i,:)  = X(i-1,:);
Y(i,:)  = Y(i-1,:);
Z(i,:)  = Z(i-1,:);
Energy(i,1)     = Energy(i-1,1);
Energy_WLC(i,1) = Energy_WLC(i-1,1);
PvMF1(i,1) = PvMF1(i-1,1);
PvMF2(i,1) = PvMF2(i-1,1);
PvMF3(i,1) = PvMF3(i-1,1);
end

elseif stepchoice > (2/3)
Xsub = X(i-1,Cen:len);
Ysub = Y(i-1,Cen:len);
Zsub = Z(i-1,Cen:len);
%
[Xs,Ys,Zs] = ArmRotate(Xsub,Ysub,Zsub);
%
X(i,:) = X(i-1,:);
Y(i,:) = Y(i-1,:);
Z(i,:) = Z(i-1,:);
%
X(i,Cen:len) = Xs;
Y(i,Cen:len) = Ys;
Z(i,Cen:len) = Zs;
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
   E1 = [[rx./Delta]', ...
         [ry./Delta]', ...
         [rz./Delta]'];
%     
   E2 = [[rx1./Delta1]', ...
         [ry1./Delta1]', ...
         [rz1./Delta1]'];
%
En1(1,1) = sum(diag(E1*E1',1));
En2(1,1) = sum(diag(E2*E2',1));
%__________________________________________________________________________
% Calculate total energy of the new state.
Xmu1=[X(i,1),Y(i,1),Z(i,1)];
Xmu2=[X(i,len),Y(i,len),Z(i,len)];
Xmu3=[X(i,Cen),Y(i,Cen),Z(i,Cen)];
PvMF1(i,1) = (kappa1 *(mu1*Xmu1'./sqrt(sum(Xmu1.^2))));
PvMF2(i,1) = (kappa2 *(mu2*Xmu2'./sqrt(sum(Xmu2.^2))));
PvMF3(i,1) = (kappa3 *(mu3*Xmu3'./sqrt(sum(Xmu3.^2))));

ll   = length(find((X(i,:).^2 + Y(i,:).^2 + Z(i,:).^2 <= (Radius+1)^2)));
E_c  = (ll)*0 + ((len) -ll)*U;
Energy_WLC(i,1) = Energy_WLC(i-1,1) + (k/a)*En2 - (k/a)*En1;
Energy(i,1) = Energy_WLC(i,1) + E_c + 3*sum(rx.^2 + ry.^2 +rz.^2)/(2*a^2) - (PvMF1(i,1)+PvMF2(i,1)+PvMF3(i,1));

%Set conditions for acceptance.
condition1 = rand(1,1);
condition2 = rand(1,1);
condition3 = rand(1,1);
condition4 = rand(1,1);
dE         = abs(Energy(i,1) - Energy(i-1,1));
%__________________________________________________________________________
%Accept or reject new state.
if ((Energy(i,1)<Energy(i-1,1)) || (condition1<exp(-beta*dE)))
else
X(i,:) = X(i-1,:);
Y(i,:) = Y(i-1,:);
Z(i,:) = Z(i-1,:);
Energy(i,1)     = Energy(i-1,1);
Energy_WLC(i,1) = Energy_WLC(i-1,1);
PvMF1(i,1)       = PvMF1(i-1,1);
PvMF2(i,1)       = PvMF2(i-1,1);
PvMF3(i,1)       = PvMF3(i-1,1);
end

%X(i,:) = X(i-1,:);
%Y(i,:) = Y(i-1,:);
%Z(i,:) = Z(i-1,:);
%Energy(i,:) = Energy(i-1,:);
%Energy_WLC(i,:) = Energy_WLC(i-1,:);

Xsub = X(i,Cen:len);
Ysub = Y(i,Cen:len);
Zsub = Z(i,Cen:len);
[Xs,Ys,Zs] = CrankshaftRotate(Xsub,Ysub,Zsub);


Xsub2 = X(i,:);
Ysub2 = Y(i,:);
Zsub2 = Z(i,:);
Xsub2(Cen:len) = Xs;
Ysub2(Cen:len) = Ys;
Zsub2(Cen:len) = Zs;

rx = diff(Xsub2);
ry = diff(Ysub2);
rz = diff(Zsub2);
rx1 = diff(X(i,:));
ry1 = diff(Y(i,:));
rz1 = diff(Z(i,:));
Delta = sqrt(rx.^2 + ry.^2 + rz.^2);
Delta1 = sqrt(rx1.^2 + ry1.^2 + rz1.^2);
%__________________________________________________________________________
   E1 = [[rx./Delta]', ...
         [ry./Delta]', ...
         [rz./Delta]'];
     
   E2 = [[rx1./Delta1]', ...
         [ry1./Delta1]', ...
         [rz1./Delta1]'];

En1(1,1) = sum(diag(E1*E1',1));
En2(1,1) = sum(diag(E2*E2',1));
%__________________________________________________________________________
% Calculate total energy of the new state.
Xmu1=[X(i,1),Y(i,1),Z(i,1)];
Xmu2=[X(i,len),Y(i,len),Z(i,len)];
Xmu3=[X(i,Cen),Y(i,Cen),Z(i,Cen)];
PvMF1(i,1) = (kappa1 *(mu1*Xmu1'./sqrt(sum(Xmu1.^2))));
PvMF2(i,1) = (kappa2 *(mu2*Xmu2'./sqrt(sum(Xmu2.^2))));
PvMF3(i,1) = (kappa3 *(mu3*Xmu3'./sqrt(sum(Xmu3.^2))));
ll   = length(find((Xsub2.^2 + Ysub2.^2 + Zsub2.^2 <= (Radius+1)^2)));
E_c  = (ll)*0 + ((len) -ll)*U;
Energy_WLC(i,1) = Energy_WLC(i-1,1) + (k/a)*En2 - (k/a)*En1;
Energy(i,1) = Energy_WLC(i,1) + E_c + 3*sum(rx.^2 + ry.^2 +rz.^2)/(2*a^2) - (PvMF1(i,1)+PvMF2(i,1)+PvMF3(i,1)) ;

%Set conditions for acceptance.
condition1 = rand(1,1);
dE         = abs(Energy(i,1) - Energy(i-1,1));
if ((Energy(i,1)<Energy(i-1,1)) || (condition1<exp(-beta*dE)))
X(i,:)  = Xsub2;
Y(i,:)  = Ysub2;
Z(i,:)  = Zsub2;
%Energy(i,1)     = Energy1;
%Energy_WLC(i,1) = Energy1_WLC;
else
X(i,:)  = X(i-1,:);
Y(i,:)  = Y(i-1,:);
Z(i,:)  = Z(i-1,:);
Energy(i,1)     = Energy(i-1,1);
Energy_WLC(i,1) = Energy_WLC(i-1,1);
PvMF1(i,1) = PvMF1(i-1,1);
PvMF2(i,1) = PvMF2(i-1,1);
PvMF3(i,1) = PvMF3(i-1,1);
end

else
Xsub                        = [X(i-1,:)',Y(i-1,:)',Z(i-1,:)'];
[Pr] = Quaternion3(thet(i,1),uu(i,:),Xsub);
X(i,:) = Pr(:,1)';
Y(i,:) = Pr(:,2)';
Z(i,:) = Pr(:,3)';

Xmu1            = [X(i,1),Y(i,1),Z(i,1)];
Xmu2            = [X(i,len),Y(i,len),Z(i,len)];
Xmu3            = [X(i,Cen),Y(i,Cen),Z(i,Cen)];
PvMF1(i,1)       = (kappa1 *(mu1*Xmu1'./sqrt(sum(Xmu1.^2))));
PvMF2(i,len)       = (kappa2 *(mu2*Xmu2'./sqrt(sum(Xmu2.^2))));
PvMF3(i,Cen)       = (kappa3 *(mu3*Xmu3'./sqrt(sum(Xmu3.^2))));

Energy_WLC(i,1) = Energy_WLC(i-1,1);
Energy(i,1)     = Energy(i-1,1) + (PvMF1(i-1,1)+PvMF2(i-1,1)+PvMF3(i-1,1)) - (PvMF1(i,1)+PvMF2(i,1)+PvMF3(i,1)) ;

%Set conditions for acceptance.
condition1 = rand(1,1);
dE         = abs(Energy(i,1) - Energy(i-1,1));
if ((Energy(i,1)<Energy(i-1,1)) || (condition1<exp(-beta*dE)))
else
X(i,:) = X(i-1,:);
Y(i,:) = Y(i-1,:);
Z(i,:) = Z(i-1,:);
Energy(i,1)     = Energy(i-1,1);
Energy_WLC(i,1) = Energy_WLC(i-1,1);
PvMF1(i,1)       = PvMF1(i-1,1);
PvMF2(i,1)       = PvMF2(i-1,1);
PvMF3(i,1)       = PvMF3(i-1,1);
end
end
if ~mod(i, floor(No/20))
        %update ~20 times in total
        if UpdateProgress(['Running Analysis    ' num2str(i) '/' num2str(No)])
            return;
        end

    end

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
