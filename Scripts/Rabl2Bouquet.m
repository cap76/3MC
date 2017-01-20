function [X,Y,Z,Energy_WLC] = Bouquet2RablPrior(Parameters,X0x,Y0x,Z0x,VMF);
%NOTE: NEEDS TO ALLOW ROTATION OF ARM.
%Description (version 3)
%
%[X,Y,Z,Energy] = BOUQUET3(Parameters,x,y,z,dx,dy,dz) generates an arbtrary 
%length Markov Chain for chromosomes in a Bouquet formation with boundary 
%conditions drawn from a von Mises-Fisher distribution. Three types of 
%update are defined: (i) a crankshaft rotation about an axis connecting 
%two randomly determined bonds; (ii) a crankshaft rotation of one arm; 
%(iii) rotation of chromosomes via random axis. Example of use: BouquetScript3.

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
No          = Parameters(1,1);    % No. steps in Markov chain.
mu1          = VMF(1,2:4);
mu2          = VMF(2,2:4);
mu3          = VMF(3,2:4);
Xmu1        = [X0x(1,1),Y0x(1,1),Z0x(1,1)];
Xmu2        = [X0x(1,len),Y0x(1,len),Z0x(1,len)];
Xmu3        = [X0x(1,Cen),Y0x(1,Cen),Z0x(1,Cen)];
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
Energy(i,1)= Energy_WLC(i,1) + E_c + 3*sum(rx.^2 + ry.^2 +rz.^2)/(2*a^2) + (PvMF1(i-1,1)+PvMF2(i-1,1)+PvMF3(i-1,1)) - (PvMF1(i,1)+PvMF2(i,1)+PvMF3(i,1));

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
Energy1_WLC = Energy_WLC(i,1) + (k/a)*En2 - (k/a)*En1;
Energy1 = Energy1_WLC + E_c + 3*sum(rx.^2 + ry.^2 +rz.^2)/(2*a^2)  + (PvMF1(i-1,1)+PvMF2(i-1,1)+PvMF3(i-1,1)) - (PvMF1(i,1)+PvMF2(i,1)+PvMF3(i,1));

%Set conditions for acceptance.
condition1 = rand(1,1);
dE         = abs(Energy1 - Energy(i,1));
if ((Energy1<Energy(i,1)) || (condition1<exp(-beta*dE)))
X(i,:)  = Xsub2;
Y(i,:)  = Ysub2;
Z(i,:)  = Zsub2;
Energy(i,1)     = Energy1;
Energy_WLC(i,1) = Energy1_WLC;
else
X(i,:)  = X(i,:);
Y(i,:)  = Y(i,:);
Z(i,:)  = Z(i,:);
Energy(i,1)     = Energy(i,1);
Energy_WLC(i,1) = Energy_WLC(i,1);
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
Energy(i,1) = Energy_WLC(i,1) + E_c + 3*sum(rx.^2 + ry.^2 +rz.^2)/(2*a^2) + (PvMF1(i-1,1)+PvMF2(i-1,1)+PvMF3(i-1,1)) - (PvMF1(i,1)+PvMF2(i,1)+PvMF3(i,1));

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
Energy1_WLC = Energy_WLC(i,1) + (k/a)*En2 - (k/a)*En1;
Energy1 = Energy1_WLC + E_c + 3*sum(rx.^2 + ry.^2 +rz.^2)/(2*a^2) + (PvMF1(i-1,1)+PvMF2(i-1,1)+PvMF3(i-1,1)) - (PvMF1(i,1)+PvMF2(i,1)+PvMF3(i,1)) ;

%Set conditions for acceptance.
condition1 = rand(1,1);
dE         = abs(Energy1 - Energy(i,1));
if ((Energy1<Energy(i,1)) || (condition1<exp(-beta*dE)))
X(i,:)  = Xsub2;
Y(i,:)  = Ysub2;
Z(i,:)  = Zsub2;
Energy(i,1)     = Energy1;
Energy_WLC(i,1) = Energy1_WLC;
else
X(i,:)  = X(i,:);
Y(i,:)  = Y(i,:);
Z(i,:)  = Z(i,:);
Energy(i,1)     = Energy(i,1);
Energy_WLC(i,1) = Energy_WLC(i,1);
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
dE         = abs(Energy(i,1) - Energy(i,1));
if ((Energy(i,1)<Energy(i,1)) || (condition1<exp(-beta*dE)))
else
X(i,:) = X(i-1,:);
Y(i,:) = Y(i-1,:);
Z(i,:) = Z(i-1,:);
PvMF1(i,1)       = PvMF1(i-1,1);
PvMF2(i,1)       = PvMF2(i-1,1);
end
end
end
%subsample       = 70000:40:109999;%linspace(50000,100000,5001);
%X = X(subsample,:); Y = Y(subsample,:); Z = Z(subsample,:); Energy_WLC = Energy_WLC(subsample,1);

%Addressx1 = (['C:\WLCData\X' num2str(Chromval) '.mat']);
%Addressy1 = (['C:\WLCData\Y' num2str(Chromval) '.mat']);
%Addressz1 = (['C:\WLCData\Z' num2str(Chromval) '.mat']);
%AddressE1 = (['C:\WLCData\E' num2str(Chromval) '.mat']);
%save(Addressx1,'X')
%save(Addressy1,'Y')
%save(Addressz1,'Z')