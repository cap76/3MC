%Example script for generating random samples from a Rabl-like
%conformation. Uses the `Rabl.m' function.
if ~isdeployed
    addpath(genapth('../'));
end

No      = 200000;                %Number of steps in the Markov Chain.
len     = 100;                  %Number of beads in the polymer chain.
Radius  = 1600;                 %Radius of the Nuclear Periphery (NP).
U       = 10000;                %Potential for beads outside of the NP.
Phi     = 10;
kappa   = 1e-10;                %Centromere clustering parameter.
k       = .5*230000*0.34/Phi;    %Persistence length.
a       = 230000*0.34/(Phi*(len-1)); %Bond vector length.
VMF     = [kappa,0,0,1];        %von Mises-Fisher parameters.
Cen     = 66;                   %Position of centromere (bead number).

Parameters1     = [No,len,k,Radius,U,a,Cen]; %Group parameters.

%Generate random walk intitial conformation.
[x0,y0,z0,theta,phi] = tMultivMF(0,0,1,kappa,Radius,1);
[x1,y1,z1] = RandFlightCon([x0,y0,z0],[x0,y0,z0],a,len);
%Ensure that bead Cen is found on the surface of nucleus.
dx = x1(1,1);
dy = y1(1,1);
dz = z1(1,1);
x1 = x1 - dx;
y1 = y1 - dy;
z1 = z1 - dz + Radius;

%Pass variables over to the ``Rabl'' function.
[net1] = Celegans(Parameters1,x1,y1,z1,VMF);
%[net2] = Celegans(Parameters1,x1,y1,z1,VMF);
return
Width   = -0.0056*Phi^2 + 1.6*Phi -3.4; 
Width   = 100;
%Plot example configuration after N steps.
PolymerPic(net1.Parameters,net1.X,net1.Y,net1.Z,Width);


[D,V] = HDM(net1,net2);
figure
imagesc(D)
xlabel('Position along Chromosome I')
ylabel('Position along Chromosome Is Homologue')
title('CDM')

[D2,V2] = LDM(net1);
figure
imagesc(D2)
xlabel('Position along Chromosome I')
ylabel('Position along Chromosome I')
title('LDM')

%Parameters for pairing model
RN   = 50;    %Size of recombination nodules
p    = 0.9;   %Binomial parameter (controls number of initialising loci)
Type = 'Centromere'; %Type of pairing: Telomere-proximal = 'Telomere' 
                   %                 Centromere-proximal = 'Centromere'
                   %                 Random = anything other string

[probpair] = Pairing(net1,net2,[RN,p]',Type,Parameters1);
disp(['Probability of pairing = ' num2str(probpair*100) '%'])