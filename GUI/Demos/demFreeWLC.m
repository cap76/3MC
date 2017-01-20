%MCMC method for generating samples from the joint distirbution of a WLC
%confined to a spherical volume.
if ~isdeployed
    addpath(genapth('../'));
end

%Initial parameters.
No      = 200000;                 %Number of steps in the Markov Chain.
len     = 50;                     %Number of beads in the polymer chain.
Radius  = 1600;                   %Radius of the Nuclear Periphery (NP).
U       = 10000;                  %Potential for beads outside of the NP.
Phi     = 10;                     %Compaction factor
k       = 0.5*230000*0.34/Phi;       %Persistence length.
a       = 230000*0.34/(Phi*(len-1)); %Bond vector length.
kappa   = 1e-10;                     %Cluster parameter
Cen     = 66;                        %Position of centromere

VMF     = [NaN,NaN,NaN,NaN];      %von Mises-Fisher parameters. Not used for unconfined WLCs (dummy variables)

Parameters     = [No,len,k,Radius,U,a,Cen];

%Generate intial configuration.
[x0,y0,z0,theta,phi] = tMultivMF(0,0,1,kappa,Radius,1);
[x1,y1,z1] = RandFlightCon([x0,y0,z0],[x0,y0,z0],a,len);

[net1] = FreeWLC(Parameters,x1,y1,z1,VMF);
[net2] = FreeWLC(Parameters,x1,y1,z1,VMF);

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
Type = 'Random'; %Type of pairing: Telomere-proximal = 'Telomere' 
                   %                 Centromere-proximal = 'Centromere'
                   %                 Random = anything other string

[probpair] = Pairing(net1,net2,[RN,p]',Type,Parameters);
disp(['Probability of pairing = ' num2str(probpair*100) '%'])