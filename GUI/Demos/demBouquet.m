%Example script for generating random samples from a Bouquet-like
%formation. Uses the `Bouquet' function.

if ~isdeployed
    addpath(genpath('../'));
end
    
No      = 500000;                  %Number of steps in the Markov Chain.
len     = 100;                    %Number of beads in the polymer chain.
Radius  = 1600;                   %Radius of the Nuclear Periphery (um).
U       = 10000;                  %Potential for beads outside of the NP.
kappa   = 5;                      %Telomere clustering parameter.
Phi     = 10;                     %Compaction factor
k       = .2*230000*0.34/Phi;           %Persistence length (um).
a       = 230000*0.34/(Phi*(len-1));   %Bond vector length.
VMF     = [kappa,0,0,1;kappa,0,0,1];   %von Mises-Fisher parameters.
Cen     = 66;                          %Centromere position
Parameters1             = [No,len,k,Radius,U,a,Cen];  %Group parameters.
[x0,y0,z0,theta,phi] = tMultivMF(0,0,1,kappa,Radius,1);
[x1,y1,z1] = RandFlightCon([x0,y0,z0],[x0,y0,z0],a,len);


tic,[net1] = Bouquet(Parameters1,x1,y1,z1,VMF);toc
tic,[net2] = BouquetOpt(Parameters1,x1,y1,z1,VMF);toc
net3{1} = net1;
net4{1} = net2;

 
Width   = -0.0056*Phi^2 + 1.6*Phi -3.4; 
%
%Plot example configuration after N steps.
PolymerPic(net1.Parameters,net1.X,net1.Y,net1.Z,Width);


[D] = CDM(net3,net4);
figure
imagesc(D)
xlabel('Position along Chromosome I')
ylabel('Position along Chromosome Is Homologue')
title('CDM')

[D2] = LDM(net1);
figure
imagesc(D2)
xlabel('Position along Chromosome I')
ylabel('Position along Chromosome I')
title('LDM')

%Parameters for pairing model
RN   = 50;    %Size of recombination nodules
p    = 0.9;   %Binomial parameter (controls number of initialising loci)
Type = 'Telomere'; %Type of pairing: Telomere-proximal = 'Telomere' 
                   %                 Centromere-proximal = 'Centromere'
                   %                 Random = anything other string

[probpair] = Pairing(net1,net2,[RN,p]',Type,Parameters1);
disp(['Probability of pairing = ' num2str(probpair*100) '%'])