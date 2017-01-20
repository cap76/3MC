%Example script for generating random samples from a Bouquet-like
%formation. Uses the `Bouquet' function.

if ~isdeployed
    addpath(genapth('../'));
end
    
No      = 200000;                 %Number of steps in the Markov Chain.
len     = 100;                    %Number of beads in the polymer chain (segments that make the chr).
Radius  = 1600;                   %Radius of the Nuclear Periphery (nm).
U       = 10000;                  %Potential for beads outside of the NP.
Phi     = 10;                     %Compaction factor
kappa   = 10;                     %Telomere clustering parameter.
Cen     = 66;                     %Position of centromere
%len > 10, No can be 1, U arbitrarily high, ought to be inf
%phi < 500, kappa, 0 < kappa <= 20
%phi - chagne to compaction factor

%These 2 don't need to be entered. maybe 230000, 1.5 and 0.34 do
%1.5 is persistence length( > 0), 230000 - chr length(bp)
k       = 1.5*230000*0.34/Phi;     %Persistence length (nm) of chr. k is rigidity const
a       = 230000*0.34/(Phi*(len-1));   %Bond vector length.l ength of each chr segment

%0 0 1 could be entered? direction of clustering (x, y, z). 3 elements, and sum of
%squares of them must be one

VMF     = [kappa,0,0,1;kappa,0,0,1];          %von Mises-Fisher parameters.

Parameters1             = [No,len,k,Radius,U,a];  %Group parameters.
%theta and phi neve rused
[x0,y0,z0,theta,phi] = tMultivMF(0,0,1,kappa,Radius,1);% last para means ONE LOOP INTERATION
[x1,y1,z1] = RandFlightCon([x0,y0,z0],[x0,y0,z0],a,len);

[net1 ] = conditionalBouquet(Parameters1,x1,y1,z1,VMF);
[net2 ] = conditionalBouquet(Parameters1,x1,y1,z1,VMF);

Width   = -0.0056*Phi^2 + 1.6*Phi -3.4; 
%
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
Type = 'Telomere'; %Type of pairing: Telomere-proximal = 'Telomere' 
                   %                 Centromere-proximal = 'Centromere'
                   %                 Random = anything other string

[probpair] = Pairing(net1,net2,[RN,p]',Type,Parameters1);
disp(['Probability of pairing = ' num2str(probpair*100) '%'])