%Example script for generating random samples from a Rabl/Bouquet-like
%formation. Uses the `Rabl2Bouquet' function.

if ~isdeployed
    addpath(genapth('../'));
end

No      = 200000;               %Number of steps in the Markov Chain.
len     = 100;                  %Number of beads in the polymer chain.
Radius  = 1600;                 %Radius of the Nuclear Periphery (NP).
U       = 10000;                %Potential for beads outside of the NP.
Phi     = 10;
kappa   = 1e-10;                %Telomere clustering parameter.
k       = .5*230000*0.34/Phi;    %Persistence length.
a       = 230000*0.34/(Phi*(len-1)); %Bond vector length.
VMF     = [kappa,0,0,1;kappa,0,0,1;kappa,0,0,1];        %von Mises-Fisher parameters.
Cen     = 66;

Parameters1     = [No,len,k,Radius,U,a,Cen];  %All parameters.

%Generate random walk intitial conformation.
[x0,y0,z0,theta,phi]    = tMultivMF(0,0,1,100,Radius,1);
[x01,y01,z01,theta,phi] = tMultivMF(0,0,1,100,Radius,1);
[x02,y02,z02,theta,phi] = tMultivMF(0,0,1,100,Radius,1);
[x1,y1,z1] = RandFlightCon([x0,y0,z0],[x01,y01,z01],a,100-Cen);
[x2,y2,z2] = RandFlightCon([x01,y01,z01],[x02,y02,z02],a,Cen+1);

x = [x1,x2(1,2:length(x2))];
y = [y1,y2(1,2:length(x2))];
z = [z1,z2(1,2:length(x2))];

%Pass variables over to the ``Rabl2BouquetPrior'' function.
[net1] = Rabl2Bouquet(Parameters1,x,y,z,VMF);
[net2] = Rabl2Bouquet(Parameters1,x,y,z,VMF);

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