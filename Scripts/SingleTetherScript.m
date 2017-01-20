addpath('./Functions','./netlab')

No      = 20000;                 %Number of steps in the Markov Chain.
len     = 50;                    %Number of beads in the polymer chain.
Radius  = 1600;                   %Radius of the Nuclear Periphery (NP).
U       = 10000;                  %Potential for beads outside of the NP.
kappa   = 1e-10;                  %Telomere clustering parameter.
k       = 0.5*230000*0.34/30; %Persistence length.
a       = 230000*0.34/(30*100);   %Bond vector length.
VMF     = [kappa,0,0,1];          %von Mises-Fisher parameters.
Cen     = 50;

Parameters1     = [No,len,k,Radius,U,a,Cen];

%Initiate random flight initial conformation fixed at the origin.
%[x0,y0,z0,theta,phi] = tMultivMF(0,0,1,kappa,Radius,1);
x0 = 0; y0 = 0; z0 = 0;
[x1,y1,z1] = RandFlightCon([x0,y0,z0],[x0,y0,z0],a,len);
%Pass variables over to the ``Bouquet'' function.

[X,Y,Z,En] = SingleTether(Parameters1,x1,y1,z1);

%Plot example configuration after N steps.
PolymerPic(X,Y,Z,Radius,No);
