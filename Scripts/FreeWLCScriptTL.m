%MCMC method for generating samples from the joint distirbution of a WLC
%confined to a spherical volume.
addpath('./Functions','./netlab')

x = linspace(0,2*pi,98);
y = sin(x) + 1.01;
y(find(y<1.2)) = 0.05;
No      = 200000;                 %Number of steps in the Markov Chain.
len     = 100;                    %Number of beads in the polymer chain.
Radius  = 1600;                   %Radius of the Nuclear Periphery (um).
U       = 0;                  %Potential for beads outside of the NP.
kappa   = 1e-10;                  %Telomere clustering parameter.
k       = y*230000*0.34/10;     %Persistence length (um).
a       = 230000*0.34/(10*100);   %Bond vector length.
VMF     = [kappa,0,0,1];          %von Mises-Fisher parameters.

Parameters             = [No,len,Radius,U,a,k];  %Group parameters.
[x0,y0,z0,theta,phi] = tMultivMF(0,0,1,kappa,Radius,1);
[x1,y1,z1] = RandFlightCon([x0,y0,z0],[x0,y0,z0],a,len);


[X,Y,Z,Energy] = FreeWLCTV(Parameters,x1,y1,z1);

%Plot example configuration after N steps.
PolymerPic(X,Y,Z,Radius,No);