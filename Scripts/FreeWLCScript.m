%MCMC method for generating samples from the joint distirbution of a WLC
%confined to a spherical volume.
addpath('./Functions','./netlab')

%Initial parameters.
No      = 20000;
len     = 100;
Radius  = 1600;
U       = 10000;
kappa   = 100;
k       = 1.5*230000*0.34/10;
a       = 230000*0.34/(10*100);
VMF     = [kappa,0,0,1];

Parameters     = [No,len,k,Radius,U,a];

%Generate intial configuration.
[x0,y0,z0,theta,phi] = tMultivMF(0,0,1,kappa,Radius,1);
[x1,y1,z1] = RandFlightCon([x0,y0,z0],[x0,y0,z0],a,len);

[X,Y,Z,Energy] = FreeWLC(Parameters,x1,y1,z1);

%Plot example configuration after N steps.
PolymerPic(X,Y,Z,Radius,No);