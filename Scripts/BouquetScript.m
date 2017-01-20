%Example script for generating random samples from a Bouquet-like
%formation. Uses the `Bouquet' function.
addpath('./Functions','./netlab')
mkdir('./Bouquet')
    
No      = 2000000;                 %Number of steps in the Markov Chain.
len     = 100;                    %Number of beads in the polymer chain.
Radius  = 1600;                   %Radius of the Nuclear Periphery (um).
U       = 10000;                  %Potential for beads outside of the NP.
kappa   = 5;%e-10;                  %Telomere clustering parameter.
k       = 1*230000*0.34/30;     %Persistence length (um).
a       = 230000*0.34/(15*100);   %Bond vector length.
VMF     = [kappa,0,0,1];          %von Mises-Fisher parameters.

Parameters1             = [No,len,k,Radius,U,a];  %Group parameters.
[x0,y0,z0,theta,phi] = tMultivMF(0,0,1,kappa,Radius,1);
[x1,y1,z1] = RandFlightCon([x0,y0,z0],[x0,y0,z0],a,len);

%Pass variables over to the ``Bouquet'' function.
%[X,Y,Z,En] = BouquetPrior(Parameters1,x1,y1,z1,VMF);
 tic
 [X,Y,Z,En] = Bouquet(Parameters1,x1,y1,z1,VMF);
 toc
   %Save trajectories
   %save(['./Bouquet/X.mat'],'X')
   %save(['./Bouquet/Y.mat'],'Y')
   %save(['./Bouquet/Z.mat'],'Z')
   %save(['./Bouquet/En.mat'],'En')

%Plot example configuration after N steps.
PolymerPic(X,Y,Z,Radius,No);