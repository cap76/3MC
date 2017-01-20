%Example script for generating random samples from a Bouquet-like
%formation. Uses the `Bouquet' function.
addpath('./Functions','./netlab')
mkdir('./Bouquet')

x = linspace(0,2*pi,98);
y = sin(x) + 1.1;
y(find(y<1)) = 0.2;
No      = 200000;                 %Number of steps in the Markov Chain.
len     = 100;                    %Number of beads in the polymer chain.
Radius  = 1600;                   %Radius of the Nuclear Periphery (um).
U       = 1e20;                  %Potential for beads outside of the NP.
kappa   = 1e-10;                  %Telomere clustering parameter.
k       = y*230000*0.34/10;     %Persistence length (um).
a       = 230000*0.34/(10*100);   %Bond vector length.
VMF     = [kappa,0,0,1];          %von Mises-Fisher parameters.

Parameters1             = [No,len,Radius,U,a,k];  %Group parameters.
[x0,y0,z0,theta,phi] = tMultivMF(0,0,1,kappa,Radius,1);
[x1,y1,z1] = RandFlightCon([x0,y0,z0],[x0,y0,z0],a,len);

%Pass variables over to the ``Bouquet'' function.
%[X,Y,Z,En] = BouquetPrior(Parameters1,x1,y1,z1,VMF);
 tic
 [X,Y,Z,En] = BouquetTV(Parameters1,x1,y1,z1,VMF);
 toc
   %Save trajectories
   %save(['./Bouquet/X.mat'],'X')
   %save(['./Bouquet/Y.mat'],'Y')
   %save(['./Bouquet/Z.mat'],'Z')
   %save(['./Bouquet/En.mat'],'En')

%Plot example configuration after N steps.
PolymerPic(X,Y,Z,Radius,No);