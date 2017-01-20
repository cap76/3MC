%Example script for generating random samples from a Bouquet-like
%formation. Uses the `Bouquet' function.
addpath('./Functions','./netlab')
mkdir('./Bouquet')
    
No      = 20000;                 %Number of steps in the Markov Chain.
len     = 50;                    %Number of beads in the polymer chain.
Radius  = 1600;                   %Radius of the Nuclear Periphery (um).
U       = 10000;                  %Potential for beads outside of the NP.
kappa   = 5;%e-10;                  %Telomere clustering parameter.
k       = 1*230000*0.34/30;     %Persistence length (um).
a       = 230000*0.34/(15*100);   %Bond vector length.
VMF     = [kappa,0,0,1];          %von Mises-Fisher parameters.

Parameters1{1}             = [No,len,0.01*k,Radius,U,a];  %Group parameters.
Parameters1{2}             = [No,len,0.1*k,Radius,U,a];
Parameters1{3}             = [No,len,0.5*k,Radius,U,a];
Parameters1{4}             = [No,len,k,Radius,U,a];


parfor m = 1:4

[x0,y0,z0,theta,phi] = tMultivMF(0,0,1,kappa,Radius,1);
[x1,y1,z1] = RandFlightCon([x0,y0,z0],[x0,y0,z0],a,len);    
    
%Pass variables over to the ``Bouquet'' function.
[Data{m}] = Bouquet(Parameters1{m},x1,y1,z1,VMF);

end

   %Save trajectories
   %save(['./Bouquet/X.mat'],'X')
   %save(['./Bouquet/Y.mat'],'Y')
   %save(['./Bouquet/Z.mat'],'Z')
   %save(['./Bouquet/En.mat'],'En')

%Plot example configuration after N steps.
%PolymerPic(X,Y,Z,Radius,No);