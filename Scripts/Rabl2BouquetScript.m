%Example script for generating random samples from a Bouquet-like
%formation. Uses the `Bouquet' function.

No      = 200000;               %Number of steps in the Markov Chain.
len     = 100;                  %Number of beads in the polymer chain.
Radius  = 1600;                 %Radius of the Nuclear Periphery (NP).
U       = 10000;                %Potential for beads outside of the NP.
kappa   = 30;%1e-9;                  %Telomere clustering parameter.
k       = .5*230000*0.34/10;   %Persistence length.
a      = 230000*0.34/(10*100);  %Bond vector length.
VMF     = [kappa,0,0,1;kappa,0,0,1;kappa,0,0,1];        %von Mises-Fisher parameters.
Cen = 20;
Parameters1     = [No,len,k,Radius,U,a,Cen];  %All parameters.

%Generate random walk intitial conformation.
[x0,y0,z0,theta,phi] = tMultivMF(0,0,1,100,Radius,1);
[x01,y01,z01,theta,phi] = tMultivMF(0,0,1,100,Radius,1);
[x02,y02,z02,theta,phi] = tMultivMF(0,0,1,100,Radius,1);
[x1,y1,z1] = RandFlightCon([x0,y0,z0],[x01,y01,z01],a,100-Cen);
[x2,y2,z2] = RandFlightCon([x01,y01,z01],[x02,y02,z02],a,Cen+1);
x = [x1,x2(1,2:length(x2))];
y = [y1,y2(1,2:length(x2))];
z = [z1,z2(1,2:length(x2))];

%Pass variables over to the ``Bouquet'' function.
%profile on
%tic
[X1,Y1,Z1,Energy1] = Rabl2Bouquet(Parameters1,x,y,z,VMF);
%toc
%profile off
%profile viewer


%Note: Centromere is NOT exploring the NP. Telomere's are.