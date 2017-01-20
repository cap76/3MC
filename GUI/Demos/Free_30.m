function [] = Free_30(Chrom,Comp);
%Example script for generating random samples from a Bouquet-like
%formation. Uses the `Bouquet' function.


switch Chrom
    case 1
        Dis = 230000;    
    case 2
        Dis = 813000;
    case 3
        Dis = 317000;
    case 4
        Dis = 1532000;
    case 5
        Dis = 577000;
    case 6
        Dis = 270000;
    case 7
        Dis = 1091000;
    case 8
        Dis = 563000;
    case 9
        Dis = 440000;
    case 10
        Dis = 746000;
    case 11
        Dis = 666000;
    case 12
        Dis = 1079000;
    case 13
        Dis = 924000;
    case 14
        Dis = 784000;
    case 15
        Dis = 109000;
    case 16
        Dis = 948000;
end
    
addpath(genpath('../'))

%addpath('../Functions','../netlab')
mkdir('../Free')
    
PhiStar = 40;
Persistence = [0.1000,0.2556, 1]; %linspace(0.1,1.5,10); 
Spread      = [1e-10,5,20];

for kk = 2%1:3
for j = [2]

No      = 100000;                 %Number of steps in the Markov Chain.
len     = 300;                    %Number of beads in the polymer chain.
Radius  = 1600;                   %Radius of the Nuclear Periphery (nm).
U       = 10000;                  %Potential for beads outside of the NP.
Phi     = Comp;
kappa   = Spread(1,kk);                 %Telomere clustering parameter.
k       = Persistence(1,j)*230000*0.34/PhiStar;  %Persistence length (nm).
a       = Dis*0.34/(Phi*(len-1));         %Bond vector length.
VMF     = [kappa,0,0,1;kappa,0,0,1];         %von Mises-Fisher parameters.

Parameters           = [No,len,k,Radius,U,a,50];  %Group parameters.
[x0,y0,z0,theta,phi] = tMultivMF(0,0,1,kappa,Radius,1);
[x1,y1,z1] = RandFlightCon([x0,y0,z0],[x0,y0,z0],a,len);
[Data   ] = FreeWLC(Parameters,x1,y1,z1,VMF);
    
    
for i = 1:25

    [Data   ] = FreeWLC(Parameters,Data.X(end,:),Data.Y(end,:),Data.Z(end,:),VMF);

    if i>15
        save(['../Free/Free_' num2str(Comp) 'nm_Chrom_' num2str(Chrom) '_' num2str(kk) '_' num2str(j) '_' num2str(i) '.mat'],'Data')
    else    
    end


end

end

end

%exit