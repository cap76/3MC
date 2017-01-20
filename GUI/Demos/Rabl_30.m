function [] = Rabl_30(Chrom,Comp);
%Example script for generating random samples from a Bouquet-like
%formation. Uses the `Bouquet' function.

rand('state',10)
randn('state',10)

switch Chrom
    case 1
        Dis = 230000;    
        Cen = 151000;
    case 2
        Dis = 813000;
        Cen = 238000;        
    case 3
        Dis = 317000;
        Cen = 114000;        
    case 4
        Dis = 1532000;
        Cen = 450000;        
    case 5
        Dis = 577000;
        Cen = 152000;        
    case 6
        Dis = 270000;
        Cen = 149000;        
    case 7
        Dis = 1091000;
        Cen = 497000;        
    case 8
        Dis = 563000;
        Cen = 106000;
    case 9
        Dis = 440000;
        Cen = 356000;
    case 10
        Dis = 746000;
        Cen = 436000;
    case 11
        Dis = 666000;
        Cen = 440000;        
    case 12
        Dis = 1079000;
        Cen = 150000;        
    case 13
        Dis = 924000;
        Cen = 268000;        
    case 14
        Dis = 784000;
        Cen = 629000;        
    case 15
        Dis = 109000;
        Cen = 327000;        
    case 16
        Dis = 948000;
        Cen = 556000;        
end    

addpath(genpath('../'))


%addpath('../Functions','../netlab')
mkdir('../Rabl')
    
PhiStar = 40;
Persistence = [0.1000,0.2556, 1]; %linspace(0.1,1.5,10); 
Spread      = [1e-10,5,20];

for kk = 1:3
for j = 1:3

No      = 100000;                 %Number of steps in the Markov Chain.
len     = 300;                    %Number of beads in the polymer chain.
Radius  = 1600;                   %Radius of the Nuclear Periphery (nm).
U       = 10000;                  %Potential for beads outside of the NP.
Phi     = Comp;
kappa   = Spread(1,kk);                 %Telomere clustering parameter.
k       = Persistence(1,j)*230000*0.34/PhiStar;  %Persistence length (nm).
a       = Dis*0.34/(Phi*(len-1));         %Bond vector length.
VMF     = [kappa,0,0,1;kappa,0,0,1];         %von Mises-Fisher parameters.
CEN     = double(int64((Cen./Dis * (len-1)) + 1)); %Centromere position

Parameters1     = [No,len,k,Radius,U,a,CEN]; %Group parameters.


%Generate random walk intitial conformation.
[x0,y0,z0,theta,phi] = tMultivMF(0,0,1,kappa,Radius,1);
[x1,y1,z1] = RandFlightCon([x0,y0,z0],[x0,y0,z0],a,len);
%Ensure that bead Cen is found on the surface of nucleus.
dx = x1(1,CEN);
dy = y1(1,CEN);
dz = z1(1,CEN);
x1 = x1 - dx;
y1 = y1 - dy;
z1 = z1 - dz + Radius;

%Pass variables over to the ``Rabl'' function.
[Data ] = Rabl(Parameters1,x1,y1,z1,VMF);
    
    
for i = 1:25

    [Data   ] = Rabl(Parameters1,Data.X(end,:),Data.Y(end,:),Data.Z(end,:),VMF);

    if i>15
        save(['../Rabl/Rabl_' num2str(Comp) 'nm_Chrom_' num2str(Chrom) '_' num2str(kk) '_' num2str(j) '_' num2str(i) '.mat'],'Data')
    else    
    end


end

end

end

%exit