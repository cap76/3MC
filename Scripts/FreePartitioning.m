PHI     = [11,30,100];
No      = 1000000;                %Number of steps in the Markov Chain.
len     = 100;                    %Number of beads in the polymer chain.
Radius  = 1600;                   %Radius of the Nuclear Periphery (um).
U       = 100000;                 %Potential for beads outside of the NP.
kappa   = [1e-10,1,5,10];%e-10;   %Telomere clustering parameter.
Cen     = 65;
kin     = linspace(0.1,1.5,5);     %Persistence length (um).
%
%theta = linspace(0,pi,5);
%mu     = [zeros(5,1),sin(theta)',cos(theta)']; %Telomere directional bias vector

    mkdir('./FreeN')


                           
                        MU1 = cell(5,4);
                        MU2 = cell(5,4);
                        MU3 = cell(5,4);  
                        
                        DIS1 = zeros(100,100);
                        DIS2 = zeros(100,100);
                        DIS3 = zeros(100,100);
                        
for ii = 1:size(PHI,2) %Permute through directions    
        disp(['Step ' num2str(ii) ' of ' num2str(size(mu,1))])        
      for jj = 1:size(kin,2) %Permute through directional variances                                        
                disp(['    Step ' num2str(jj) ' of ' num2str(size(kappa,2))])

                

             
                        aa = 230000*0.34/(len*PHI(,ii));
                        kk = aa*len*kin(1,jj);                


         Parameters     = [No,len,kk,Radius,U,aa];

%Generate intial configuration.
[x1,y1,z1] = RandFlightCon([0,0,0],[0,0,0],aa,len);

[X,Y,Z,Energy] = FreeWLC(Parameters,x1,y1,z1);       

load(['./FreeN/X1_l=' num2str(ii) '_j=' num2str(jj) '.mat'],'X')
load(['./FreeN/Y1_l=' num2str(ii) '_j=' num2str(jj) '.mat'],'Y')
load(['./FreeN/Z1_l=' num2str(ii) '_j=' num2str(jj) '.mat'],'Z')
load(['./FreeN/E1_l=' num2str(ii) '_j=' num2str(jj) '.mat'],'En')

                                
            end
end
        
save('./RablEctopicMap1.mat','MU1')
save('./RablAllelicMap1.mat','MU2')
save('./RablAllelicMap2.mat','MU3')
