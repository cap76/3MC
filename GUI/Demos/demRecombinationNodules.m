if ~isdeployed
    addpath(genpath('../'));
end


%Randomly position TF factories
X   = rand(1000,1)*3200 - 1600;
Y   = rand(1000,1)*3200 - 1600;
Z   = rand(1000,1)*3200 - 1600;
Dis = X.^2+Y.^2+Z.^2;
Ind = find(Dis<(1600-50).^2);

%Pick only X,Y and Z coordinates that lie in nucleus
X = X(Ind(1:21));
Y = Y(Ind(1:21));
Z = Z(Ind(1:21));


FactoryRad = 50; %Radius of TF factory
No      = 500000;                  %Number of steps in the Markov Chain.
len     = 100;                    %Number of beads in the polymer chain.
Radius  = 1600;                   %Radius of the Nuclear Periphery (um).
U       = 10000;                  %Potential for beads outside of the NP.
kappa   = 5;                      %Telomere clustering parameter.
Phi     = 30;                     %Compaction factor

k       = .2*230000*0.34/Phi;           %Persistence length (um).
a       = 230000*0.34/(Phi*(len-1));   %Bond vector length.
VMF     = [kappa,0,0,1;kappa,0,0,1];   %von Mises-Fisher parameters.
Cen     = 66;                          %Centromere position
Parameters1             = [No,len,k,Radius,U,a,Cen];  %Group parameters.
[x0,y0,z0,theta,phi] = tMultivMF(0,0,1,kappa,Radius,1);
[x1,y1,z1] = RandFlightCon([x0,y0,z0],[x0,y0,z0],a,len);

%Generate some chromosome trajectories
NoConfigs = 5;
for i = 1:NoConfigs
    net{i} = Bouquet(Parameters1,x1,y1,z1,VMF);
end

%Calculate how often each locus falls within 50nm of a TF factory
N = cell(NoConfigs,1);
for j = 1:NoConfigs
    for factory = 1:21
        for i = 1:size(net{j}.X,2)
   
        x0   = repmat(X(factory,1),size(net{j}.X,1),1);
        y0   = repmat(Y(factory,1),size(net{j}.X,1),1);
        z0   = repmat(Z(factory,1),size(net{j}.X,1),1);         
        Dis1 = sqrt((net{j}.X(:,i)-x0).^2 + (net{j}.Z(:,i)-z0).^2 + (net{j}.Y(:,i)-y0).^2);

        N{j}(factory,i) = length(find(Dis1<(FactoryRad+30)));
        
        
        end
    end
end