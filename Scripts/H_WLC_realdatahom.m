function Dat = H_WLC_realdatahom(chromosome);

addpath(genpath('./'))
addpath(genpath('../Pseudotime_Algorithms'))

warning off all

%Load existing data already generated
[Xreal1,Yreal1,Zreal1] = textread('./YeastGenome/ChromTraj.txt','%f %f %f');

%load(['./HomogeneousRerun_Sept2012B/HomogenousPolymer_' num2str(chromosome) '_FullTraj.mat'])
%M = median(log(Dat.MCMC.Y(1000:100:end)));
%S = std(log(Dat.MCMC.Y(1000:100:end)));
%V = var(log(Dat.MCMC.Y(1000:100:end)));
%lens1 = length((Dat.MCMC.Y(1000:100:end)));

%index of individual chromosomes
Indexes{1}  = 1:450; 
Indexes{2}  = 450+1:2059;
Indexes{3}  = 2059+1:2669;
Indexes{4}  = 2669+1:5719;
Indexes{5}  = 5719+1:6849;
Indexes{6}  = 6849+1:7379;
Indexes{7}  = 7379+1:9549;
Indexes{8}  = 9549+1:10659;
Indexes{9}  = 10659+1:11509;
Indexes{10} = 11509+1:12979;
Indexes{11} = 12979+1:14289;
Indexes{12} = 14289+1:19119;
Indexes{13} = 19119+1:20949;
Indexes{14} = 20949+1:22499;
Indexes{15} = 22499+1:24669;
Indexes{16} = 24669+1:26538;

%__________________________________________________________________________
%X, Y and Z coordinates of chromsomes
XXX = Xreal1(Indexes{chromosome});
YYY = Yreal1(Indexes{chromosome});
ZZZ = Zreal1(Indexes{chromosome});
stepsize = 1;%double(int64(length(XXX)/100));
Xtrue = XXX(5:stepsize:end-5);
Ytrue = YYY(5:stepsize:end-5);
Ztrue = ZZZ(5:stepsize:end-5);

%Rescale polymer to unit length
dX    = diff(Xtrue);
dY    = diff(Ytrue);
dZ    = diff(Ztrue);
Delta = sqrt(dX.^2 + dY.^2 + dZ.^2);
dX    = dX./Delta;
dY    = dY./Delta;
dZ    = dZ./Delta;
Xtrue = [Xtrue(1,1),Xtrue(1,1)+cumsum(dX')];%/length(Xtrue);
Ytrue = [Ytrue(1,1),Ytrue(1,1)+cumsum(dY')];%/length(Xtrue);
Ztrue = [Ztrue(1,1),Ztrue(1,1)+cumsum(dZ')];%/length(Xtrue);

%load('./Heterologous WLC/CorrectPriorTuneN/Run2/CorrectPrior_Tunep_3.mat')

%__________________________________________________________________________
%Parameters for the polymer statistic model
Nochains = 4;
No       = 200000;                 %Number of steps in the Polymer Markov Chain.
No2      = 300000;                 %Number steps in Marko Chain
len      = length(Xtrue);                    %Number of beads in the polymer chain.
kappa    = 1e-10;                  %Von Mises Fisher clustering of telomeres
Radius   = 1600;                   %Radius of the Nuclear Periphery (um).
U        = 0;                      %Potential for beads outside of the NP.
a        = 1;%/(99);                %Bond vector length.
%__________________________________________________________________________
%Parameters for the GP model
covfunc = {'covSum',{'covSEiso','covNoise'}};  %Covariance function
N       = len-2;
%Htrue       = [log(.1); log(0.3); log(1e-5)]; %Hyperparameters of covariance function
%H       = [log(.05); log(0.5); log(1e-5)]; %Hyperparameters of covariance function
x       = linspace(0,1,N);
%K       = feval(covfunc{:}, H, x'); %Covariance matrix at observation points
n       = length(x);
sig     = 0.2;

Dat.Polymer.X  = Xtrue;
Dat.Polymer.Y  = Ytrue;
Dat.Polymer.Z  = Ztrue;
    
for nochains = 1%:Nochains
    rng(nochains)
    
%__________________________________________________________________________
%Create arrays/matrices for storing data
p    = 0.1;
Y    = zeros(No2,1);
Po    = zeros(No2,1);
%GPL    = zeros(No2,1);
Zan    = zeros(No2,1);
PolyL = zeros(No2,1);
%Phyp = zeros(No2,1);
%Hmcmc = zeros(No2,3);
%Hmcmc(1,:) = H';
%NoSwitches = zeros(No2,1);

%__________________________________________________________________________
%Make intial guess in MCMC chain
Y(1,1)     = randn(1,1);%real(gsamp(zeros(1,N)',K,1)); %Generate sample function represent true kappa
k          = exp(Y(1,1) -1); %set current kappa to Yguess
Parameters = [No,len,Radius,U,a,k];  %Group parameters.

%Estimate partition function
Zan(1,1)   = EstimatePartition(Parameters);


for j = 1%:size(PInd,2)
Delta      = sqrt(diff(Xtrue(j,:)).^2 + diff(Ytrue(j,:)).^2 + diff(Ztrue(j,:)).^2);
E1         = [(diff(Xtrue(j,:))./Delta)',(diff(Ytrue(j,:))./Delta)',(diff(Ztrue(j,:))./Delta)'];
En         = sum(k'.*diag(E1*E1',1));
PolyLA(j,1) = (1/a)*En;
end
PolyL(1,1) = sum(PolyLA);

%Phyp(1,1) = loggauss([-1,-1],diag(ones(1,2))*0.5,Hmcmc(1,1:2));


%K = feval('covSEiso',H(1:2),x'); 
%L = chol(K/exp(2*H(3))+eye(n));
%alpha   = solve_chol(L,Y(1,:)')/exp(2*H(3));
%nlZ     = (Y(1,:)')'*alpha/2 + sum(log(diag(L))) + n*log(2*pi*exp(2*H(3)))/2;  %-log marg lik

%GPL(j,1)   = -nlZ;%loggauss(zeros(1,N),K+diag(ones(1,N))*1e-5,Y(1,:));

%NoSwitches(1,1) = p;

Po(1,1) = PolyL(1,1) - N*Zan(1,1);

for i = 2:No2
   
        %Print output
        if double(int64(i/10000))==(i/10000)
        Accep1 = (length(unique(Po(i-9999:i-1,1)))-1)./(length(Po(i-9999:i-1,1))-1);
        disp(['Step ' num2str(i) '. Acceptence rate = ' num2str(Accep1) '. No = ' num2str(No)])
        end
    
    if double(int64(i/100))==(i/100)
    %disp(['Step ' num2str(i)])    
    %update p
    Accep = (length(unique(Po(i-99:i-1,1)))-1)./(length(Po(i-99:i-1,1))-1);    
        if abs(0.25-Accep)<0.1          
        else
            if Accep>0.25
            p = p*0.9;
            p = max(p,0.05);
            else
            p = p*1.1;
            p = min(p,0.95);
            end
        end    
    else
    end
    
    %Tune hyperparams
     if double(int64(i/10))==(i/10)            
       Accep = length(unique(Po(i-9:i-1,1)))./length(Po(i-9:i-1,1));    
         if abs(0.25-Accep)<0.1          
         else
             if abs(0.2-Accep)<0.25
               sig = sig*0.9;
                else
               sig = sig*1.1;
             end
               sig = max(sig,0.1);
               sig = min(sig,0.5);
         end        
     end
        
    
    Y(i,1) = Y(i-1,1)+randn(1,1)*sig;    
    k      = exp(Y(i,1) -1); %set current kappa to Yguess
            
    Parameters = [No,len,Radius,U,a,k];  %Group parameters.

    %Estimate partition function
    Zan(i,1)   = EstimatePartition(Parameters);


for j = 1%size(PInd,2)
Delta      = sqrt(diff(Xtrue(j,:)).^2 + diff(Ytrue(j,:)).^2 + diff(Ztrue(j,:)).^2);
E1         = [(diff(Xtrue(j,:))./Delta)',(diff(Ytrue(j,:))./Delta)',(diff(Ztrue(j,:))./Delta)'];
En         = sum(k'.*diag(E1*E1',1));
PolyLA(j,1) = (1/a)*En;
end
PolyL(i,1) = sum(PolyLA);
    
Po(i,1) = PolyL(i,1) - N*Zan(i,1);%*length(PInd);
deltaE = abs(Po(i,1) - Po(i-1,1));

    if (Po(i,1))>(Po(i-1,1)) || rand(1,1)<(exp(-deltaE))
    
    else
        Y(i,:) = Y(i-1,:);
        Zan(i,1) = Zan(i-1,1);
        Po(i,1) = Po(i-1,1);
        PolyL(i,1) = PolyL(i-1,1);
        %GPL(i,1) = GPL(i-1,1);
        %NoSwitches(i,1) = NoSwitches(i-1,1);
       %Hmcmc(i,:) =  Hmcmc(i-1,:);
       %H = Hmcmc(i-1,:);
       %Phyp(i,1) = Phyp(i-1,1);
    end
end


Posterior{nochains} = Po;
%GPLike{nochains} = GPL;
PolymerLike{nochains} = PolyL;
Yinferred{nochains} = Y;
Partition{nochains} = Zan;
%SwitchProb{nochains} = NoSwitches;

end



Dat.MCMC.Po    = Posterior;
%Dat.MCMC.GPL   = GPLike;
Dat.MCMC.PolyL = PolymerLike;
Dat.MCMC.Y     = Yinferred;
Dat.MCMC.Z     = Partition;
%Dat.MCMC.p     = SwitchProb;

%Dat.True.Z     = Zantrue;
%Dat.True.GPL   = GPLtrue;
%Dat.True.PolyL = PolyLtrue;
%Dat.True.Ylat  = Ylat;

%Dat.Polymer.X  = Xtrue;
%Dat.Polymer.Y  = Ytrue;
%Dat.Polymer.Z  = Ztrue;

mkdir('./HWLC_r')
save(['./HWLC_r/InCorrectPrior_Tunep_Tunehyp_' num2str(chromosome) '_hom.mat'],'Dat')

%keyboard
%Ind0 = linspace(1,100000,100000);
%Ind1 = Ind0(find(Po==max(Po)));
%Ind1 = Ind1(1,1);

%plot((Y(1,:)),'g.')
%hold on
%plot(Ylat,'r.')
%plot((Y(Ind1,:)),'b.')
%errorbar(mean(Y(15000:5:end,:),1),2*sqrt(var(Y(15000:5:end,:))),'k')
