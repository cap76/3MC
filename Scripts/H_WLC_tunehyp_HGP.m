function Dat = H_WLC_tunehyp_HGP(Arc);

%load('/Users/christopherpenfold/Documents/MATLAB/3MC/HWLC/CorrectPrior_Tunep_3.mat')

addpath(genpath('./'))
addpath(genpath('../GitHub'))
addpath(genpath('../Pseudotime_Algorithms'))

%Kstar = 
%ell     = exp(hyp.cov(1));    %Characteristic length scale
%sf2     = exp(2*hyp.cov(2));  %Signal variance
%sn2     = exp(2*hyp.lik);      %Noise variance of likGauss
%n       = length(x);
%K       = exp(2*H(2))*exp(- Kstar/exp(H(1))^2 /2);

%Load existing data already generated
%load('./Heterologous WLC/CorrectPriorTuneN/Run2/CorrectPrior_Tunep_3.mat')

addpath(genpath('./'))
addpath(genpath('../Pseudotime_Algorithms'))

warning off all

indvals = [1,2,4,8,16]; %No. chromosome trajectories
%__________________________________________________________________________
%Parameters for the polymer statistic model
Nochains = 4;
No       = 200000;                 %Number of steps in the Polymer Markov Chain.
No2      = 300000;                 %Number steps in Marko Chain
len      = 100;                    %Number of beads in the polymer chain.
kappa    = 1e-10;                  %Von Mises Fisher clustering of telomeres
Radius   = 1600;                   %Radius of the Nuclear Periphery (um).
U        = 0;                      %Potential for beads outside of the NP.
a        = 1/(99);                %Bond vector length.
%__________________________________________________________________________
%Parameters for the GP model
covfunc = {'covSum',{'covBranchingProcess_3A','covNoise'}};  %Covariance function
N       = len-2;

try
    load('HWLC_HGP/Run1.mat')
    Ylat = Output.Ylat;
    Xtrue = Output.Xtrue;
    Ytrue = Output.Ytrue;
    Ztrue = Output.Ztrue;
    Htrue = Output.Htrue;
    H = Output.H;
    x       = [linspace(0,1,N)',ones(1,N)';linspace(0,1,N)',2*ones(1,N)';linspace(0,1,N)',3*ones(1,N)']';
    K       = feval(covfunc{:}, Htrue, x'); %Covariance matrix at observation points
    Kalt       = feval(covfunc{:}, H, x'); %Covariance matrix at observation points
    n       = length(x);
    sig     = 0.2;
catch
    Htrue  = [-10000;1;-10000;1; log(.1); log(0.3); log(.1); log(0.3); log(.1); log(0.3); log(1e-5)]
    H      = [-10000;1;-10000;1; log(.05); log(0.5); log(.05); log(0.5); log(.05); log(0.5); log(1e-5)]; %Hyperparameters of covariance function
x       = [linspace(0,1,N)',ones(1,N)';linspace(0,1,N)',2*ones(1,N)';linspace(0,1,N)',3*ones(1,N)']';
K       = feval(covfunc{:}, Htrue, x'); %Covariance matrix at observation points
Kalt       = feval(covfunc{:}, H, x'); %Covariance matrix at observation points
n       = length(x);
sig     = 0.2;
%Generate the latent GP
Ylat    = real(gsamp(zeros(1,3*N)',K,1)); %Generate sample function represent true kappa

%__________________________________________________________________________
%Remaining polymer statistic model that depend on GP model
for reploop = 1:3 %3 hierachical GPs
ktrue           = exp(Ylat((reploop-1)*N+1:reploop*N) - 1);                      %Persistence length (um).
Parameters  = [No,len,Radius,U,a,ktrue];  %Group parameters.
Zantrue(1,reploop)   = EstimatePartition(Parameters);
%__________________________________________________________________________
%Simulate initial configuration
[x0,y0,z0,theta,phi] = tMultivMF(0,0,1,kappa,Radius,1);
[x1,y1,z1]           = RandFlightCon([x0,y0,z0],[x0,y0,z0],a,len);
%__________________________________________________________________________
%Generate WLC trajectories
[X,Y,Z,Energy]       = FreeWLCTV(Parameters,x1,y1,z1);
%Save data for later use
PInd  = 1000:50:(size(X,1));
Xtrue{reploop} = X(PInd,:);
Ytrue{reploop} = Y(PInd,:);
Ztrue{reploop} = Z(PInd,:);
end
    Output.Xtrue = Xtrue;
    Output.Ytrue = Ytrue;
    Output.Ztrue = Ztrue;
    Output.Htrue = Htrue;
    Output.H = H;
    Output.Ylat = Ylat;
    save('HWLC_HGP/Run1.mat','Output')
end

%%Etrue = mean(Energy(1000:end));
%Ylattrue = Dat.True.Ylat;
%Hyptrue = Dat.GP.H;

clear Dat

Dat.True.Ylat = Ylat;
Dat.True.H = Htrue;

Dat.Polymer.X  = Xtrue;
Dat.Polymer.Y  = Ytrue;
Dat.Polymer.Z  = Ztrue;


for jki = Arc%length(indvals)

    
for nochains = 1%:Nochains
    rng(nochains)
    
%__________________________________________________________________________
%Create arrays/matrices for storing data
p    = 0.1;
Y    = zeros(No2,3*N);
Po    = zeros(No2,1);
GPL    = zeros(No2,1);
Zan    = zeros(No2,3);
PolyL = zeros(No2,3);
Phyp = zeros(No2,1);
Hmcmc = zeros(No2,11);
Hmcmc(1,:) = H';
NoSwitches = zeros(No2,1);

%__________________________________________________________________________
%Make intial guess in MCMC chain
Y(1,:)     = real(gsamp(zeros(1,3*N)',K,1)); %Generate sample function represent true kappa

for reploop = 1:3
k          = exp(Y(1,(reploop-1)*N+1:reploop*N) -1); %set current kappa to Yguess
Parameters = [No,len,Radius,U,a,k];  %Group parameters.
%Estimate partition function
Zan(1,reploop)   = EstimatePartition(Parameters);
for j = 1:indvals(1,jki)%:size(PInd,2)
Delta      = sqrt(diff(Xtrue{reploop}(j,:)).^2 + diff(Ytrue{reploop}(j,:)).^2 + diff(Ztrue{reploop}(j,:)).^2);
E1         = [(diff(Xtrue{reploop}(j,:))./Delta)',(diff(Ytrue{reploop}(j,:))./Delta)',(diff(Ztrue{reploop}(j,:))./Delta)'];
En         = sum(k'.*diag(E1*E1',1));
PolyLA(j,1) = (1/a)*En;
end
PolyL(1,reploop) = sum(PolyLA);
end

K = feval('covBranchingProcess_3A',H(1:end-1),x'); 
L = chol(K/exp(2*H(end))+eye(n));
alpha   = solve_chol(L,Y(1,:)')/exp(2*H(end));
nlZ     = (Y(1,:)')'*alpha/2 + sum(log(diag(L))) + n*log(2*pi*exp(2*H(end)))/2;  %-log marg lik

GPL(1,1)   = -nlZ;%loggauss(zeros(1,N),K+diag(ones(1,N))*1e-5,Y(1,:));
NoSwitches(1,1) = p;


Phyp(1,1) = loggauss([-1,-1,-1,-1,-1,-1],diag(ones(1,6))*0.5,Hmcmc(1,5:end-1));

Po(1,1) = sum(PolyL(1,:)) + GPL(1,1) + Phyp(1,1) - sum(indvals(1,jki)*Zan(1,:));%*length(PInd)

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
            Hmcmc(i,:) = Hmcmc(i-1,:);
            Hmcmc(i,5:end-1) = Hmcmc(i,5:end-1) + sig*randn(1,6);  
        else
            if abs(0.2-Accep)<0.25
              sig = sig*0.9;
               else
              sig = sig*1.1;
            end
              sig = max(sig,0.1);
              sig = min(sig,0.5);
             Hmcmc(i,:) = Hmcmc(i-1,:);
             Hmcmc(i,5:end-1) = Hmcmc(i,5:end-1) + sig*randn(1,6);   
        end        
        else
             Hmcmc(i,:) = Hmcmc(i-1,:);
    end
    H = Hmcmc(i,:)';
        
    
    %Update Y
    No = binornd(N-2,p)+1;
    NoSwitches(i,1) = p;
    Ind = 1:1:N;
    Ind2= randperm(N);
    Ind2= Ind2(1:No);
    Ind2=sort(Ind2); 
    
    uchoice = randi([1 3],1,1);
    Ind2 = (Ind2-1) + (uchoice-1)*N+1;
    
    set2 = 1:1:length(x);
    set2 = setdiff(set2,(uchoice-1)*N+1:uchoice*N);
    Ind2 = [Ind2,set2];
    
    xtrain = x(Ind2);
    ytrain = Y(i-1,Ind2);
    Indc = setdiff(1:1:length(x),Ind2);

    %Update the GP at subsampled locations.
    K1     = feval(covfunc{:}, H, x');    
    L1     = chol(K1(Ind2,Ind2))';                        % cholesky factorization of the covariance
    alpha1 = solve_chol1(L1',ytrain');
    Kstar  = K1(Ind2,Indc);    
    Mu     = Kstar' * alpha1;                                      % predicted means
    alphaprime  = solve_chol(L1',Kstar);
    VAR    = K1(Indc,Indc) - Kstar' * alphaprime; 
    Y(i,:) = Y(i-1,:);
    Y(i,Indc) = real(gsamp(Mu,VAR,1));
    k          = exp(Y(i,:) -1); %set current kappa to Yguess
        
    %Get the FT, abd BT
    FT = loggauss(Mu,VAR,Y(i,Indc));
    BT = loggauss(Mu,VAR,Y(i-1,Indc));
     
    %Get the ML of the full GP   
    try
    K = feval('covBranchingProcess_3A',H(1:end-1),x'); 
    L = chol(K/exp(2*H(end))+eye(n));
    alpha   = solve_chol(L,Y(i,:)')/exp(2*H(end));
    GPL(i,1)     = - ( (Y(i,:)')'*alpha/2 + sum(log(diag(L))) + n*log(2*pi*exp(2*H(end)))/2 );  %-log marg lik
    catch
        keyboard
    end
    
   
    %Estimate partition function
    for reploop = 1:3
    k          = exp(Y(i,(reploop-1)*N+1:reploop*N) -1); %set current kappa to Yguess
    Parameters = [No,len,Radius,U,a,k];  %Group parameters.
    Zan(i,reploop)   = EstimatePartition(Parameters);

    for j = 1:indvals(1,jki)%size(PInd,2)
    Delta      = sqrt(diff(Xtrue{reploop}(j,:)).^2 + diff(Ytrue{reploop}(j,:)).^2 + diff(Ztrue{reploop}(j,:)).^2);
    E1         = [(diff(Xtrue{reploop}(j,:))./Delta)',(diff(Ytrue{reploop}(j,:))./Delta)',(diff(Ztrue{reploop}(j,:))./Delta)'];
    En         = sum(k'.*diag(E1*E1',1));
    PolyLA(j,1) = (1/a)*En;
    end
    PolyL(i,reploop) = sum(PolyLA);
    end

    

    Phyp(i,1) = loggauss([-1,-1,-1,-1,-1,-1],diag(ones(1,6))*0.5,Hmcmc(i,5:end-1));
    %Phyp(i,1) = loggauss([-1,-1],diag(ones(1,2))*0.5,Hmcmc(i,1:2));
    %Po(i,1) = PolyL(i,1) + GPL(i,1) + Phyp(i,1) - indvals(1,jki)*Zan(i,1);%*length(PInd);

    Po(i,1) = sum(PolyL(i,:)) + GPL(i,1) + Phyp(i,1) - sum(indvals(1,jki)*Zan(i,:));%*length(PInd)
            
    deltaE = abs(Po(i,1) + (BT-FT) - Po(i-1,1));
    
  
    %keyboard
    if (Po(i,1)+BT)>(Po(i-1,1)+FT) || rand(1,1)<(exp(-deltaE))
    
    else
        Y(i,:) = Y(i-1,:);
        Zan(i,:) = Zan(i-1,:);
        Po(i,1) = Po(i-1,1);
        PolyL(i,:) = PolyL(i-1,:);
        GPL(i,1) = GPL(i-1,1);
        NoSwitches(i,1) = NoSwitches(i-1,1);
       Hmcmc(i,:) =  Hmcmc(i-1,:);
       H = Hmcmc(i-1,:);
       Phyp(i,1) = Phyp(i-1,1);
    end
end


Posterior{nochains} = Po;
GPLike{nochains} = GPL;
PolymerLike{nochains} = PolyL;
Yinferred{nochains} = Y;
Partition{nochains} = Zan;
SwitchProb{nochains} = NoSwitches;

end



Dat.MCMC.Po    = Posterior;
Dat.MCMC.GPL   = GPLike;
Dat.MCMC.PolyL = PolymerLike;
Dat.MCMC.Y     = Yinferred;
Dat.MCMC.Z     = Partition;
Dat.MCMC.p     = SwitchProb;

%Dat.True.Z     = Zantrue;
%Dat.True.GPL   = GPLtrue;
%Dat.True.PolyL = PolyLtrue;
%Dat.True.Ylat  = Ylat;

Dat.Polymer.X  = Xtrue;
Dat.Polymer.Y  = Ytrue;
Dat.Polymer.Z  = Ztrue;

%mkdir('./HWLC_r')
save(['./HWLC_HGP//InCorrectPrior_Tunep_Tunehyp_HGP_' num2str(jki) '.mat'],'Dat')

end

%keyboard
%Ind0 = linspace(1,100000,100000);
%Ind1 = Ind0(find(Po==max(Po)));
%Ind1 = Ind1(1,1);

%plot((Y(1,:)),'g.')
%hold on
%plot(Ylat,'r.')
%plot((Y(Ind1,:)),'b.')
%errorbar(mean(Y(15000:5:end,:),1),2*sqrt(var(Y(15000:5:end,:))),'k')
