%Define Parameters.
No      = 200000;
len     = 100;
Radius  = 1600;
U       = 10000;
kappa   = 100;
k       = (0.5)*230000*0.34/10;
a       = 230000*0.34/(10*100);
Cen     = 1;

%Vector containing all parameters.
Parameters1     = [No,len,k,Radius,U,a,Cen]; VMF = [0,0,0,0];

%Generate an example wormlike chain.
[x0,y0,z0,theta,phi] = tMultivMF(0,0,1,kappa,Radius,1);
[x1,y1,z1]      = RandFlightCon([x0,y0,z0],[x0,y0,z0],a,len);
[X,Y,Z,Energy]  = SingleTether(Parameters1,x1,y1,z1);
%Start and end positions.
XB = X(100000,1); 
XE = X(100000,100);
YB = Y(100000,1);
YE = Y(100000,100);
ZB = Z(100000,1); 
ZE = Z(100000,100);

clear X Y Z Energy

%Array holding probabilities.
Prob = zeros(1,10);

for i = 1:10

Parameters2     = [No,len,(1/i)*230000*0.34/10,Radius,U,a,Cen];
    
%Calculate the Partition function for a WLC fixed at one end.    
[X,Y,Z,Energy]  = SingleTether(Parameters2,x1,y1,z1);
subsample       = 100000:50:200000;
Z1              = Energy(subsample,1); 
Partition       = sum(exp(Z1 - max(Z1)))/length(Z1)
figure(i)
plot(Energy)
clear X Y Z Energy

%Calculate initial conformation satisfying the end positions of the originally generated WLC.
[Pos]   = BBRandFlight(100,a,[XB,YB,ZB],[XE,YE,ZE]);
X0x     = Pos(:,1)';
Y0x     = Pos(:,2)';
Z0x     = Pos(:,3)';

%Calculate unnormalised probability of the such conformation under
%certain parameters.
[X,Y,Z,Energy]  = FixedEnd(Parameters2,X0x,Y0x,Z0x);
Z0              = Energy; 

hold on
plot(Energy,'r')

%Normalise probability.
Prob(1,i)            = sum((exp(Z0 - max(Z1))/Partition))/length(Z0);

clear X Y Z Energy
disp([num2str(i)])
end

%plot(Prob1,'o-')