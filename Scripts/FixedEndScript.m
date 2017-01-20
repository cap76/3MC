%Example script for plotting a chromosome configuration with both ends
%tethered.
addpath('./Functions','./netlab')


No      = 20000;
len     = 100; %Number of links in chain.
a       = 4;   %Link size.
k       = 400;
Radius  = 150;
U       = 10000;
%Define start and end positions. Note: the distance between these two
%vectors must not be greater than the length of the chain.
x0      = [-150,0,0];
xn      = [150,0,0];

%Generate initial conformation.
[Pos] = BBRandFlight(len,a,x0,xn);

%Generate a random flight conformation as BBRandFlight is less than ideal.
Parameters1 = [No,len,1,Radius,U,a];
X0x = Pos(:,1)';
Y0x = Pos(:,2)';
Z0x = Pos(:,3)';
[X1,Y1,Z1,Energy] = FixedEnd(Parameters1,X0x,Y0x,Z0x);

%Change the rigigidity parameter to the one we're interested in. Note above
%we set it to a very small value in order to force the initial conformation
%to correspond to a random flight state. Now we are interested in
%exploreing chromosomes with persistence length = k.
Parameters  = [No,len,k,Radius,U,a];
X0x = X1(length(X1),:);
Y0x = Y1(length(X1),:);
Z0x = Z1(length(X1),:);

[X,Y,Z,Energy] = FixedEnd(Parameters,X0x,Y0x,Z0x);

%Plot example path from MCMC chain.
PolymerPic(X,Y,Z,Radius,No);