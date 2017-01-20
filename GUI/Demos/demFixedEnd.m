%Example script for plotting a chromosome configuration with both ends
%tethered.

if ~isdeployed
    addpath(genapth('../'));
end

No      = 200000;
len     = 100; %Number of links in chain.
a       = 4;   %Link size.
k       = 400;
Radius  = 150;
U       = 10000;
%Define start and end positions. Note: the distance between these two
%vectors must not be greater than the length of the chain.
x0      = [-150,0,0];
xn      = [150,0,0];
Phi     = 30;

%Generate initial conformation.
[Pos] = BBRandFlight(len,a,x0,xn);

%Generate a random flight conformation as BBRandFlight is less than ideal.
Parameters1 = [No,len,1,Radius,U,a];
X0x = Pos(:,1)';
Y0x = Pos(:,2)';
Z0x = Pos(:,3)';
[net] = FixedEnd(Parameters1,X0x,Y0x,Z0x);

%Change the rigigidity parameter to the one we're interested in. Note above
%we set it to a very small value in order to force the initial conformation
%to correspond to a random flight state. Now we are interested in
%exploreing chromosomes with persistence length = k.
Parameters  = [No,len,k,Radius,U,a];
X0x = net.X(length(net.X),:);
Y0x = net.Y(length(net.X),:);
Z0x = net.Z(length(net.X),:);

[net] = FixedEnd(Parameters,X0x,Y0x,Z0x);

%Plot example path from MCMC chain.

Width   = -0.0056*Phi^2 + 1.6*Phi -3.4; 

%Plot example configuration after N steps.
PolymerPic(net.Parameters,net.X,net.Y,net.Z,Width,No);


%Fit a beta distribution to a set of loci of choice, in this case loci 3,
%10 and 15

Loci = [3,10,15]';

[val] = RDF(net.Parameters,net.X,net.Y,net.Z,Loci);
x     = linspace(0,1,100);
figure(2)
for j = 1:size(val,2)
    Y = betapdf(x,val{j}(1,1),val{j}(1,2));
    hold on
    plot(x,Y,'-')
end