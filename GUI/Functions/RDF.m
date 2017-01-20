function [val] = RDF(net,Loci);

X = net.X;
Y = net.Y;
Z = net.Z;
Parameters = net.Parameters;

%Calculate the raidal distribution for a set of loci and approximate
%according to a Beta distribution

Radius = Parameters(1,4); % Radius of NP.

for i = 1:length(Loci)
   
    R = sqrt(X(:,Loci(i,1)).^2 + Y(:,Loci(i,1)).^2 + Z(:,Loci(i,1)).^2)/Radius;
    try
    [val{i}] = betafit(R);
    catch
        warning('Some values greater than 1')
       R = R(find(R<1)); 
   [val{i}] = betafit(R);
    end
end
