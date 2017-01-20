function [val] = LocusLocusDistance(varargin)

switch nargin
    
    case 0
        warning('Wrong number of input parameters')
    case 1
        warning('Wrong number of input parameters')    
    case 2
        warning('Wrong number of input parameters')    
    case 3
        net1  = varargin{1};
        LociA = varargin{2};
        LociB = varargin{3};
        Parameters = net1.Parameters;
        Radius = Parameters(1,4);
        X1 = net1.X;
        X2 = net1.X;
        Y1 = net1.Y;
        Y2 = net1.Y;                
        Z1 = net1.Z;
        Z2 = net1.Z;        
    case 4
        net1  = varargin{1};
        net2  = varargin{2};
        LociA = varargin{3};
        LociB = varargin{4};
        Parameters = net1.Parameters;
        Radius = Parameters(1,4);
        X1 = net1.X;
        X2 = net2.X;
        Y1 = net1.Y;
        Y2 = net2.Y;                
        Z1 = net1.Z;
        Z2 = net2.Z;
end
        
       

%Calculate the distance between two loci and approximate accoring to a beta
%distribution



for i = 1:length(LociA)
   
    R = sqrt((X1(:,LociA(i,1)) - X2(:,LociB(i,1))).^2 + ... 
    (Y1(:,LociA(i,1)) - Y2(:,LociB(i,1))).^2 + ... 
    (Z1(:,LociA(i,1)) - Z2(:,LociB(i,1))).^2)/(2*Radius);    

    try
        [val{i}] = betafit(R);
    catch
        warning('Some values greater than 1 or less than 0')
        R        = R(find(R<1)); 
        R        = R(find(R>0)); 
        [val{i}] = betafit(R);
    end
end
