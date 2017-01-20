function [D] = CDM(net);
%Function to genetate the chromosome distance map given two sets of trajectories.
%Note that if only a single set of trajectories is given, the chromosomes
%will be treated as homologues and the dataset split in two for calculations.

%switch nargin
%case 1
%
%    net  = varargin{1};    
                
    warning('Only one set of trajectories given. Treating as homologoues and partitioning data')
                
        X1 = []; Y1 = []; Z1 = [];
        for i = 1:length(net)        
            X=[X1;net{i}.X];
            Y=[Y1;net{i}.Y];
            Z=[Z1;net{i}.Z];
        end
         
        %Randomly permute the net
        Ind1 = randperm(size(X,1));
        Ind2 = randperm(size(X,1));        
        %Take one half of the net
 
        X1 = X(Ind1,:);
        Y1 = Y(Ind1,:);
        Z1 = Z(Ind1,:);    
        X2 = X(Ind2,:);
        Y2 = Y(Ind2,:);
        Z2 = Z(Ind2,:);              
        
        Parameters1 = net{1}.Parameters;
        R1 = Parameters1(1,4);
        R2 = R1;
        
        BS = 2;
        
%case 2
%    data1 = varargin{1};    
%    data2 = varargin{2}; 
%     
%         X1 = []; Y1 = []; Z1 = [];
%         for i = 1:length(data1)        
%             X1=[X1;data1{i}.X];
%             Y1=[Y1;data1{i}.Y];
%             Z1=[Z1;data1{i}.Z];
%         end
%         X2 = []; Y2 = []; Z2 = [];
%         for i = 1:length(data2) 
%             X2=[X2;data2{i}.X];
%             Y2=[Y2;data2{i}.Y];
%             Z2=[Z2;data2{i}.Z];
%         end
%         
%         if size(X1,1)~=size(X2,1)
%                 warning('Different number of samples for the two architectures')                
%             l = min(size(X1,1),size(X2,1));
%         else
%             l = size(X1,1);
%         end
%         
%         Parameters1 = data1{1}.Parameters;
%         Parameters2 = data2{1}.Parameters;
%         R1 = Parameters1(1,4);
%         R2 = Parameters2(1,4);
%         
%         BS = 2;
%     
% case 3
% 
%         data1 = varargin{1};    
%         data2 = varargin{2}; 
%     
%         X1 = []; Y1 = []; Z1 = [];
%         for i = 1:length(data1)        
%             X1=[X1;data1{i}.X];
%             Y1=[Y1;data1{i}.Y];
%             Z1=[Z1;data1{i}.Z];
%         end
%         X2 = []; Y2 = []; Z2 = [];
%         for i = 1:length(data2) 
%             X2=[X2;data2{i}.X];
%             Y2=[Y2;data2{i}.Y];
%             Z2=[Z2;data2{i}.Z];
%         end
%         
%         if size(X1,1)~=size(X2,1)
%                 warning('Different number of samples for the two architectures')                
%             l = min(size(X1,1),size(X2,1));
%         else
%             l = size(X1,1);
%         end
%         
%         Parameters1 = data1{1}.Parameters;
%         Parameters2 = data2{1}.Parameters;
%         R1 = Parameters1(1,4);
%         R2 = Parameters2(1,4);
%         
%         BS = varargin{3};
% 
% end

if R1==R2
    R = R1;
else
    warning('Radii are different. Not normalising by radius')
    R = 1;
end



l1 = size(X1,2);
l2 = size(X2,2);

D = zeros(l1,l2,BS);

for BtStrap = 1:BS

    Ind1 = randperm(l1);
    Ind2 = randperm(l2);
    
    X01 = X1(Ind1,:);
    Y01 = Y1(Ind1,:);
    Z01 = Z1(Ind1,:);
    X02 = X2(Ind2,:);
    Y02 = Y2(Ind2,:);
    Z02 = Z2(Ind2,:);
    
for i = 1:l1
    for j = 1:l2
        
        D(i,j,BtStrap) = mean(sqrt((X01(:,i)-X02(:,j)).^2 + ...
                 (Y01(:,i)-Y02(:,j)).^2 + ...
                 (Z01(:,i)-Z02(:,j)).^2)./(2*R));
             
% V(i,j) = std(sqrt((X1(:,i)-X2(:,j)).^2 + ...
%                 (Y1(:,i)-Y2(:,j)).^2 + ...
%                 (Z1(:,i)-Z2(:,j)).^2)./(2*R));             
             
     end
end


end

D = (mean(D,3));

