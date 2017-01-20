function [D] = LDM(net)

    
    
    Parameters1 = net{1}.Parameters;
    radius      = Parameters1(1,4);    
    
    X1 = [];
    Y1 = [];
    Z1 = [];
    
    %Combine different seeds
    for i = 1:length(net);
        X1=[X1;net{i}.X];
        Y1=[Y1;net{i}.Y];
        Z1=[Z1;net{i}.Z];
    end

l1 = size(X1,2);

D = zeros(l1,l1);
    
for i = 1:l1
    for j = 1:l1
        
 %Calculate the mean distance between all locus-locus distances       
 D(i,j) = mean(sqrt((X1(:,i)-X1(:,j)).^2 + ...
                 (Y1(:,i)-Y1(:,j)).^2 + ...
                 (Z1(:,i)-Z1(:,j)).^2)./(2*radius));
 
 %Calculate the variance between all locus-locus distances            
 %V(i,j) = std(sqrt((X1(:,i)-X1(:,j)).^2 + ...
 %                (Y1(:,i)-Y1(:,j)).^2 + ...
 %                (Z1(:,i)-Z1(:,j)).^2)./(2*R));             
             
    end
end