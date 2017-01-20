function [D,V] = HDM(net1,net2);


X1=net1.X;
Y1=net1.Y;
Z1=net1.Z;

Parameters1 = net1.Parameters;


X2=net2.X;
Y2=net2.Y;
Z2=net2.Z;

Parameters2 = net2.Parameters;

R1 = Parameters1(1,4);
R2 = Parameters1(1,4);

if R1==R2
    R = R1;
else
    warning('Radii are different. Not normalising by radius')
    R = 1;
end

l1 = size(X1,2);
l2 = size(X2,2);

for i = 1:l1
    for j = 1:l2
        
        D(i,j) = mean(sqrt((X1(:,i)-X2(:,j)).^2 + ...
                 (Y1(:,i)-Y2(:,j)).^2 + ...
                 (Z1(:,i)-Z2(:,j)).^2)./(2*R));
             
 V(i,j) = std(sqrt((X1(:,i)-X2(:,j)).^2 + ...
                 (Y1(:,i)-Y2(:,j)).^2 + ...
                 (Z1(:,i)-Z2(:,j)).^2)./(2*R));             
             
     end
end
