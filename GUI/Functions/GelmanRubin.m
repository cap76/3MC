function [GR] = GelmanRubin(net,net1);


I  = length(net.X(1,:));
l  = length(net.X(:,1));
l2 = double(int64(l/2));
m  = 2;

for i = 1:I    

Xhat  = mean(net.X(l2:end,i));
Xhat1 = mean(net1.X(l2:end,i));
n     = linspace(1,length(net.X)-l2,length(net.X)-l2);


Wx  = (sum((net.X(l2:end,i) - Xhat).^2) + sum((net1.X(l2:end,i) - Xhat1).^2))./((l-l2-1)*m);
Bx  = ((l-l2)/(1)) * (    (mean(net.X(l2:end,i))-mean([net.X(l2:end,i);net1.X(l2:end,i)])).^2  +   (mean(net1.X(l2:end,i))-mean([net.X(l2:end,i);net1.X(l2:end,i)])).^2    ) ;
V   = (1 - 1/(l-l2))*Wx + Bx/(l-l2);

GR(i,1) = sqrt(V./Wx);
end