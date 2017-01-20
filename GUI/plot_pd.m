function plot_pd(net)

%plot(net{1}.E);
%hold on
%plot(net{2}.E,'r');
%[D] = Recombination(net);

if iscell(net)==1
net1 = net{1};
R    = net1.Parameters(1,4);
else
net1 = net;
R    = net1.Parameters(1,4);
end



l = size(net1.X,1);
Ind = randperm(l);
Ind = Ind(1,1);

plot3(net1.X(Ind,:)./R,net1.Y(Ind,:)./R,net1.Z(Ind,:)./R,'.-');
xlim([-1.1 1.1])
ylim([-1.1 1.1])
zlim([-1.1 1.1])
axis square
set(gca,'XTick',[],'YTick',[],'ZTick',[])


hold off
