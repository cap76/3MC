
load('../InferPaths/Heterologous WLC/CorrectPriorTuneN/Run2/CorrectPrior_Tunep_3.mat')
obs{1,1} = 1:1:10; 
obs{1,2} = 20:1:100;

obs{2,1} = 1:1:25;
obs{2,2} = 30:1:100;

obs{3,1} = 1:1:42;
obs{3,2} = 50:1:100;

obs{4,1} = 1:1:70;
obs{4,2} = 80:1:100;

deltax = 0.31;
deltay = 0.31;
deltaz = 0.31;

for j = 1:4
subplot(4,2,j);
%plot3(Dat.Polymer.X(j,:),Dat.Polymer.Y(j,:),Dat.Polymer.Z(j,:),'k.-'),hold on,
PolymerPic([0,1,1,1],Dat.Polymer.X(j,[obs{j,1},obs{j,2}]),Dat.Polymer.Y(j,[obs{j,1},obs{j,2}]),Dat.Polymer.Z(j,[obs{j,1},obs{j,2}]),0.01),set(gca,'XTick',[],'YTick',[],'ZTick',[])
xlim([mean(Dat.Polymer.X(j,:))-deltax, mean(Dat.Polymer.X(j,:))+deltax])
ylim([mean(Dat.Polymer.Y(j,:))-deltay, mean(Dat.Polymer.Y(j,:))+deltay])
zlim([mean(Dat.Polymer.Z(j,:))-deltaz, mean(Dat.Polymer.Z(j,:))+deltaz])
%title('Sample trajectorie')
end



load HWLC_r/InCorrectPrior_Tunep_Tunehyp_2_missing.mat

aa=prctile(exp(Dat.MCMC.Y{1}(50000:2:end,:)-1),[0.5 99.5]);
EB1 = aa(:,1:98);
M = mean(exp(Dat.MCMC.Y{1}(50000:2:end,:)-1),1);
M1 = M(:,1:98);
T = exp(Dat.True.Ylat-1);
T1 = T(:,1:98);

subplot(4,2,[5 6]);
z = 1:1:98;
f = [[M1+(EB1(1,:)-M1)]'; flipdim([M1-(M1-EB1(2,:))]',1)];
fill([z'; flipdim(z',1)], f, [7 7 7]/8)
hold on; plot(z, M1, 'k-', 'LineWidth', 2); plot(z, T1, 'r-', 'LineWidth', 2)
ylabel('\kappa/L_c')
title('N = 2')
set(gca,'XTick',[]),ylim([0 1])


load HWLC_r/InCorrectPrior_Tunep_Tunehyp_3_missing.mat

aa=prctile(exp(Dat.MCMC.Y{1}(50000:2:end,:)-1),[0.5 99.5]);
EB1 = aa(:,1:98);
M = mean(exp(Dat.MCMC.Y{1}(50000:2:end,:)-1),1);
M1 = M(:,1:98);
T = exp(Dat.True.Ylat-1);
T1 = T(:,1:98);

subplot(4,2,[7 8]);
z = 1:1:98;
f = [[M1+(EB1(1,:)-M1)]'; flipdim([M1-(M1-EB1(2,:))]',1)];
fill([z'; flipdim(z',1)], f, [7 7 7]/8)
hold on; plot(z, M1, 'k-', 'LineWidth', 2); plot(z, T1, 'r-', 'LineWidth', 2)
ylabel('\kappa/L_c')
title('N = 4')
set(gca,'XTick',[]),ylim([0 1])