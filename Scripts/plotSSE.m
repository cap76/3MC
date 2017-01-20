addpath('/Users/christopher_penfold/Downloads/violin/')

load HWLC_r/InCorrectPrior_Tunep_Tunehyp_1_missing.mat
D(:,1) = ((repmat(Dat.True.Ylat,40001,1) - Dat.MCMC.Y{1}(100000:5:end,:)).^2,3920098,1);


load HWLC_r/InCorrectPrior_Tunep_Tunehyp_2_missing.mat
D(:,2) = ((repmat(Dat.True.Ylat,40001,1) - Dat.MCMC.Y{1}(100000:5:end,:)).^2,3920098,1);

load HWLC_r/InCorrectPrior_Tunep_Tunehyp_3_missing.mat
D(:,3) = ((repmat(Dat.True.Ylat,40001,1) - Dat.MCMC.Y{1}(100000:5:end,:)).^2,3920098,1);


load Results/HWLC/CorrectPrior_Tunep_1.mat
D(:,4) = reshape((repmat(Dat.True.Ylat,40001,1) - Dat.MCMC.Y{1}(100000:5:end,:)).^2,3920098,1);

load Results/HWLC/CorrectPrior_Tunep_2.mat
D(:,5) = reshape((repmat(Dat.True.Ylat,40001,1) - Dat.MCMC.Y{1}(100000:5:end,:)).^2,3920098,1);

load Results/HWLC/CorrectPrior_Tunep_3.mat
D(:,6) = reshape((repmat(Dat.True.Ylat,40001,1) - Dat.MCMC.Y{1}(100000:5:end,:)).^2,3920098,1);

load Results/HWLC/CorrectPrior_Tunep_4.mat
D(:,7) = reshape((repmat(Dat.True.Ylat,40001,1) - Dat.MCMC.Y{1}(100000:5:end,:)).^2,3920098,1);

load Results/HWLC/InCorrectPrior_Tunep_Tunehyp_1_simdata.mat
D(:,8) = reshape((repmat(Dat.True.Ylat,40001,1) - Dat.MCMC.Y{1}(100000:5:end,:)).^2,3920098,1);

load Results/HWLC/InCorrectPrior_Tunep_Tunehyp_2_simdata.mat
D(:,9) = reshape((repmat(Dat.True.Ylat,40001,1) - Dat.MCMC.Y{1}(100000:5:end,:)).^2,3920098,1);

load Results/HWLC/InCorrectPrior_Tunep_Tunehyp_3_simdata.mat
D(:,10) = ((repmat(Dat.True.Ylat,40001,1) - Dat.MCMC.Y{1}(100000:5:end,:)).^2,3920098,1);


load Results/HWLC/InCorrectPrior_Tunep_Tunehyp_4_simdata.mat
D(:,11) = ((repmat(Dat.True.Ylat,40001,1) - Dat.MCMC.Y{1}(100000:5:end,:)).^2,3920098,1);


load Results/HWLC/InCorrectPrior_Tunep_1.mat
D(:,12) = ((repmat(Dat.True.Ylat,40001,1) - Dat.MCMC.Y{1}(100000:5:end,:)).^2,3920098,1);

load Results/HWLC/InCorrectPrior_Tunep_2.mat
D(:,13) = ((repmat(Dat.True.Ylat,40001,1) - Dat.MCMC.Y{1}(100000:5:end,:)).^2,3920098,1);

load Results/HWLC/InCorrectPrior_Tunep_3.mat
D(:,14) = ((repmat(Dat.True.Ylat,40001,1) - Dat.MCMC.Y{1}(100000:5:end,:)).^2,3920098,1);

load Results/HWLC/InCorrectPrior_Tunep_4.mat
D(:,15) = ((repmat(Dat.True.Ylat,40001,1) - Dat.MCMC.Y{1}(100000:5:end,:)).^2,3920098,1);

violin(D(1::10:end,[4:end,1,2,3]))