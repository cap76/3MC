meanfunc={'meanSum',{'meanConst','meanLinear'}}

X = []; Y = [];
for i = 1%:5
    load(['../Results/Partition/EmulatePartition_' num2str(i) '.mat'])
    X = [X;sort(Part.Ylat,2)];
    Y = [Y; Part.Zantrue];
end
%maxX = max(X')'/100;

Yall = Y;
Xall = X;
ind = find(isinf(Y)==0 & Y< 5e5);
X = X(ind,:);
Y = Y(ind,:);
%me = mean(Y);
%st = std(Y-me);
%Yhat = (Y-me)./st;
%Yall = (Yall-me)./st;
%me = min(Y);
%%Yhat = log(Y-me+1);
%Yhat = Y;

%[Yhat2 inds2] = sort(Yhat,'descend');
%[Yhat3 inds3] = sort(Yhat);

Yhat = Y;

%inds = unique([inds2(1:100)',inds3(1:100)',1:1:1000]);
inds = unique([1:1:1000]);

%Hyperparameters
hyp.cov   = log([1;1]);
hyp.lik   = log(.1);
hyp.mean  = [ones(98,1);mean(Yhat(inds,1))/1000];
hyp2      = feval(@minimize,hyp,@gp,-10000, @infExact, meanfunc, @covSEiso, @likGauss, X(inds,:),  Yhat(inds,1)/1000);
[Zan2 s2 m2pl s3] = gp(hyp2, @infExact, meanfunc, @covSEiso, @likGauss, X(inds,:),  Yhat(inds,1)/1000, X);


%Hyperparameters
hyp3.cov  = [1;1];
hyp3.lik  = log(.1);
hyp3.mean = mean(Yhat(inds,1)/1000); 
meanfunc2 = 'meanConst';
hyp3 = feval(@minimize,hyp3,@gp,-10000, @infExact, @meanConst, @covSEiso, @likGauss, X(inds,:),  Yhat(inds,1)/1000);
[Zan3 s2 m2pl s3] = gp(hyp3, @infExact, @meanConst, @covSEiso, @likGauss, X(inds,:),  Yhat(inds,1)/1000, X);


%[Zan s2 m2pl s3] = gp(hyp2, @infExact, meanfunc, @covSEiso, @likGauss, X(inds,:),  Yhat(inds,1), X);

EmulatParams.X = X(inds,:);
EmulatParams.Y = Yhat(inds,:);
EmulatParams.hyp1 = hyp2;
EmulatParams.hyp2 = hyp3;
EmulatParams.meanfunc1 = meanfunc;
EmulatParams.meanfunc2 = meanfunc2;
EmulatParams.scale = 1000;
%EmulatParams.me = me;
%EmulatParams.st = st;

save('EmulateParams.mat','EmulatParams')
