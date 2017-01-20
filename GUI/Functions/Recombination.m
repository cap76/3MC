function [PR] = Recombination(net);


    warning('Only one set of trajectories given. Treating as homologoues and partitioning data')
                
        X1 = []; Y1 = []; Z1 = [];
        for i = 1:length(net)        
            X=[X1;net{i}.X];
            Y=[Y1;net{i}.Y];
            Z=[Z1;net{i}.Z];
        end
         
        %Randomly permute the data
        Ind1 = randperm(size(X,1));
        Ind2 = randperm(size(X,1));        
        %Take one half of the data
 
        X1 = X(Ind1,:);
        Y1 = Y(Ind1,:);
        Z1 = Z(Ind1,:);    
        X2 = X(Ind2,:);
        Y2 = Y(Ind2,:);
        Z2 = Z(Ind2,:);              
        
        Parameters1 = net{1}.Parameters;
        R1 = Parameters1(1,4);
        R2 = R1;
        
        BS = 10;
        proximity = 50;                

%X1 = []; Y1 = []; Z1 = [];
%for i = 1:length(net1)        
%    X1=[X1;net1{i}.X];
%    Y1=[Y1;net1{i}.Y];
%    Z1=[Z1;net1{i}.Z];
%end

        Parameters = net{1}.Parameters;
        R1 = Parameters(1,4);


%X2 = []; Y2 = []; Z2 = [];
%for i = 1:length(net2) 
%    X2=[X2;net2{i}.X];
%    Y2=[Y2;net2{i}.Y];
%    Z2=[Z2;net2{i}.Z];
%end


if size(X1,1)<50000 || BS<10
   warning('Number of samples and bootstrapped samples is low. Increase for greater accuracy') 
end

load('./Output/ssDNA.mat','Data')
 
ssDNA   = Data;
ssDNA.X = ssDNA.X(:,:)/(2*R1);
ssDNA.Y = ssDNA.Y(:,:)/(2*R1);
ssDNA.Z = ssDNA.Z(:,:)/(2*R1);
 
[l1 l2] = size(X1);
 
%Empty cells for storage
 PR   = zeros(BS,l2);
 maxD = (proximity/(2*R1))^2;

     for loopno = 1:BS %Bootstrap samples by permuting trajectories

        I1  = randperm(l1);
        I2  = randperm(l1);

        XX1 = X1(I1,:)./(2*R1);
        XX2 = X2(I2,:)./(2*R1);
        YY1 = Y1(I1,:)./(2*R1);
        YY2 = Y2(I2,:)./(2*R1);
        ZZ1 = Z1(I1,:)./(2*R1);
        ZZ2 = Z2(I2,:)./(2*R1);

        %SS = zeros(size(ssDNA.X,1),100);
        
        pTel = zeros(l1,l2);
        
        for MHstep = 1:size(X1,1)
            
            %if double(int64(MHstep/1000))==(MHstep/1000)
            %disp(['Step ' num2str(MHstep)])            
            %end          
    
            Dis1  = [(XX2(MHstep,:)-XX1(MHstep,:));(YY2(MHstep,:)-YY1(MHstep,:));(ZZ2(MHstep,:)-ZZ1(MHstep,:))];        
            Dis2 = repmat([ssDNA.X(MHstep,end);ssDNA.Y(MHstep,end);ssDNA.Z(MHstep,end)],1,l2);
    
            pro = sum(sqrt((Dis1-Dis2).^2),1);
            pTel(MHstep,:) = (pro<proximity/3200);
    
        end
    
        
        for i = 1:l2      
            [phat,pci] = binofit(length(find(pTel(:,i)>0)),length(pTel(:,i)));
            PR(loopno,i) = phat;
        end
                
     end