
for chrom = 1:9

    X1 = [];
    Y1 = [];
    Z1 = [];
    X2 = [];
    Y2 = [];
    Z2 = [];    
    X3 = [];
    Y3 = [];
    Z3 = [];        
    
    for step = 16:25
load(['/Users/christopherpenfold/Desktop/Trunk/Bouquet/Bouquet_40nm_Chrom_' num2str(chrom) '_2_1_' num2str(step) '.mat'])

X1 = [X1;Data.X(:,[1,300])];
Y1 = [Y1;Data.Y(:,[1,300])];
Z1 = [Z1;Data.Z(:,[1,300])];

load(['/Users/christopherpenfold/Desktop/Trunk/Bouquet/Bouquet_40nm_Chrom_' num2str(chrom) '_2_3_' num2str(step) '.mat'])

X2 = [X2;Data.X(:,[1,300])];
Y2 = [Y2;Data.Y(:,[1,300])];
Z2 = [Z2;Data.Z(:,[1,300])];

load(['/Users/christopherpenfold/Desktop/Trunk/Free/Free_40nm_Chrom_' num2str(chrom) '_2_1_' num2str(step) '.mat'])

X3 = [X3;Data.X(:,[1,300])];
Y3 = [Y3;Data.Y(:,[1,300])];
Z3 = [Z3;Data.Z(:,[1,300])];

    end

    
    D(1,chrom) = mean( sqrt( ( X1(:,1)-X1(:,end) ).^2 + ( Y1(:,1)-Y1(:,end) ).^2 + ( Z1(:,1)-Z1(:,end) ).^2 ) );
    sD(1,chrom) = std( sqrt( ( X1(:,1)-X1(:,end) ).^2 + ( Y1(:,1)-Y1(:,end) ).^2 + ( Z1(:,1)-Z1(:,end) ).^2 ) );

    D(2,chrom) = mean( sqrt( ( X2(:,1)-X2(:,end) ).^2 + ( Y2(:,1)-Y2(:,end) ).^2 + ( Z2(:,1)-Z2(:,end) ).^2 ) );
    sD(2,chrom) = std( sqrt( ( X2(:,1)-X2(:,end) ).^2 + ( Y2(:,1)-Y2(:,end) ).^2 + ( Z2(:,1)-Z2(:,end) ).^2 ) );

    D(3,chrom) = mean( sqrt( ( X3(:,1)-X3(:,end) ).^2 + ( Y3(:,1)-Y3(:,end) ).^2 + ( Z3(:,1)-Z3(:,end) ).^2 ) );
    sD(3,chrom) = std( sqrt( ( X3(:,1)-X3(:,end) ).^2 + ( Y3(:,1)-Y3(:,end) ).^2 + ( Z3(:,1)-Z3(:,end) ).^2 ) );    
    
end

x = linspace(1,100,16);
scale = x(1,2) - x(1,1);
scale = scale/6;

errorbar(x,D(3,:),sD(3,:),'s')
hold on
errorbar(x+scale,D(1,:),sD(1,:),'o')
errorbar(x+2*scale,D(2,:),sD(2,:),'^')
set(gca,'XTick',x+0.2,'XTickLabel',{'I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI'})