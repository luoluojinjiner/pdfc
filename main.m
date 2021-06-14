clc
clear
close all
    I=imread(['F:\remote_sensing_image\remote_sensing_image\seg_train\train\src\','5.png']); %I was the original, add Image Path.
    I1=imread(['F:\remote_sensing_image\remote_sensing_image\seg_train\train\label\','5.png']);%I1 is GT figure, add Image Path.

m=20;
figure 
imshow(I);
I1=double(I1);

[data,lable] = supmean(I);
I1=imresize(uint8(I1),[400,400]);
ben=length(data(:,1));
data=data+1;
for i=1:ben  
    huh2(i) =I1(data(i,1),data(i,2));
end
huh2=huh2+1;
supmean1=data;
[m1 n1]=size(lable);label=lable;
data=double(data(:,3:5));
ben=length(data(:,1));

figure 
plot(data(:,1),data(:,2),'ob','MarkerFaceColor','r','markersize',6);
set(gcf,'PaperType','a3');
xlabel('R value','FontSize',18,'Fontname','Times newman');
ylabel('G value','FontSize',18,'Fontname','Times newman');
set(gca,'FontName','Times New Roman','FontSize',18);

hold on;

num_data = size(data,1);%Number of samples
num_datax=num_data;
num_d = size(data,2);%Sample dimension
num_d=2;
c=color(1);
c=c/255;
D1=supmean1(:,3:5);
IDX=length(D1(:,1));

he=distance1(data);

he1=he;
[ xyh1, index ] = sort(he,2);%Sort each row of data from smallest to largest, resulting in an index set of the nearest point for each point¡£
index2=index;
mz=30;
table=zeros(256,256);
table1=zeros(256,256);table2=zeros(256,256);table3=zeros(256,256);
D1=double(D1);
XS1=zeros(2,mz);
okol= (1:ben);
okol=okol';
supmeanx1=[supmean1 okol];
supmeanx1=sortrows(supmeanx1,5);

%The distance between each point and the surrounding DD point. This super parameter can be set artificially. This is set to 1, which means that all of the points are considered core points, and no edge points¡£
xyh=xyh1;
mm=mean(xyh(:,m));%The mean of the distance from the m-th point
zhong=zeros(2,num_data);
zhong1=zeros(1,num_data);
for i=1:num_data%Find the core point, Zhong is the core point, Zhong1 is the edge point, if % is greater than mm must be the edge point, if < mm£¬
      kj1=find(xyh(i,:) > mm, 1, 'first'); %Then the number of inner points in the circle with the radius of mm for each point should be calculated as the weight of the core points¡£
     quan1(i) =kj1;
if xyh(i,m)<=mm
zhong(1,i)=i;%Zhong is the core point.
kj=find(xyh(i,:) > mm, 1, 'first');
zhong(2,i)=kj; 
else
zhong1(i)=i;%Zhong1 is the edge point.
end
end
data1=zeros(num_data,num_d);
for i=1:length(zhong)
if zhong(1,i)~=0
data1(i,1)=data(zhong(1,i),1);
data1(i,2)=data(zhong(1,i),2);
quan(i)=zhong(2,i);
end
end

ind1 = find(data1(:,1)~=0) ;
 quan(find(quan==0))=[];
data1(all(data1==0,2),:)=[];
plot(data1(:,1),data1(:,2),'ok','MarkerFaceColor','g','markersize',6); 
axis([0 245 0 190]);
set(gcf,'PaperType','a3');
xlabel('R value','FontSize',18,'Fontname','Times newman');
ylabel('G value','FontSize',18,'Fontname','Times newman');
set(gca,'FontName','Times New Roman','FontSize',18);
hold off;
poi=1:length(ind1);
A=[poi',ind1];%The original location that makes up the core of the map.
zhong2=zhong1;
hs=zhong(1,:)';
hg=zhong(2,:)';hg1=hg;hg=hg/max(hg);
sdata=[data hg];
for i=1:255
    for j=1:255
        for r=1:ben
   table(data(r,1),data(r,2))=hg(r);%Scatter data of the first dimension
   table1(data(r,1),data(r,2))=hg1(r);    
   table2(data(r,1),data(r,2))=r;
   table3(data(r,1),data(r,2))=hs(r);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%% Draw heat map %%%%%%%%%%%%%%%%%%%%%%%%%%%
AQ=round(255/10);
h = fspecial('gaussian',AQ,6);%Gaussian filtering
N2=imfilter(table,h);
%N2=table;
ind = sub2ind(size(N2),data1(:,1),data1(:,2));
% ind = sub2ind(size(N2),data(:,1),data(:,2));
col = N2(ind);
figure
scatter(data1(:,1),data1(:,2),30,col,'filled');
box on;
%axis([5 245 1 235]);%island_042
%axis([15 210 50 220]);
%axis([0 232 0 230]);%lake_016
%axis([10 70 15 82]);%meadow_076
axis([0 245 0 190]);%desert_511
set(gcf,'PaperType','a3');
xlabel('R value','FontSize',18,'Fontname','Times newman');
ylabel('G value','FontSize',18,'Fontname','Times newman');
set(gca,'FontName','Times New Roman','FontSize',18);
N2=N2/max(max(N2));

%»æÖÆpcolorÍ¼ 
figure
h=pcolor(N2');shading interp
set(h,'edgecolor','none','facecolor','interp');
axis([5 250 25 215]);%desert_511
set(gcf,'PaperType','a3');
xlabel('R value','FontSize',18,'Fontname','Times newman');
ylabel('G value','FontSize',18,'Fontname','Times newman');
set(gca,'FontName','Times New Roman','FontSize',18);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Finish the heat map %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B = imregionalmax(N2);
B=double(B);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% peaks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[row, col] = find( B ~= 0 );
DSW=[row col];
DSW1=distance(DSW);
for i=1:length(row)
    for j=1:length(row)
       if DSW1(i,j)~=0
        DSW1(j,i)=0;
       end
    end
end
for i=1:length(row)
    for j=1:length(row)
    if DSW1(i,j)<25&&i~=j&&DSW1(i,j)~=0
    B(DSW(i,1),DSW(i,2))=0;  
    end
    end
end
zo=sum(sum(B));
%%%%%%%%%%%%%%% Make sure that each peak point remains the number of clusters %%%%%%%%%%%
X=[1 1 1 1 1 
    1 1 1 1 1
     1 1 1 1 1
     1 1 1 1 1
     1 1 1 1 1];
B=imfilter(B,X);
for i=1:255
    for j=1:255
   if table1(i,j)==0
       table1(i,j)=1;
   end
   end
end
table1=table1.*B;
tabley=table1;
[row, col] = find( table1 ~= 0 ); 
rows=cell(zo,1);cols=cell(zo,1);
   for i=1:zo
       sds= max(max(table1));
   [rows{i},cols{i}]=find(sds==(table1));rows{i}=rows{i}(1);cols{i}=cols{i}(1);
 %%%%%%%%%%%%%% 26 neighborhood %%%%%%%%%%%%% 
 %The number of neighborhood points is determined
 linyu=9;
   for x=1:linyu
       for y=1:linyu
         oping=rows{i}-((linyu+1)/2)+x;oping1=cols{i}-((linyu+1)/2)+y;
         if  oping<=0
             oping=1;
         end
         if  oping>256
             oping=256;
         end
           if  oping1<=0
             oping1=1;
         end
         if  oping1>256
             oping1=256;
         end
       table1(oping,oping1)=0;
       end
   end
   end
for i=1:zo
    table1(rows{i},cols{i})=1;
end
   table1=table1.*tabley;
[row, col] = find( table1 ~= 0 ); 
%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%
ss=sum(sum(table1));
center1=cell(2,length(row));
for sd=1:length(row)
    center1{2,sd}=table1(row(sd),col(sd)) ;%sf is the peak point density weight
    sf(sd)=table1(row(sd),col(sd)) ;
    se(sd)=table1(row(sd),col(sd))/ss ;%se is the normalized value of the peak point
    sx(sd)=table2(row(sd),col(sd));
    center1{1,sd}=table2(row(sd),col(sd));
end
DATA=center1(1,:);
for i=1:m1
    for j=1:n1
        for x=1:length(DATA)
            if label(i,j)==center1{1,x}
        IM(i,j,1)=255;
        IM(i,j,2)=0;
        IM(i,j,3)=0;
            end
        end
    end
end
supim1=uint8(cat(3,IM(:,:,1),IM(:,:,2),IM(:,:,3)));%IM is the three-channel image of the superpixel segmentation graph



hold on

hold off
figure,plot(data(:,1),data(:,3),'ok','MarkerFaceColor','b','markersize',4);%draw it
set(gca,'FontName','Times New Roman','FontSize',18);
hold on
for j=1:length(DATA)
for i=1:ben
    if i==center1{1,j}
    c = num2str(i);%Êý×Ö×ª×Ö·û
    plot(data(i,1),data(i,3),'ob','MarkerFaceColor','r','markersize',6);%draw it
    end
end
end
hold off
%%%%%%%%%%%%%%%%%%%%%% End of the drawing %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% Calculate the gravity from each core point to the peak %%%%%%%%%%%%%%%%

for sd=1:length(DATA)
for i=1:length(ind1)
       rou=zhong(2,ind1(i));
       dij=(he(center1{1,sd},zhong(1,ind1(i))));
       F(sd,i)=(rou+center1{2,sd})/dij.^2;
end
end
for i=1:length(poi)
       [sl,p(i)]=max(F(:,i));
end
xvc=cell(1,length(DATA));xva=cell(1,length(DATA));
for sd=1:length(DATA)
    for c=1:length(p)
        if p(c)==sd
    xvc{sd}=[xvc{sd};zhong(1,ind1(c))];%sx is the data number of the peak 
    xva{sd}=[xva{sd};zhong(2,ind1(c))];
        end
    end
end

for i=1:length(DATA)
   means(i)=mean(xva{i});
end
k_value=length(DATA);
asds=cell2mat(center1);
cen=data(asds(1,:),:);%cen is the peak point, the initial center point¡£
datah=data(ind1,:);%datah is a set of core points
huh=PDFC(datah,cen,quan,means,k_value);

centerx=cell(1,length(DATA));
for i=1:length(huh) 
   for j=1: length(DATA)
    if j==huh(i)
        centerx{j}=[centerx{j};ind1(i)];
    end
    end
end

ind2 = find(zhong2~=0) ;
L=0;
while 1
ol=length(ind2);
        for j=1:ol
        if isempty(ind2)
        break;
        end
for y=2:length(huh)
            for i=1:k_value
            C=centerx{:,i};
            B=index2(ind2(j),y);%Put several points into each cluster at a time, using the nearest point as the criterion
            tf=ismember(B,C);%Places the boundary points to be allocated into the corresponding class
            if tf==1
              L=L+1;
              centerx{i}=[centerx{i};ind2(j)];
              ind2(j)=0; 
              break;
            end
            end
            if tf==1 
               break;
            end
end
            if all(ind2==0)
               break;
            end  
        end  
   if all(ind2==0)
      break;
   end      
end

%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Data points are returned to the graph to display the final result graph %%%%%%%%%%%%
sd=color(1);
centerq=cell(1,length(DATA));
for i=1:k_value
    for j=1:num_datax
        if ismember(j,centerx{i})
           centerq{i}=[centerq{i};quan1(j)];
           huh(j)=i;
    supmean1(j,3)=sd(i,1);supmean1(j,4)=sd(i,2);supmean1(j,5)=sd(i,3);
        end
    end
end

IM=zeros(m1,n1,3);
for i=1:m1
    for j=1:n1
        IM(i,j,1)=uint8(supmean1(label(i,j),3));
        IM(i,j,2)=uint8(supmean1(label(i,j),4));
        IM(i,j,3)=uint8(supmean1(label(i,j),5));
    end
end

supim=uint8(cat(3,IM(:,:,1),IM(:,:,2),IM(:,:,3)));%IM is the three-channel image of the superpixel segmentation graph
figure;imshow(supim);
huh1=zeros(ben,1);
for i=1:ben
    for j=1:length(centerx) 
    if ismember(i,centerx{j})
    huh1(i,1)=j;
    end
    end
end



IM=zeros(m1,n1,3);
for i=1:m1
    for j=1:n1
        IM(i,j,1)=uint8(supmean1(label(i,j),3));
        IM(i,j,2)=uint8(supmean1(label(i,j),4));
        IM(i,j,3)=uint8(supmean1(label(i,j),5));
    end
end

%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% Draw a scatter diagram of clustering %%%%%%%%%%%%%
c=color(1);
c=c/255;
Clusters=zeros(length(xvc),num_datax);

for i=1:k_value
for j=1:length(huh1)
    if i==huh1(j)
Clusters(i,j)=j;
    else
Clusters(i,j)=0;
    end
end
end
figure
for i=1:k_value
    ClusK=setdiff(Clusters(i,:),0);
    plot(data(ClusK,1),data(ClusK,2),'.','color',c(i,:),'markersize',18);%
    box on;
    hold on;
end
axis([0 245 25 215]);
set(gcf,'PaperType','a3');
xlabel('R value','FontSize',18,'Fontname','Times newman');
ylabel('G value','FontSize',18,'Fontname','Times newman');
set(gca,'FontName','Times New Roman','FontSize',18);
%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% indicators %%%%%%%%%%%%%%%%
std1=0;std2=0;
% sdd=centerq{i};
 %ms=centerq{2}/max(centerq{2});
for i=1:k_value
p = imhist(centerq{i}/max(centerq{i}));
p(p==0) = [ ];
p = p ./ numel(centerq{i}/max(centerq{i}));
H1 = -sum(p.*log2(p));

std2=std2+H1;
std1=std1+std(centerq{i}/max(centerq{i}));
end
std1=std1/k_value;%
P4 = figure;clf;
[silh3,h3] = silhouette(D1,huh,'sqeuclidean');%The contour coefficient is plotted
scindex=mean(silh3);

xs1(m)=std2/scindex;
IIindex=std1*scindex;
xs(m)=IIindex;
m=0;quan=0;means=0;std1=0;F=0;p=0;p1=0;scindex=0;%clear¡£

MIhat=nmi(huh',huh2);%Normalized Mutual information
[AR,RI,MI,HI]=RandIndex(huh',huh2);%Rand index
fprintf('%f  ',AR,RI,MI,HI);
[TP,FN,FP,TN]=New_index(huh',huh2);
fprintf('\n')
fprintf('%f  ',TP,FN,FP,TN); 
fprintf('\n')
ACC= (TP+TN)/(TP+TN+FN+FN);
P=(TP+FN);
N=(FP+TN);
TPR=TP/P;FPR=FP/N;
P=TP/(TP+FP);R=TP/(TP+FN);
FM=sqrt((TP/(TP+FP))*(TP/(TP+FN))); 
huh2=[];

