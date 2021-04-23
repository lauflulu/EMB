% measure DIB area and droplet volumes

clear all
close all
clc

%%
datafiles={'190925_full_100mM\data_full_100mM',...
        '190927_full_10mM\data_full_10mM',...
        '190928_full_1mM\data_full_1mM',...
        '190929_full_100uM\data_full_100uM',...
        '190930_ctrl1_100mM\data_ctrl1_100mM',...
        '190930_ctrl1_10mM\data_ctrl1_10mM',...
        '191002_ctrl1_1mM\data_ctrl1_1mM',...
        '191002_ctrl1_100uM\data_ctrl1_100uM',...
        '191010_ctrl2_100mM\data_ctrl2_100mM',...
        '191010_ctrl2_10mM\data_ctrl2_10mM'};
%% main loop
for p=3%:size(datafiles,2)
    
    load([datafiles{p},'.mat'])
    ID={data.name}; N=length(ID);
    T=length(data(1).t);
    for n=25%1:N
        
        Area=zeros(T,6)./0; %initialize nan array
        DIB=zeros(T,5)./0; %initialize nan array
        for t=49%1:T
            %% measure area
            try
                imageName=split(datafiles{p},'\');
                imageName=[imageName{1},'\o9-size_videos\',data(n).name,'-B2.tif'];
                I=imbinarize(imread(imageName, t),0.5);
                areaImage=I;
                A = regionprops(areaImage,'Area','centroid','PixelList');
                xy = cat(1,A.Centroid);
                area = cat(1,A.Area);
                X=data(n).x(t,:)'; X=X(~isnan(X));
                Y=data(n).y(t,:)'; Y=Y(~isnan(Y));
                R=data(n).r(t,:)'; R=R(~isnan(R));
                D=data(n).d(t,:)'; D=D(~isnan(D));

                for i=1:length(X)
                    centerDistance=sqrt((xy(:,1)-X(i)).^2+((xy(:,2)-Y(i)).^2))./R(i); 
                    Area(t,i)=area(((area<50000).*(centerDistance<1))==1);
                end

    %             figure(1)
    %             imshow(areaImage)
    %             hold on
    %             for i=1:length(A)
    %                 plot(xy(i,1),xy(i,2),'b*')
    %                 text(xy(i,1),xy(i,2),sprintf('%d',A(i).Area))
    %             end
    %             for i=1:length(X)
    %                 plot(X(i,1),Y(i,1),'r*')
    %             end
    %             hold off
    %%
                 figure(2)
                 imshow(I)

                bilayerImage=1-I;
                bilayerImage = bwmorph(bilayerImage,'skel',Inf);
                bilayerImage = bwmorph(bilayerImage,'branchpoints');
                s = regionprops(bilayerImage,'centroid');
                bps = cat(1,s.Centroid);
                for i=1:(length(X)-1)
                    minD1=ones(size(bps,1),1);
                    minD2=ones(size(bps,1),1);
                    for j=1:size(bps,1)
                        centerDistance=sqrt((xy(:,1)-X(i)).^2+((xy(:,2)-Y(i)).^2))./R(i); 
                        idx1=[1:length(A)]*(((area<50000).*(centerDistance<1))==1);
                        centerDistance=sqrt((xy(:,1)-X(i+1)).^2+((xy(:,2)-Y(i+1)).^2))./R(i+1);
                        idx2=[1:length(A)]*(((area<50000).*(centerDistance<1))==1);
                        D1=sqrt((A(idx1).PixelList(:,1)-bps(j,1)).^2+((A(idx1).PixelList(:,2)-bps(j,2)).^2));
                        minD1(j,1)=min(D1)/R(i);
                        D2=sqrt((A(idx2).PixelList(:,1)-bps(j,1)).^2+((A(idx2).PixelList(:,2)-bps(j,2)).^2));
                        minD2(j,1)=min(D2)/R(i+1);
                    end
                    idx=(((minD2<0.1).*(minD1<0.1))==1);
                    dibEdges=bps(idx,:);
                    
                    hold on
                    plot(dibEdges(:,1),dibEdges(:,2),'x','Color',i/5*[1 0 0],'LineWidth',3')
                    
                    dibLength=sqrt((dibEdges(1,1)-dibEdges(2,1)).^2+((dibEdges(1,2)-dibEdges(2,2)).^2));
                    DIB(t,i)=dibLength;
                end

            end
        data(n).area=Area;
        data(n).DIB=DIB;
        end
    end
    save([datafiles{p},'_v2.mat'],'data');
end

