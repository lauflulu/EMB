clear all
close all

%%
datafiles={'raw_data\data_full_100mM_v2.mat',...
        'raw_data\data_full_10mM_v2.mat',...
        'raw_data\data_full_1mM_v2.mat',...
        'raw_data\data_full_100uM_v2.mat',...
        'raw_data\data_ctrl1_100mM_v2.mat',...
        'raw_data\data_ctrl1_10mM_v2.mat',...
        'raw_data\data_ctrl1_1mM_v2.mat',...
        'raw_data\data_ctrl1_100uM_v2.mat',...
        'raw_data\data_ctrl2_100mM_v2.mat',...
        'raw_data\data_ctrl2_10mM_v2.mat'};
    
%% FI(t), with sender
fusionTime=91;
s=3;
load(datafiles{s});

a={data.r}; a=cat(3, a{:});
excluder=ones(size(a,3),1);
for n=1:length(data)
    if ~(sum(sum(isnan(a(1:fusionTime,:,n))))==0)
        excluder(n,1)=0;
    end
end
data=data(excluder==1);

% shuffle around data
a={data.YFP}; a=cat(3, a{:});
b={data.RFP}; b=cat(3, b{:});
c={data.Cy5}; c=cat(3, c{:});

% normalization by reference dye
FIy=a./c; FIy=FIy(1:fusionTime,1:6,:); FIy=permute(FIy,[3 2 1]);
FIy=FIy/max(mean(FIy,1),[],'all');
FIr=b./c; FIr=FIr(1:fusionTime,1:6,:); FIr=permute(FIr,[3 2 1]);
FIr=FIr/max(mean(FIr,1),[],'all');

% write in g
N=sum(excluder); I=2; T=fusionTime; X=6;
g=zeros(N,I,X,T);
g(:,1,:,:)=FIy; g(:,2,:,:)=FIr;

% mean g
meanG=EMB_g2meanG(g);
meanY=squeeze(meanG(1,1,:,:));
meanR=squeeze(meanG(1,2,:,:));

% cov g
covG=EMB_g2covG(g);
stdY=squeeze(covG(1,1,:,:)).^0.5;
stdR=squeeze(covG(2,2,:,:)).^0.5;
    
t=(0:90)/60*5;
figure(1)
    for i=1:6
    subplot(2,6,i)
        hold all
        plot(t,squeeze(g(:,1,i,:)),'-m','LineWidth',1)
        plot(t,meanY(i,:),'-','LineWidth',2,'Color',[0 0 1])
        plot(t,meanY(i,:)+stdY(i,:),'--','LineWidth',1,'Color',[0 0 1])
        plot(t,meanY(i,:)-stdY(i,:),'--','LineWidth',1,'Color',[0 0 1])
        ylim([0,inf]); ylabel('YFP (au)');box('on')
        xlabel('Time (h)');xlim([0,7.5]);
        ylim([-0.1,1.4]);
        xticks(0:2.5:7.5);

    subplot(2,6,i+6)
        hold all
        plot(t,squeeze(g(:,2,i,:)),'-c','LineWidth',1)
        plot(t,meanR(i,:),'-','LineWidth',2,'Color',[1 0 0])
        plot(t,meanR(i,:)+stdR(i,:),'--','LineWidth',1,'Color',[1 0 0])
        plot(t,meanR(i,:)-stdR(i,:),'--','LineWidth',1,'Color',[1 0 0])
        ylim([0,inf]); ylabel('RFP (au)')
        xlabel('Time (h)');xlim([0,7.5]); box('on')
        ylim([-0.1,1.8]);
        xticks(0:2.5:7.5);
    end
        
%% overview figure normalized traces of pos19_3
figure(2);
    
    for i=1:6
        w=(6-i)/5;
        subplot(1,2,1)
            hold all
            plot(t,squeeze(g(14,1,i,:)),'-c','LineWidth',2,'Color',(w)*[0 1 1])
            ylim([0,inf]); ylabel('RFP (au)')
            xlabel('Time (h)'); box('on')
            xticks(0:2.5:7.5);xlim([0,7.5]);

        subplot(1,2,2)
            hold all
            plot(t,squeeze(g(14,2,i,:)),'-r','LineWidth',2,'Color',(1-w)*[1 0 0])
            ylim([0,inf]); ylabel('RFP (au)')
            xlabel('Time (h)'); box('on')
            xticks(0:2.5:7.5);xlim([0,7.5]);

    end