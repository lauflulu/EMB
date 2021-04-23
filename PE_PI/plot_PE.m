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
    

%% PE
time=[0,6,12,24,90]; 
load(datafiles{3});
fusionTime=91;
g=EMB_data2g(data,fusionTime);
g=g(:,:,:,time+1);
time=time*5/60;
[N,I,X,T]=size(g);

if N>12
    m=flip(unique(ceil((6:11)/12*N)));
else
    m=flip(7:N);
end
tic
%[meanPE,stdPE] = EMB_extrapolatePE(@EMB_g2gaussPE,g,100,m,true);
[PE,meanPE,stdPE] = EMB_bootstrap(@(g)EMB_extrapolatePE(@EMB_g2gaussPE,g,100,m,false),g,500);
toc
%[PE,meanPE,stdPE] = EMB_bootstrap(@EMB_g2gaussPE,g,1000);
[Pcorr,stdPcorr]=EMB_PE2Pcorr(meanPE,stdPE);


meanG=EMB_g2meanG(g);
meanY=squeeze(meanG(1,1,:,T));
meanR=squeeze(meanG(1,2,:,T));

covG=EMB_g2covG(g);
stdY=squeeze(covG(1,1,:,T)).^0.5;
stdR=squeeze(covG(2,2,:,T)).^0.5;
%%
figure(5)
    subplot(2,3,1)
        hold all
        plot(1:5,squeeze(meanY),'o-b','LineWidth',2)
        plot(1:5,squeeze(meanY+stdY),'--b','LineWidth',2)
        plot(1:5,squeeze(meanY-stdY),'--b','LineWidth',2)
        plot(1:5,squeeze(meanR),'o-r','LineWidth',2)
        plot(1:5,squeeze(meanR+stdR),'--r','LineWidth',2)
        plot(1:5,squeeze(meanR-stdR),'--r','LineWidth',2)
        ylim([0,1.5]);xlim([0.5,5.5]);
        xlabel('Droplet')
        ylabel('FI (au)')
        box('on');

    subplot(2,3,2)
        hold all
        plot(1:5,squeeze(meanPE(1,:,T)),'o-b','LineWidth',2)
        plot(1:5,squeeze(meanPE(2,:,T)),'o-r','LineWidth',2)
        plot(1:5,squeeze(meanPE(3,:,T)),'o-k','LineWidth',2)
        plot(1:5,squeeze(meanPE(1,:,T)+stdPE(1,:,T)),'--b','LineWidth',2)
        plot(1:5,squeeze(meanPE(2,:,T)+stdPE(2,:,T)),'--r','LineWidth',2)
        plot(1:5,squeeze(meanPE(3,:,T)+stdPE(3,:,T)),'--k','LineWidth',2)
        plot(1:5,squeeze(meanPE(1,:,T)-stdPE(1,:,T)),'--b','LineWidth',2)
        plot(1:5,squeeze(meanPE(2,:,T)-stdPE(2,:,T)),'--r','LineWidth',2)
        plot(1:5,squeeze(meanPE(3,:,T)-stdPE(3,:,T)),'--k','LineWidth',2)
        plot(1:5,0.5*ones(1,5),'--k')
        plot(1:5,5*ones(1,5),'--k')
        %ylim([0,8]);
        xlim([0.5,5.5]);
        xlabel('Droplet')
        ylabel('\sigma_x (droplets)')
        box('on');
        
    subplot(2,3,3)
        hold all
        plot(1:5,Pcorr(1,:,T),'.-b','MarkerSize',20)
        plot(1:5,Pcorr(2,:,T),'.-r','MarkerSize',20)
        plot(1:5,Pcorr(3,:,T),'.-k','MarkerSize',20)
        plot(1:5,Pcorr(1,:,T)+stdPcorr(1,:,T),'--b')
        plot(1:5,Pcorr(2,:,T)+stdPcorr(2,:,T),'--r')
        plot(1:5,Pcorr(3,:,T)+stdPcorr(3,:,T),'--k')
        plot(1:5,Pcorr(1,:,T)-stdPcorr(1,:,T),'--b')
        plot(1:5,Pcorr(2,:,T)-stdPcorr(2,:,T),'--r')
        plot(1:5,Pcorr(3,:,T)-stdPcorr(3,:,T),'--k')
        xlim([0.5,5.5]);ylim([0,1]);box('on');
        xlabel('Droplet')
        ylabel('P_{corr}')
    
    for p=1:3
        subplot(2,3,3+p)
            hold all
            for t=1:T
                color=[0,0,0];
                color(p)=t/T;
                plot(1:5,squeeze(Pcorr(p,:,t))','.-',...
                    'Color',color,'LineWidth',1,'MarkerSize',20)
                plot(1:5,squeeze(Pcorr(p,:,t)+stdPcorr(p,:,t))','--',...
                    'Color',color,'LineWidth',1)
                plot(1:5,squeeze(Pcorr(p,:,t)-stdPcorr(p,:,t))','--',...
                    'Color',color,'LineWidth',1)
            end
            xlim([0.5,5.5]);
            box('on'); ylim([0,1]);
            xlabel('Time (h)')
            ylabel('P_{corr}')
    end