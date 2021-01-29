close all
clear all

%%
datafiles={'data\data_full_100mM_v2.mat',...
        'data\data_full_10mM_v2.mat',...
        'data\data_full_1mM_v2.mat',...
        'data\data_full_100uM_v2.mat',...
        'data\data_ctrl1_100mM_v2.mat',...
        'data\data_ctrl1_10mM_v2.mat',...
        'data\data_ctrl1_1mM_v2.mat',...
        'data\data_ctrl1_100uM_v2.mat',...
        'data\data_ctrl2_100mM_v2.mat',...
        'data\data_ctrl2_10mM_v2.mat'};
    
S=length(datafiles);

%%
load(datafiles{3});
fusionTime=91;
g=EMB_data2g(data,fusionTime);
cov=EMB_g2covG(g);
[N,I,X,T]=size(g);

pearson=cov(1,2,:,:)./cov(1,1,:,:).^0.5./cov(2,2,:,:).^0.5;

%%
T=1;
figure(1)
    for x=1:X
        subplot(1,X,x)
            scatter(g(:,1,x,T),g(:,2,x,T))
            xlim([0,1.4]);ylim([0,1.8]);
            xlabel('g1');ylabel('g2');
            box('on'); title(sprintf('p_1_2 = %.3f', pearson(1,1,x,T)))
    end