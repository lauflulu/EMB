clear all
close all
%%

load('20-04-21_Figure2.mat')
g = EMB_sim2g(y_R,[8,5]);
g=g(:,:,:,2:end);
time=t(2:end);

meanG=EMB_g2meanG(g);
covG=EMB_g2covG(g);

meanY=squeeze(meanG(1,1,:,:));
meanR=squeeze(meanG(1,2,:,:));
stdY=squeeze(covG(1,1,:,:)).^0.5;
stdR=squeeze(covG(1,2,:,:)).^0.5;

[PIsga,perm] = EMB_g2piSGA(g,100); % without extrapolation since N=1000
%%
figure(1)
    subplot(3,1,1)
        hold all
        plot(time, meanY, '-');set(gca,'ColorOrderIndex',1)
        plot(time, meanY+stdY, '--');set(gca,'ColorOrderIndex',1)
        plot(time, meanY-stdY, '--');set(gca,'ColorOrderIndex',1)
    subplot(3,1,2)
        hold all
        plot(time, meanR, '-');set(gca,'ColorOrderIndex',1)
        plot(time, meanR+stdR, '--');set(gca,'ColorOrderIndex',1)
        plot(time, meanR-stdR, '--');set(gca,'ColorOrderIndex',1)
    subplot(3,1,3)
        hold all
        plot(time, PIsga', '-');set(gca,'ColorOrderIndex',1)
        %plot(time, PIsga+stdPIsga, '--');set(gca,'ColorOrderIndex',1)
        %plot(time, PIsga-stdPIsga, '--');
        box('on');% xlim([0,7.5]);ylim([-0.1,1.5]);
        xlabel('Time (h)');ylabel('I_{SGA}')
        legend('g1','g2','joint');%xticks([0:2.5:7.5]);
        title(sprintf('N = %d',N));