clear all
close all
%%
load('20-04-21_Figure2.mat')
g = EMB_sim2g(y_R,[8,5,1]);
g=g(:,:,:,2:end);
time=t(2:end)/60;

meanG=EMB_g2meanG(g);
covG=EMB_g2covG(g);

meanY=squeeze(meanG(1,1,:,:));
meanR=squeeze(meanG(1,2,:,:));
stdY=squeeze(covG(1,1,:,:)).^0.5;
stdR=squeeze(covG(2,2,:,:)).^0.5;

[PIsga,perm] = EMB_g2piDIR(g(:,1:2,:,:),10); % without extrapolation since N=1000
%%
T=12;
close all
figure(1)
    for x=1:5
    subplot(4,5,x)
        hold all
        plot(time, squeeze(g(:,1,x,:)),'-k');
        plot(time, meanY(x,:), '-c','LineWidth',2);
        plot(time, meanY(x,:)+stdY(x,:), '--c','LineWidth',2);
        plot(time, meanY(x,:)-stdY(x,:), '--c','LineWidth',2);
        ylim([-0.1,1.2]); box('on')
    subplot(4,5,5+x)
        hold all
        plot(time, squeeze(g(:,2,x,:)),'-k')
        plot(time, meanR(x,:), '-r','LineWidth',2);
        plot(time, meanR(x,:)+stdR(x,:), '--r','LineWidth',2);
        plot(time, meanR(x,:)-stdR(x,:), '--r','LineWidth',2);
        ylim([-0.1,1.2]); box('on')
    subplot(4,5,10+x)
        hold all
        plot(time(1:T), squeeze(g(:,3,x,1:T)),'-k')
        plot([0,time(T)],[0.45e-3,0.45e-3],'--r','LineWidth',2)
        plot([0,time(T)],[4.5e-3,4.5e-3],'-r','LineWidth',2)
        plot([0,time(T)],[45e-3,45e-3],'--r','LineWidth',2)
        box('on')
    end
    subplot(4,1,4)
        hold all
        plot(time, PIsga', '-');set(gca,'ColorOrderIndex',1)
        %plot(time, PIsga+stdPIsga, '--');set(gca,'ColorOrderIndex',1)
        %plot(time, PIsga-stdPIsga, '--');
        box('on');% xlim([0,7.5]);ylim([-0.1,1.5]);
        xlabel('Time (h)');ylabel('I_{SGA}')
        legend('g1','g2','joint');%xticks([0:2.5:7.5]);