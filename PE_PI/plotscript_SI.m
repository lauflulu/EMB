clear all
close all

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

%% calculations, cf. function descriptions for details
b=0.05; % bin size for binPDF in percent/100, scaled to max(meanG)
t=49; % time point for profile plots etc.
B=20; % number of bins for gaussPDF
Boot=100; % number of bootstrap samples should be >1000 for good quality

load(datafiles{3});

time=data(1).t/60; time=time(1,1:91);

g=EMB_data2g(data,fusionTime);
[binPDF,maxGbin]=EMB_g2binPDF(g,b);

bPDFy=squeeze(sum(binPDF(:,:,:,t),2));
bPDFr=squeeze(sum(binPDF(:,:,:,t),1));

[gaussPDF,maxGgauss]=EMB_g2gaussPDF(g,B);
nPDFy=squeeze(sum(gaussPDF(:,:,:,t),2)); nPDFy=nPDFy./max(max(nPDFy))*0.5;
nPDFr=squeeze(sum(gaussPDF(:,:,:,t),1)); nPDFr=nPDFr./max(max(nPDFr))*0.5;

meanG=EMB_g2meanG(g);
meanY=squeeze(meanG(1,1,:,:));
meanR=squeeze(meanG(1,2,:,:));

covG=EMB_g2covG(g);
stdY=squeeze(covG(1,1,:,:)).^0.5;
stdR=squeeze(covG(2,2,:,:)).^0.5;

[PE,perm] = EMB_g2gaussPE(g);
P=size(perm,1);
[PI,perm] = EMB_g2cheatPI(g,B);

[~,meanPE,stdPE] = EMB_bootstrap(@EMB_g2gaussPE,g,Boot);

[~,meanPI,stdPI] = EMB_bootstrap(@(g)EMB_g2cheatPI(g,B),g,Boot); % this may take a while!

%% binPDFs 
figure(4)
    hy=subplot(2,4,2);
        cmap=zeros(256,3);
        for i=1:256
            w=(256-i)/255;
            cmap(i,:)=w*[1 1 1] + (1-w)*[0 0.51 0.79];
        end
        hold all
        imagesc(bPDFy,'XData', [1 5], 'YData', [0 maxGbin(1,1)])
        set(gca, 'Ydir', 'normal')
        colormap(hy,cmap)
        cby=colorbar;

        for i=1:5
            plot(nPDFy(:,i)+i,linspace(0,maxGgauss(1,1),B),'-b')
        end
        box('on');
        xlabel('Droplet Index');ylabel(cby,'p(g|x)')
        ylim([0,maxGgauss(1,1)]);

    subplot(2,4,1)
        hold all
        pdf_g1=sum(bPDFy,2); pdf_g1=pdf_g1/sum(pdf_g1)*length(pdf_g1)/maxGbin(1,1);
        yBins=linspace(0,maxGbin(1,1),27);
        
        pdf_G1=sum(nPDFy,2); pdf_G1=pdf_G1/sum(pdf_G1)/maxGgauss(1,1)*B;
 
        barh(yBins,pdf_g1)
        plot(pdf_G1,linspace(0,maxGgauss(1,1),B))
        ylim([0,maxGgauss(1,1)])
        ylabel('FI (au)'); box('on');
    
    subplot(2,4,6)
        hold all
        px=sum(bPDFy,1); px=px/sum(px);
        p_x=sum(nPDFy,1); p_x=p_x/sum(p_x);
        
        bar(1:5,px)
        plot(1:5,p_x)
        box('on');

    
    hr=subplot(2,4,4);
        cmap=zeros(256,3);
        for i=1:256
            w=(256-i)/255;
            cmap(i,:)=w*[1 1 1] + (1-w)*[0.86 0.03 0];
        end
        hold all
        imagesc(bPDFr,'XData', [1 5], 'YData', [0 maxGbin(1,2)])
        set(gca, 'Ydir', 'normal')
        colormap(hr,cmap)
        cbr=colorbar;
        for i=1:5
            plot(nPDFr(:,i)+i,linspace(0,maxGgauss(2,1),B),'-r')
        end
        
        box('on'); ylim([0 maxGgauss(2,1)]);
        ylabel('FI (au)');xlabel('Droplet Index');
        ylabel(cbr,'p(g|x)')
        
    subplot(2,4,3)
        hold all
        pdf_g2=sum(bPDFr,2); pdf_g2=pdf_g2/sum(pdf_g2)*length(pdf_g2)/maxGbin(1,2);
        rBins=linspace(0,maxGbin(1,2),34);
        pdf_G2=sum(nPDFr,2); pdf_G2=pdf_G2/sum(pdf_G2)/maxGgauss(2,1)*B;
        barh(rBins,pdf_g2)
        plot(pdf_G2,linspace(0,maxGgauss(2,1),B))
        ylim([0,maxGgauss(2,1)])
        ylabel('FI (au)');box('on');
    
    subplot(2,4,8)
        hold all
        px=sum(bPDFr,1); px=px/sum(px);
        p_x=sum(nPDFr,1); p_x=p_x/sum(p_x);
        bar(1:5,px)
        plot(1:5,p_x)
        box('on');

%% PE
figure(5)
    subplot(3,1,1)
        hold all
        yyaxis left
        plot(1:5,squeeze(meanY(:,t)),'o-b','LineWidth',2)
        plot(1:5,squeeze(meanY(:,t)+stdY(:,t)),'--b','LineWidth',2)
        plot(1:5,squeeze(meanY(:,t)-stdY(:,t)),'--b','LineWidth',2)
        ylim([0,0.08]);
        yyaxis right
        plot(1:5,squeeze(meanR(:,t)),'o-r','LineWidth',2)
        plot(1:5,squeeze(meanR(:,t)+stdR(:,t)),'--r','LineWidth',2)
        plot(1:5,squeeze(meanR(:,t)-stdR(:,t)),'--r','LineWidth',2)
        ylim([0,0.16]);xlim([0.5,5.5]);
        xlabel('Droplet')
        ylabel('FI (au)')
        box('on');

    subplot(3,1,[2,3])
        hold all
        plot(1:5,squeeze(meanPE(1,1,:,t)),'o-b','LineWidth',2)
        plot(1:5,squeeze(meanPE(1,2,:,t)),'o-r','LineWidth',2)
        plot(1:5,squeeze(meanPE(1,3,:,t)),'o-k','LineWidth',2)
        plot(1:5,squeeze(meanPE(1,1,:,t)+stdPE(1,1,:,t)),'--b','LineWidth',2)
        plot(1:5,squeeze(meanPE(1,2,:,t)+stdPE(1,2,:,t)),'--r','LineWidth',2)
        plot(1:5,squeeze(meanPE(1,3,:,t)+stdPE(1,3,:,t)),'--k','LineWidth',2)
        plot(1:5,squeeze(meanPE(1,1,:,t)-stdPE(1,1,:,t)),'--b','LineWidth',2)
        plot(1:5,squeeze(meanPE(1,2,:,t)-stdPE(1,2,:,t)),'--r','LineWidth',2)
        plot(1:5,squeeze(meanPE(1,3,:,t)-stdPE(1,3,:,t)),'--k','LineWidth',2)
        plot(1:5,0.5*ones(1,5),'--k')
        plot(1:5,5*ones(1,5),'--k')
        ylim([0,8]);xlim([0.5,5.5]);
        xlabel('Droplet')
        ylabel('\sigma_x (droplets)')
        box('on');
    
%% PI(t)

figure(6)
    hold all
    plot(time,squeeze(meanPI(1,1,:)),'-b','MarkerSize',5,'LineWidth',2)
    plot(time,squeeze(meanPI(1,1,:)+stdPI(1,1,:)),'--b','MarkerSize',5,'LineWidth',2)
    plot(time,squeeze(meanPI(1,1,:)-stdPI(1,1,:)),'--b','MarkerSize',5,'LineWidth',2)
    plot(time,squeeze(meanPI(1,2,:)),'-r','MarkerSize',5,'LineWidth',2)
    plot(time,squeeze(meanPI(1,2,:)+stdPI(1,1,:)),'--r','MarkerSize',5,'LineWidth',2)
    plot(time,squeeze(meanPI(1,2,:)-stdPI(1,1,:)),'--r','MarkerSize',5,'LineWidth',2)
    plot(time,squeeze(meanPI(1,3,:)),'-k','MarkerSize',5,'LineWidth',2)
    plot(time,squeeze(meanPI(1,3,:)+stdPI(1,1,:)),'--k','MarkerSize',5,'LineWidth',2)
    plot(time,squeeze(meanPI(1,3,:)-stdPI(1,1,:)),'--k','MarkerSize',5,'LineWidth',2)
    xlabel('Time (h)')
    ylabel('PI (bits)')
    ylim([0,1.5]);
    xlim([0,7.5]);xticks(0:2.5:7.5);
    box('on')
        
%% cheat PI vs. B

B=[2,5,10,20,50,100,1000];
PI91=zeros(3,length(B));
for i=1:length(B)
    [PI,perm] = EMB_g2cheatPI(g,B(i));
    PI91(:,i)=PI(:,91);
end
figure(11)
    hold all
    plot(B,PI91,'x-')      

%% droplet and bilayer size
[Rb,Rd,Ab,Vd,rl,rr,d] = EMB_data2geometry(data);

meanAb=mean(Ab,2:3,'omitnan');
stdAb=std(Ab,[],2:3,'omitnan');

meanVd=mean(Vd,2:3,'omitnan');
stdVd=std(Vd,[],2:3,'omitnan');
CVVd=stdVd./meanVd;

meanVd_intra=mean(Vd,2,'omitnan');
stdVd_intra=std(Vd,[],2,'omitnan');
CVVd_intra=mean(stdVd_intra./meanVd_intra,3,'omitnan');

meanVd_inter=mean(Vd,2,'omitnan');
stdVd_inter=std(meanVd_inter,[],3,'omitnan');
CVVd_inter=stdVd_inter./mean(meanVd_inter,3,'omitnan');

theta=(asin(rl)+asin(rr))/2/pi*180;
theta2=180-acos(-(d.^2-(Rb.*rl.^(-1)).^2-(Rb.*rr.^(-1)).^2)./(Rb.*rr.^(-1))./(Rb.*rl.^(-1)))/pi*180;
D=Rb.*((rl.^(-2)-1).^0.5 + (rr.^(-2)-1).^0.5);
%%
figure(7)     
    subplot(3,2,1)
        hold all
        plot(time,meanAb,'-r','LineWidth',3)
        plot(time,meanAb+stdAb,'-b','LineWidth',3)
        plot(time,meanAb-stdAb,'-b','LineWidth',3)
        box('on');xlabel('Time (h)');ylabel('DIB Area (?m^2)');
    
    subplot(3,2,2)
        hold all
        plot(time,meanVd,'-r','LineWidth',3)
        plot(time,meanVd+stdVd,'-b','LineWidth',3)
        plot(time,meanVd-stdVd,'-b','LineWidth',3)
        box('on');xlabel('Time (h)');ylabel('Droplet Volume (?m^3)');
        
    subplot(3,2,3)
        hold all
        plot(squeeze(rl(1,:,:)),squeeze(rr(1,:,:)),'ob')
        plot(squeeze(rl(13,:,:)),squeeze(rr(49,:,:)),'or')
        plot(squeeze(rl(49,:,:)),squeeze(rr(91,:,:)),'ok')
        box('on');xlabel('Rb/Rl');ylabel('Rb/Rr');
    subplot(3,2,4)
        hold all
        plot(squeeze(Rb(1,:,:)),squeeze(theta(1,:,:)),'ob');
        plot(squeeze(Rb(13,:,:)),squeeze(theta(13,:,:)),'or');
        plot(squeeze(Rb(49,:,:)),squeeze(theta(49,:,:)),'ok');
        box('on');xlabel('Rb');ylabel('theta');
    subplot(3,2,5)
        hold all
        plot(squeeze(d(1,:,:)),squeeze(D(1,:,:)),'ob')
        plot(squeeze(d(13,:,:)),squeeze(D(49,:,:)),'or')
        plot(squeeze(d(49,:,:)),squeeze(D(91,:,:)),'ok')
        plot(0:140,0:140,'-g','LineWidth',3)
        box('on');xlabel('d');ylabel('D');
    subplot(3,2,6)
        hold all
        plot(squeeze(theta(1,:,:)),squeeze(theta2(1,:,:)),'ob');
        plot(squeeze(theta(13,:,:)),squeeze(theta2(13,:,:)),'or');
        plot(squeeze(theta(49,:,:)),squeeze(theta2(49,:,:)),'ok');
        plot(0:180,0:180,'-g','LineWidth',3)
        box('on');xlabel('Rb');ylabel('theta2');

%% fit cumulative pdf
fusionTime=91;

load(datafiles{3});
g=EMB_data2g(data,fusionTime);

t=49;

fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0],...
               'StartPoint',[0.005,0.05]);
ft = fittype('0.5*(erf((x-b)./(sqrt(2).*a))+1)','options',fo);

% YFP
Xy=squeeze(sort(g(:,1,:,t),1)); % cumulative pdf
Yy=[1:size(g,1)]./size(g,1); Yy=Yy';

legy=cell(5,1);
xy=linspace(0,0.1,1000)';
yy=zeros(1000,5);
for i=1:5
    [fy,gofy]=fit(Xy(:,i),Yy,ft);
    legy{i,1}=sprintf('R^2 = %1.3f\nCVfit/CVs = %2.1f/%2.1f=%1.2f',...
        gofy.rsquare,fy.a/fy.b*100,std(Xy(:,i))/mean(Xy(:,i))*100,...
        fy.a/fy.b*100/(std(Xy(:,i))/mean(Xy(:,i))*100));
    yy(:,i)=0.5*(erf((xy - fy.b)/(2^0.5*fy.a))+1);
end

% RFP
Xr=squeeze(sort(g(:,2,:,t),1)); % cumulative pdf
Yr=[1:size(g,1)]./size(g,1); Yr=Yr';

legr=cell(5,1);
xr=linspace(0,0.2,1000)';
yr=zeros(1000,5);
for i=1:5
    [fr,gofr]=fit(Xr(:,i),Yr,ft);
    legr{i,1}=sprintf('R^2 = %1.3f\nCVfit/CVs = %2.1f/%2.1f=%1.2f',...
        gofr.rsquare,fr.a/fr.b*100,std(Xr(:,i))/mean(Xr(:,i))*100,...
        fr.a/fr.b*100/(std(Xr(:,i))/mean(Xr(:,i))*100));
    yr(:,i)=0.5*(erf((xr - fr.b)/(2^0.5*fr.a))+1);
end

% figure
figure(3)
    subplot(2,1,1)
        hold all
        for i=1:5
            w=(5-i)/4;
            col=w*[0 0 1]+(1-w)*[0 1 1];

            plot(Xy(:,i),Yy,'.','MarkerEdgeColor',col) 
            hy(i)=plot(xy,yy(:,i),'-','Color',col);
        end
        box('on');ylabel('cumulative pdf');xlabel('FI (au)'); ylim([0 1.1])
        legend(hy,legy)

    subplot(2,1,2)
        hold all
        for i=1:5
            w=(5-i)/4;
            col=w*[1 0 0]+(1-w)*[1 1 0];
            plot(Xr(:,i),Yr,'.','MarkerEdgeColor',col) 
            hr(i)=plot(xr,yr(:,i),'-','Color',col);
        end
        box('on');ylabel('cumulative pdf');xlabel('FI (au)');ylim([0 1.1])
        legend(hr,legr)
    
%% FI(t), with sender
% load(datafiles{3});
% a={data.r}; a=cat(3, a{:});
% excluder=ones(size(a,3),1);
% for n=1:length(data)
%     if ~(sum(sum(isnan(a(1:fusionTime,:,n))))==0)
%         excluder(n,1)=0;
%     end
% end
% data=data(excluder==1);
% 
% % shuffle around data
% a={data.YFP}; a=cat(3, a{:});
% b={data.RFP}; b=cat(3, b{:});
% c={data.Cy5}; c=cat(3, c{:});
% 
% % normalization by reference dye
% FIy=a./c; FIy=FIy(1:fusionTime,1:6,:); FIy=permute(FIy,[3 2 1]);
% FIr=b./c; FIr=FIr(1:fusionTime,1:6,:); FIr=permute(FIr,[3 2 1]);
% 
% % write in g
% N=sum(excluder); I=2; T=fusionTime; X=6;
% g=zeros(N,I,X,T);
% g(:,1,:,:)=FIy; g(:,2,:,:)=FIr;
% 
% % mean g
% meanG=EMB_g2meanG(g);
% meanY=squeeze(meanG(1,1,:,:));
% meanR=squeeze(meanG(1,2,:,:));
% 
% % cov g
% covG=EMB_g2covG(g);
% stdY=squeeze(covG(1,1,:,:)).^0.5;
% stdR=squeeze(covG(2,2,:,:)).^0.5;
%     
% t=(0:90)/60*5;
% figure(2)
%     for i=1:6
%     subplot(2,6,i)
%         hold all
%         plot(t,squeeze(g(:,1,i,:)),'-m','LineWidth',1)
%         plot(t,meanY(i,:),'-','LineWidth',2,'Color',[0 0 1])
%         plot(t,meanY(i,:)+stdY(i,:),'--','LineWidth',1,'Color',[0 0 1])
%         plot(t,meanY(i,:)-stdY(i,:),'--','LineWidth',1,'Color',[0 0 1])
%         ylim([0,inf]); ylabel('YFP (au)');box('on')
%         xlabel('Time (h)');xlim([0,7.5]);
%         ylim([0,max(max(max(g(:,1,:,:),[],1),[],3),[],4)]);
%         xticks(0:2.5:7.5);
% 
%     subplot(2,6,i+6)
%         hold all
%         plot(t,squeeze(g(:,2,i,:)),'-c','LineWidth',1)
%         plot(t,meanR(i,:),'-','LineWidth',2,'Color',[1 0 0])
%         plot(t,meanR(i,:)+stdR(i,:),'--','LineWidth',1,'Color',[1 0 0])
%         plot(t,meanR(i,:)-stdR(i,:),'--','LineWidth',1,'Color',[1 0 0])
%         ylim([0,inf]); ylabel('RFP (au)')
%         xlabel('Time (h)');xlim([0,7.5]); box('on')
%         ylim([0,max(max(max(g(:,2,:,:),[],1),[],3),[],4)]);
%         xticks(0:2.5:7.5);
%     end
%         
% %% overview figure normalized traces of pos19_3
% figure(1);
%     
%     for i=1:6
%         w=(6-i)/5;
%         subplot(1,2,1)
%             hold all
%             plot(t,squeeze(g(14,1,i,:)),'-c','LineWidth',2,'Color',(w)*[0 1 1])
%             ylim([0,inf]); ylabel('RFP (au)')
%             xlabel('Time (h)'); box('on')
%             xticks(0:2.5:7.5);xlim([0,7.5]);
% 
%         subplot(1,2,2)
%             hold all
%             plot(t,squeeze(g(14,2,i,:)),'-r','LineWidth',2,'Color',(1-w)*[1 0 0])
%             ylim([0,inf]); ylabel('RFP (au)')
%             xlabel('Time (h)'); box('on')
%             xticks(0:2.5:7.5);xlim([0,7.5]);
% 
%     end