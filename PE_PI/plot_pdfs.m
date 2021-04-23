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

%% calculations, cf. function descriptions for details
b=0.1; % bin size for binPDF in percent/100, scaled to max(meanG)
t=91; % time point for profile plots etc.
B=100; % number of bins for gaussPDF
Boot=1000; % number of bootstrap samples should be >=1000 for good quality
fusionTime=91;

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

%[PE,perm] = EMB_g2gaussPE(g);
%P=size(perm,1);
%[PI,perm] = EMB_g2cheatPI(g,B);

%[~,meanPE,stdPE] = EMB_bootstrap(@EMB_g2gaussPE,g,Boot);

%[~,meanPI,stdPI] = EMB_bootstrap(@(g)EMB_g2cheatPI(g,B),g,Boot); % this may take a while!

%% binPDFs single genes
figure(4)
    hy=subplot(1,4,2);
        cmap=zeros(256,3);
        for i=1:256
            w=(256-i)/255;
            cmap(i,:)=w*[1 1 1] + (1-w)*[0 210/255 188/255];
        end
        hold all
        imagesc(bPDFy,'XData', [1 5], 'YData', [0+b/2 maxGbin(1,1)-b/2])
        set(gca, 'Ydir', 'normal')
        colormap(hy,cmap)
        cby=colorbar;

        for i=1:5
            plot(nPDFy(:,i)+i,linspace(0,maxGgauss(1,1),B),'-b')
        end
        box('on');
        xlabel('Droplet Index');ylabel(cby,'p(g|x)')
        xlim([0.5,5.5]);ylim([0,1.4]);

    subplot(1,4,1)
        hold all
        pdf_g1=sum(bPDFy,2); pdf_g1=pdf_g1/sum(pdf_g1)*length(pdf_g1)/maxGbin(1,1);
        yBins=linspace(b/2,maxGbin(1,1)-b/2,size(binPDF,1));
        
        pdf_G1=sum(nPDFy,2); pdf_G1=pdf_G1/sum(pdf_G1)/maxGgauss(1,1)*B;
 
        barh(yBins,pdf_g1,1)
        plot(pdf_G1,linspace(0,maxGgauss(1,1),B))
        ylim([0,1.4])
        ylabel('FI (au)'); box('on');

    hr=subplot(1,4,4);
        cmap=zeros(256,3);
        for i=1:256
            w=(256-i)/255;
            cmap(i,:)=w*[1 1 1] + (1-w)*[255/255 76/255 0];
        end
        hold all
        imagesc(bPDFr,'XData', [1 5], 'YData', [0+b/2 maxGbin(1,2)-b/2])
        set(gca, 'Ydir', 'normal')
        colormap(hr,cmap)
        cbr=colorbar;
        for i=1:5
            plot(nPDFr(:,i)+i,linspace(0,maxGgauss(2,1),B),'-r')
        end
        
        box('on'); xlim([0.5,5.5]);ylim([0 1.8]);
        ylabel('FI (au)');xlabel('Droplet Index');
        ylabel(cbr,'p(g|x)')
        
    subplot(1,4,3)
        hold all
        pdf_g2=sum(bPDFr,2); pdf_g2=pdf_g2/sum(pdf_g2)*length(pdf_g2)/maxGbin(1,2);
        rBins=linspace(b/2,maxGbin(1,2)-b/2,size(binPDF,2));
        pdf_G2=sum(nPDFr,2); pdf_G2=pdf_G2/sum(pdf_G2)/maxGgauss(2,1)*B;
        barh(rBins,pdf_g2,1)
        plot(pdf_G2,linspace(0,maxGgauss(2,1),B))
        ylim([0,1.8])
        ylabel('FI (au)');box('on');
    
%% binPDF joint
figure(14)
    cmap=zeros(256,3);
    for i=1:256
        w=(256-i)/255;
        cmap(i,:)=w*[1 1 1] + (1-w)*[0 0 0];
    end
    contheights=[];
    for x=1:5
    hb=subplot(1,6,x+1);
        hold all
        bPDF=binPDF(:,:,x,91)*size(binPDF,1)*size(binPDF,2)/maxGbin(1,1)/maxGbin(1,2);
        imagesc(bPDF','XData', [0+b/2 maxGbin(1,1)-b/2], 'YData', [0+b/2 maxGbin(1,2)-b/2]);
        set(gca, 'Ydir', 'normal')
        
        g1=linspace(0,maxGgauss(2,1),B);g2=linspace(0,maxGgauss(1,1),B);
        [G1,G2] = meshgrid(g1,g2);
        gPDF=gaussPDF(:,:,x,91);
        sigma1 = mvnpdf([meanG(1,1,x,91)+1*covG(1,1,x,91)^0.5, meanG(1,2,x,91)],meanG(:,:,x,91),covG(:,:,x,91));
        sigma2 = mvnpdf([meanG(1,1,x,91)+2*covG(1,1,x,91)^0.5, meanG(1,2,x,91)],meanG(:,:,x,91),covG(:,:,x,91));
        sigma3 = mvnpdf([meanG(1,1,x,91)+3*covG(1,1,x,91)^0.5, meanG(1,2,x,91)],meanG(:,:,x,91),covG(:,:,x,91));
        contour(G2,G1,gPDF,[sigma1,sigma2,sigma3],'LineWidth',2,'LineColor','k');
        
        %plot(g(:,2,x,91),g(:,1,x,91),'or')
        
        colormap(hb,cmap);caxis([0,max(binPDF(:,:,:,91)*size(binPDF,1)*size(binPDF,2)/maxGbin(1,1)/maxGbin(1,2),[],'all')]);
        xlim([0 maxGbin(1,1)]);ylim([0 maxGbin(1,2)]);
        xticks(0:b:maxGbin(1,1)); yticks(0:b:maxGbin(1,2))
        box('on');
        xlabel('g1');
    end
    cby=colorbar;
    ylabel(cby,'p(g1,g2|x,t)')
    htot=subplot(1,6,1);
        hold all
        bPDF=sum(binPDF(:,:,:,91),3)*size(binPDF,1)*size(binPDF,2)/maxGbin(1,1)/maxGbin(1,2)/5;
        imagesc(bPDF','XData', [0+b/2 maxGbin(1,1)-b/2], 'YData', [0+b/2 maxGbin(1,2)-b/2]);
        set(gca, 'Ydir', 'normal')
        
       
        gPDF=sum(gaussPDF(:,:,:,91),3)/5;
        
        contour(G2,G1,gPDF,[0.02,0.2,2],'LineWidth',2,'LineColor','k');
        
        %plot(reshape(g(:,2,:,91),[],1),reshape(g(:,1,:,91),[],1),'or')
        colormap(htot,cmap);caxis([0,max(bPDF,[],'all')]);
        xlim([0 maxGbin(1,1)]);ylim([0 maxGbin(1,2)]);
        xticks(0:b:maxGbin(1,1)); yticks(0:b:maxGbin(1,2))
        box('on');
        xlabel('g1');ylabel('g2')
    
    cby=colorbar;
    ylabel(cby,'p(g1,g2|x,t)')
    
    
%% fit cumulative pdf
fusionTime=91;

load(datafiles{3});
g=EMB_data2g(data,fusionTime);

t=91;

fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0],...
               'StartPoint',[0.1,0.5]);
ft = fittype('0.5*(erf((x-b)./(sqrt(2).*a))+1)','options',fo);

% YFP
Xy=squeeze(sort(g(:,1,:,t),1)); % cumulative pdf
Yy=[1:size(g,1)]./size(g,1); Yy=Yy';

legy=cell(5,1);
xy=linspace(0,1.4,1000)';
yy=zeros(1000,5);
for i=1:5
    [fy,gofy]=fit(Xy(:,i),Yy,ft);
    legy{i,1}=sprintf('R^2 = %1.3f\nCVfit/CVs = %2.1f/%2.1f=%1.2f',...
        gofy.rsquare,fy.a/fy.b*100,std(Xy(:,i))/mean(Xy(:,i))*100,...
        fy.a/fy.b*100/(std(Xy(:,i))/mean(Xy(:,i))*100));
    yy(:,i)=fy(xy);
end

% RFP
Xr=squeeze(sort(g(:,2,:,t),1)); % cumulative pdf
Yr=[1:size(g,1)]./size(g,1); Yr=Yr';

legr=cell(5,1);
xr=linspace(0,1.8,1000)';
yr=zeros(1000,5);
for i=1:5
    [fr,gofr]=fit(Xr(:,i),Yr,ft);
    legr{i,1}=sprintf('R^2 = %1.3f\nCVfit/CVs = %2.1f/%2.1f=%1.2f',...
        gofr.rsquare,fr.a/fr.b*100,std(Xr(:,i))/mean(Xr(:,i))*100,...
        fr.a/fr.b*100/(std(Xr(:,i))/mean(Xr(:,i))*100));
    yr(:,i)=fr(xr);
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
        box('on');ylabel('cumulative pdf');xlabel('FI (au)'); ylim([-.05 1.05])
        legend(hy,legy)

    subplot(2,1,2)
        hold all
        for i=1:5
            w=(5-i)/4;
            col=w*[1 0 0]+(1-w)*[1 1 0];
            plot(Xr(:,i),Yr,'.','MarkerEdgeColor',col) 
            hr(i)=plot(xr,yr(:,i),'-','Color',col);
        end
        box('on');ylabel('cumulative pdf');xlabel('FI (au)');ylim([-.05 1.05])
        legend(hr,legr)