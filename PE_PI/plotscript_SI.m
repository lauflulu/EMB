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
    
%% PI DIR I_shuffle verification
% a B*=5 is a valid choice for estimating single gene PI for
% all data sets except 100 mM full (N=9) and 1 mM ctrl1 (N=12)
Bstar=zeros(10,3);

S=length(datafiles);
PI=cell(S,1);
stdPI=cell(S,1);

for sample=1:length(datafiles)
    sample
    
    load(datafiles{sample});

    fusionTime=91;
    g=EMB_data2g(data,fusionTime);
    g=g(:,:,:,91);
    [N,I,X,T]=size(g);
    
    if N>12
        m=flip(unique(ceil((6:11)/12*N)));
    else
        m=flip(7:N);
    end
    
    K=100;
    tol=0.1;
    binNumber=2:20;
    B=length(binNumber);
    PIb=zeros(B-2,3);
    stdPIb=zeros(B-2,3);

    for i=3:B
        [PIb(i-2,:),stdPIb(i-2,:)] = EMB_extrapolatePIv1(@EMB_g2piDIR,g,K,m,binNumber(1:i),true,false);
    end
    PI{S,1}=PIb; stdPI{S,1}=stdPIb;
    figure(3)
        subplot(5,2,sample)
        hold all
        errorbar(binNumber(3:B)'*ones(1,3),PIb,stdPIb,'--o')
        plot(binNumber,zeros(1,B),'--k')
        plot(binNumber,0.1*ones(1,B),'--k')
        ylim([-0.1,1.7]); box('on'); xlim([0,21]);
        xlabel('B*');ylabel('PI (bits)');
    
      
bb=binNumber(3:B)'*ones(1,3);
Bmax=bb.*((PIb(:,:)+stdPIb(:,:))>tol); 
Bmax(Bmax==0)=nan;
Bmax=min(Bmax,[],1)-1;
Bmax(isnan(Bmax))=21;
Bmax(Bmax==3)=nan;
Bstar(sample,:)=Bmax;
end
%% PI DIR/ SGA, extrapolation for t=91
% a B*=5 is a valid choice for estimating single gene PI for
% all data sets except 100 mM full (N=9) and 1 mM ctrl1 (N=12)
for sample=3%1:length(datafiles)
    load(datafiles{sample});

    fusionTime=91;
    g=EMB_data2g(data,fusionTime);
    g=g(:,:,:,[1:12:49,91]);
    [N,I,X,T]=size(g);
    m=ceil([0.95,0.9,0.85,0.8,0.75,0.5]*N);
    % number of iterations for each (bin,m)
    K=100;
    % shuffle?
    binNumber=2:6;

    [PIdir,stdPIdir] = EMB_extrapolatePI(@EMB_g2piDIR,g,K,m,binNumber,false,true);
    
    [PIsga,stdPIsga] = EMB_extrapolatePI(@EMB_g2piSGA,g,K,m,100,false,true);
    
    figure(11)
        hold all
        for p=1:3
            errorbar(PIsga(p,:),PIdir(p,:),...
                stdPIdir(p,:),stdPIdir(p,:),stdPIsga(p,:),stdPIsga(p,:),...
                '.-','MarkerSize',20)
        end
        %errorbar(PIdir,PIsga,stdPIsga,'.b','MarkerSize',20)
        plot(linspace(0,1.2),linspace(0,1.2),'--k')
        box('on'); xlim([-0.1,1.3]);ylim([-0.1,1.3]);
        xlabel('I_{SGA}');ylabel('I_{DIR}')
        legend('g1','g2','joint')
end

%% PI DIR/ SGA, time traces
% a B*=5 is a valid choice for estimating single gene PI for
% all data sets except 100 mM full (N=9) and 1 mM ctrl1 (N=12)
time=0:90; time=time*5/60;

for sample=3%1:length(datafiles)
    load(datafiles{sample});

    fusionTime=91;
    g=EMB_data2g(data,fusionTime);
    [N,I,X,T]=size(g);
    
    m=ceil([0.95,0.9,0.85,0.8,0.75,0.5]*N);
    % number of iterations for each (bin,m)
    K=100;
    binNumber=2:6;

    [PIdir,stdPIdir] = EMB_extrapolatePI(@EMB_g2piDIR,g,K,m,binNumber,false,false);
    
    [PIsga,stdPIsga] = EMB_extrapolatePI(@EMB_g2piSGA,g,K,m,100,false,false);
    
    figure(11)
        subplot(1,2,1)
            hold all
            plot(time, PIdir, '-');set(gca,'ColorOrderIndex',1)
            plot(time, PIdir+stdPIdir, '--');set(gca,'ColorOrderIndex',1)
            plot(time, PIdir-stdPIdir, '--');
            box('on'); xlim([0,7.5]);ylim([0,1.2]);
            xlabel('Time (h)');ylabel('I_{DIR}')
            legend('g1','g2','joint')
            xticks([0:2.5:7.5]);
            
        subplot(1,2,2)
            hold all
            plot(time, PIsga, '-');set(gca,'ColorOrderIndex',1)
            plot(time, PIsga+stdPIsga, '--');set(gca,'ColorOrderIndex',1)
            plot(time, PIsga-stdPIsga, '--');
            box('on'); xlim([0,7.5]);ylim([0,1.2]);
            xlabel('Time (h)');ylabel('I_{SGA}')
            legend('g1','g2','joint');xticks([0:2.5:7.5]);
end
%% PE
time=0:12:90; time=time*5/60;
load(datafiles{3});
fusionTime=91;
g=EMB_data2g(data,fusionTime);
g=g(:,:,:,1:12:91);
[N,I,X,T]=size(g);
m=ceil([0.95,0.9,0.85,0.8,0.75,0.7,0.6,0.5]*N);

[meanPE,stdPE] = EMB_extrapolatePE(@EMB_g2gaussPE,g,100,m,true);
%[PE,meanPE,stdPE] = EMB_bootstrap(@(g)EMB_extrapolatePE(@EMB_g2gaussPE,g,100,m,false),g,10);
%[PE,meanPE,stdPE] = EMB_bootstrap(@EMB_g2gaussPE,g,1000);
[Pcorr,lb,ub]=EMB_PE2Pcorr(meanPE,stdPE);

T=5;
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
        plot(1:5,Pcorr(1,:,T),'-.b','MarkerSize',20)
        plot(1:5,Pcorr(2,:,T),'-.r','MarkerSize',20)
        plot(1:5,Pcorr(3,:,T),'-.k','MarkerSize',20)
        plot(1:5,lb(1,:,T),'--b')
        plot(1:5,lb(2,:,T),'--r')
        plot(1:5,lb(3,:,T),'--k')
        plot(1:5,ub(1,:,T),'--b')
        plot(1:5,ub(2,:,T),'--r')
        plot(1:5,ub(3,:,T),'--k')
        xlim([0.5,5.5]);ylim([0,1]);box('on');
        xlabel('Droplet')
        ylabel('P_{corr}')
    
    for p=1:3
        subplot(2,3,3+p)
            plot(time,squeeze(Pcorr(p,:,:)))
            xlim([0,7.5]);xticks(0:2.5:7.5);
            box('on'); ylim([0,1]);
            xlabel('Time (h)')
            ylabel('P_{corr}')
    end
    


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
    
%% FI(t), with sender
load(datafiles{3});
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
figure(2)
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
figure(1);
    
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