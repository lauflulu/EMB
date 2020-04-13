clear all
close all
%%
datafiles={'190925_full_100mM\data_full_100mM.mat',...
        '190927_full_10mM\data_full_10mM.mat',...
        '190928_full_1mM\data_full_1mM.mat',...
        '190929_full_100uM\data_full_100uM.mat',...
        '190930_ctrl1_100mM\data_ctrl1_100mM.mat',...
        '190930_ctrl1_10mM\data_ctrl1_10mM.mat',...
        '191002_ctrl1_1mM\data_ctrl1_1mM.mat',...
        '191002_ctrl1_100uM\data_ctrl1_100uM.mat',...
        '191010_ctrl2_100mM\data_ctrl2_100mM.mat',...
        '191010_ctrl2_10mM\data_ctrl2_10mM.mat'};  
    
%%
T_filter=91;
x=linspace(0,0.3,1000);
bootX=10000;
T=97;


for p=3%1:size(datafiles,2)
    load(datafiles{p});
    a={data.YFP};
    b={data.RFP};
    c={data.Cy5};
    d={data.dataSize};
    a=cat(3, a{:});
    b=cat(3, b{:});
    c=cat(3, c{:});
    d=squeeze(cat(3, d{:}));
    d=(min(d,[],1)>=T_filter);
    
    FIy=a./c;
    FIy=FIy(:,2:6,d);
    meanYFP=nanmean(FIy,3);
    stdYFP=nanstd(FIy,[],3);
    FIr=b./c;
    FIr=FIr(:,2:6,d);
     
    bootN=size(FIr,3);
    bootRFP=zeros(bootN,size(FIr,2));
    bootYFP=zeros(bootN,size(FIr,2));
    
    I_Y=zeros(1,90,bootX);
    J_Y=zeros(5,90,bootX);
    I_R=zeros(1,90,bootX);
    J_R=zeros(5,90,bootX);
    
    
    for t=1:T
        t
        for j=1:bootX
            %rfp
            draw=randi(size(FIr,3),bootN,1);
            for jj=1:bootN
                bootRFP(jj,:)=squeeze(FIr(t,:,draw(jj,1)))';
            end
            mu_R=nanmean(bootRFP,1); sigma_R=nanstd(bootRFP,[],1);
            for i=1:5
                pdf_R(:,i)=normpdf(x,mu_R(1,i),sigma_R(1,i));
            end
            pdf_R=pdf_R/nansum(nansum( pdf_R));
            J_R(:,t,j)=ddx_f(mu_R).^2./sigma_R.^2 + ...
                            2*ddx_f(sigma_R).^2./sigma_R.^2;

            iR=log2( pdf_R./(sum(pdf_R,2)*sum(pdf_R,1))); 
            iR(iR==inf)=0;
            I_R(1,t,j)=nansum(nansum( pdf_R .* iR));
            
            %yfp
            for jj=1:bootN
                bootYFP(jj,:)=squeeze(FIy(t,:,draw(jj,1)))';
            end
            mu_Y=nanmean(bootYFP,1); sigma_Y=nanstd(bootYFP,[],1);
            for i=1:5
                pdf_Y(:,i)=normpdf(x,mu_Y(1,i),sigma_Y(1,i));
            end
            pdf_Y=pdf_Y/nansum(nansum( pdf_Y));
            J_Y(:,t,j)=ddx_f(mu_Y).^2./sigma_Y.^2 + ...
                            2*ddx_f(sigma_Y).^2./sigma_Y.^2;

            iY=log2( pdf_Y./(sum(pdf_Y,2)*sum(pdf_Y,1))); 
            iY(iY==inf)=0;
            I_Y(1,t,j)=nansum(nansum( pdf_Y .* iY));  
        end
    end

%%
figure(2);
    %subplot(5,2,p)
        hold all
        plot([0:T]*5/60,mean(I_Y,3),'-c','LineWidth',2)
        plot([0:T]*5/60,mean(I_Y,3)+std(I_Y,[],3),'--c','LineWidth',2)
        plot([0:T]*5/60,mean(I_Y,3)-std(I_Y,[],3),'--c','LineWidth',2)
        plot([0:T]*5/60,mean(I_R,3),'-m','LineWidth',2)
        plot([0:T]*5/60,mean(I_R,3)+std(I_R,[],3),'--m','LineWidth',2)
        plot([0:T]*5/60,mean(I_R,3)-std(I_R,[],3),'--m','LineWidth',2)
        box('on')
        xlabel('Time (h)')
        ylabel('PI (bits)')
        ylim([0,1.5]); xlim([0,8]);
        header=split(datafiles{p},'\');
        title([header{1},sprintf(', N = %d',sum(d))],'Interpreter','none')
end
%%
function ddxpdf=ddx_f(f)
    [~,Lx]=size(f);
    leftEdge=f(:,2)-f(:,1);
    rightEdge=f(:,Lx)-f(:,Lx-1);
    center = f(:,3:Lx)/2 - f(:,1:Lx-2)/2; %central derivative
    ddxpdf =[leftEdge,center,rightEdge];
end