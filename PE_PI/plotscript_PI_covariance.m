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
T=90;
BinX=1000;
g=linspace(0,0.3,BinX);

for p=1:size(datafiles,2)
    load(datafiles{p});
    a={data.YFP};
    b={data.RFP};
    c={data.Cy5};
    d={data.dataSize};
    a=cat(3, a{:});
    b=cat(3, b{:});
    c=cat(3, c{:});
    d=squeeze(cat(3, d{:}));
    d=(min(d,[],1)>T);
    N=sum(d);
    
    FIy=a./c;
    FIy=FIy(:,2:6,d);
    meanYFP=nanmean(FIy,3);
    stdYFP=nanstd(FIy,[],3);
    FIr=b./c;
    FIr=FIr(:,2:6,d);
    meanRFP=nanmean(FIr,3);
    stdRFP=nanstd(FIr,[],3);
    
    I_Y=zeros(1,90);
    J_Y=zeros(5,90);
    I_R=zeros(1,90);
    J_R=zeros(5,90);
    
    J_YR=zeros(5,90);
    I_YR=zeros(1,90);
    
    for t=1:90
        mu_R=meanRFP(t,:); sigma_R=stdRFP(t,:);
        for x=1:5
            pdf_R(:,x)=normpdf(g,mu_R(1,x),sigma_R(1,x));
        end
        pdf_R=pdf_R/nansum(nansum( pdf_R));
        J_R(:,t)=ddx_f(mu_R).^2./sigma_R.^2 + ...
                        2*ddx_f(sigma_R).^2./sigma_R.^2;
        
        iR=log2( pdf_R./(sum(pdf_R,2)*sum(pdf_R,1))); 
        iR(iR==inf)=0;
        I_R(1,t)=nansum(nansum( pdf_R .* iR));
        
        mu_Y=meanYFP(t,:); sigma_Y=stdYFP(t,:);
        for i=1:5
            pdf_Y(:,i)=normpdf(g,mu_Y(1,i),sigma_Y(1,i));
        end
         pdf_Y=pdf_Y/nansum(nansum( pdf_Y));
        J_Y(:,t)=ddx_f(mu_Y).^2./sigma_Y.^2 + ...
                        2*ddx_f(sigma_Y).^2./sigma_Y.^2;
        
        iY=log2( pdf_Y./(sum(pdf_Y,2)*sum(pdf_Y,1))); 
        iY(iY==inf)=0;
        I_Y(1,t)=nansum(nansum(pdf_Y .* iY));
        %covariance matrix
        
        FI_YR=zeros(5,N,2);
        FI_YR(:,:,1)=FIy(t,:,:); FI_YR(:,:,2)=FIr(t,:,:);
        
%         pdf_YR=zeros(BinX,5,2);
%         pdf_YR(:,:,1)=pdf_Y; pdf_YR(:,:,2)=pdf_R;
        
        mu_YR=zeros(5,2);
        mu_YR(:,1)=mu_Y;mu_YR(:,2)=mu_R;
        
%         sigma_YR=zeros(1,5,2);
%         sigma_YR(:,:,1)=sigma_Y;sigma_YR(:,:,2)=sigma_R;

        C=zeros(2,2,5);
        for x=1:5
            for i=1:2
                for j=1:2
                    C(i,j,x) = (nanmean(FI_YR(x,:,i).*FI_YR(x,:,j)) ...
                             - mu_YR(x,i).*mu_YR(x,j))*N/(N-1);
                end
            end
            
        end
        dgdx=ddx_f(mu_YR');
        dCdx=zeros(size(C));
        for i=1:2
            for j=1:2
                dCdx(i,j,:)=ddx_f(C(i,j,:));
            end
        end
        for x=1:5
            J_YR(x,t)=dgdx(:,x)'*inv(squeeze(C(:,:,x)))*dgdx(:,x) ...
                     + 0.5*trace(inv(C(:,:,x))*dCdx(:,:,x)*inv(C(:,:,x))*dCdx(:,:,x));
        end
        % bivariate conditional pdf
        pdf_YR=zeros(BinX,BinX,5);
        [G1,G2] = meshgrid(g,g);
        G = [G1(:) G2(:)];
        for x=1:5
            pdf_YR(:,:,x)=reshape(mvnpdf(G,mu_YR(x,:),C(:,:,x)),BinX,BinX)'; %rows Y, columns R
        end
        pdf_YR=pdf_YR/nansum(nansum(nansum( pdf_YR)));
        pdf_x=sum(sum(pdf_YR,1),2);
        int_YR=zeros(5,1);
        for x=1:5
            i_YR=log2( pdf_YR(:,:,x)./(sum(pdf_YR,3)*pdf_x(1,1,x))); 
            i_YR(i_YR==inf)=0;
            int_YR(x,1)=nansum(nansum(pdf_YR(:,:,x)/pdf_x(1,1,x).*i_YR));
        end
        I_YR(1,t)=nansum(squeeze(pdf_x) .* int_YR);
    end
    
    %% plots
    sigmaY_x=sqrt(1./J_Y);
    sigmaR_x=sqrt(1./J_R);
    sigmaYR_x=sqrt(1./J_YR);

%     figure(1);
%     %subplot(5,2,p)
%         hold all
%         plot(sigmaY_x(:,49),'o-b','LineWidth',2)
%         plot(sigmaR_x(:,49),'o-r','LineWidth',2)
%         plot(sigmaYR_x(:,49),'o-g','LineWidth',2)
%         xlabel('Droplet')
%         ylabel('PE (droplets)')
%         %ylim([0,1.5]);
%         header=split(datafiles{p},'\');
%         title([header{1},sprintf(', N = %d',sum(d))],'Interpreter','none')


    figure(2);
    subplot(5,2,p)
        hold all
        plot([1:90]*5/60,I_Y,'o-b','LineWidth',2)
        plot([1:90]*5/60,I_R,'o-r','LineWidth',2)
        plot([1:90]*5/60,I_YR,'o-g','LineWidth',2)

        xlabel('Time (h)')
        ylabel('PI (bits)')
        ylim([0,1.5]);
        header=split(datafiles{p},'\');
        title([header{1},sprintf(', N = %d',sum(d))],'Interpreter','none')
    
    
    
    
end

%% testing only
%

function ddxpdf=ddx_f(f)
    [~,Lx]=size(f);
    leftEdge=f(:,2)-f(:,1);
    rightEdge=f(:,Lx)-f(:,Lx-1);
    center = f(:,3:Lx)/2 - f(:,1:Lx-2)/2; %central derivative
    ddxpdf =[leftEdge,center,rightEdge];
end