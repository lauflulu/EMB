function [PI,stdPI] = EMB_extrapolatePIv3(fun,g,K,m,binNumber,shuffle,plotON)
%EMB_EXTRAPOLATEPI error bars are SD of K individual fits
M=size(m,2);
B=size(binNumber,2);
[N,I,X,T]=size(g);

%all combinations of genes
perm= unique(nchoosek([zeros(1,I-1),1:I],I),'rows');
P=size(perm,1);
pi=zeros(K,M,B,P,T);


for k=1:K
    for n=1:M
        for b=1:B
            draw=randsample(N,m(n),false); %true/false=w/wo replacement
            h=g(draw,:,:,:);
            
            % shuffle pdf as sanity check, pdf_x and pdf_g stay as is
            if shuffle
                for j=1:m(n)
                    shuffle=randsample(X,X,false); 
                    h(j,:,:,:)=h(j,:,shuffle,:);
                end
            end
            pi(k,n,b,:,:)=fun(h,binNumber(b));
%               pi(k,n,b,:,:)=EMB_g2piDIR(h,bin(b));
        %   pi(k,n,b,:,:)=EMB_g2piFGA(h,bin(b));
        end
    end
end

% extrapolate PI
PI=zeros(P,T);
stdPI=zeros(P,T);

x=1./m';
if plotON
    figure
end
for p=1:P
    D=sum(perm(p,:)~=0);
    y=1./binNumber.^(D+1)';
    for t=1:T
        fMs=cell(B,K);
        fBs=cell(1,K);
        
        if B>1
            piK=zeros(K,1);
            piB=zeros(B,K);
            for k=1:K
                
                for b=1:B
                    [fM,~] = fit(x,squeeze(pi(k,:,b,p,t))','poly1');
                    piB(b,k)=fM(0);
                    fMs{b,k}=fM;
                end      
                [fB,~] = fit(y,piB(:,k),'poly1');
                piK(k,1)=fB(0);
                fBs{1,k}=fB;
            end
            PI(p,t)=mean(piK,1);
            stdPI(p,t) = std(piK,[],1);
            
        else
            piK=zeros(K,1);
            for k=1:K
                [fM,~] = fit(x,squeeze(pi(k,:,:,p,t))','poly1');
                piK(k,1)=fM(0);
                fMs{1,k}=fM;
            end
 
            PI(p,t)=mean(piK,1);
            stdPI(p,t) = std(piK,[],1);
        end
        
    end
    %plots
    if plotON
        if B>1
            subplot(P,2,2*p-1)
                hold all
                for b=1:B
                    line=zeros(K,length(x)+1);
                    for k=1:K
                        hl=plot([0;x],fMs{b,k}([0;x]),'-k');
                        hl.Color(4) = 0.2;
                        line(k,:)=fMs{b,k}([0;x]);
                    end
                    set(gca,'ColorOrderIndex',b)
                    plot([0;min(x)],mean(line(:,1:2),1),'--','LineWidth',2); set(gca,'ColorOrderIndex',b)
                    plot(x,mean(line(:,2:end),1),'-','LineWidth',2); set(gca,'ColorOrderIndex',b)
                    errorbar(x,squeeze(mean(pi(:,:,b,p,t),1)),squeeze(std(pi(:,:,b,p,t),[],1)),'.','MarkerSize',20)
                end
                errorbar(zeros(B,1),mean(piB,2),std(piB,[],2),'.r','MarkerSize',20)
                xlabel('1/m'); ylabel('PI naive')
                limx=max(x);
                box('on'); xlim([-0.05*limx,1.05*limx])
        
            subplot(P,2,2)
                hold all
                line=zeros(K,length(y)+1);
                for k=1:K
                    hl=plot([0;y],fBs{1,k}([0;y]),'-k');
                    hl.Color(4) = 0.2;
                    line(k,:)=fBs{1,k}([0;y]);
                end
                set(gca,'ColorOrderIndex',p)
                plot([0;min(y)],mean(line(:,1:2),1),'--','LineWidth',2);set(gca,'ColorOrderIndex',p)
                plot(y,mean(line(:,2:end),1),'-','LineWidth',2);set(gca,'ColorOrderIndex',p)
                errorbar(y,mean(piB,2),std(piB,[],2),'.','MarkerSize',20)
                
                errorbar(0,PI(p,t),stdPI(p,t),'.b','MarkerSize',20)
                xlabel('1/B^{D+1}'); ylabel('PI naive')
                limx=max(1./binNumber.^(2)');
                box('on');xlim([-0.05*limx,1.05*limx])
        else
                hold all
                line=zeros(K,length(x)+1);
                for k=1:K
                    hl=plot([0;x],fMs{1,k}([0;x]),'-k');
                    hl.Color(4) = 0.2;
                    line(k,:)=fMs{1,k}([0;x]);
                end
                errorbar(x,squeeze(mean(pi(:,:,:,p,t),1))',squeeze(std(pi(:,:,:,p,t),[],1))','.','MarkerSize',20)
                plot([0;min(x)],mean(line(:,1:2),1),'--r','LineWidth',2)
                plot(x,mean(line(:,2:end),1),'-r','LineWidth',2)
                errorbar(zeros(B,1),PI(p,t),stdPI(p,t),'.r','MarkerSize',20)
                xlabel('1/M'); ylabel('PI naive');
                limx=max(x);
                xlim([-0.05*limx,1.05*limx])
                box('on')
        end
        
    end
end

end

