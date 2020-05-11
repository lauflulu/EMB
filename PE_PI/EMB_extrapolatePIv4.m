function [PI,stdPI] = EMB_extrapolatePIv4(fun,g,K,m,binNumber,shuffle,plotON)
%EMB_EXTRAPOLATEPI error bar is SD of MI estimate of smallest subsample
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
        end
    end
end

% extrapolate PI
PI=zeros(P,T);
stdPI=zeros(P,T);

conf=0.68;

x=1./m';
if plotON
    figure
end
for p=1:P
    D=sum(perm(p,:)~=0);
    y=1./binNumber.^(D+1)';
    meanI0=zeros(size(y));
    stdI0=zeros(size(y));
    for t=1:T
        
        meanI=squeeze(mean(pi(:,:,:,p,t),1));
        stdI=squeeze(std(pi(:,:,:,p,t),[],1));
        fMs=cell(B,1);
        
        if B>1
            for b=1:B
                [fM,gofM] = fit(x,meanI(:,b),'poly1');
                fMs{b,1}=fM;
                meanI0(b,1)=fM(0);
                stdI0(b,1) = max(stdI(:,b));
            end
            [fB,gofB] = fit(y,meanI0,'poly1');
            PI(p,t)=fB(0);
            stdPI(p,t) = max(stdI0);
            
        else
            [fM,gofM] = fit(x,meanI','poly1');
            fMs{1,1}=fM;
            PI(p,t)=fM(0);
            stdPI(p,t) = max(stdI);
        end
        
    end
    %plots
    if plotON
        
        if B>1
            subplot(P,2,2*p-1)
                hold all
                for b=1:B
                    errorbar(x,meanI(:,b),stdI(:,b),'--.','MarkerSize',20)
                    plot([0;x],fMs{b,1}([0;x]),'-k')
                end
                errorbar(zeros(B,1),meanI0,stdI0,'.r','MarkerSize',20)
                xlabel('1/m'); ylabel('PI naive')
                limx=max(x);
                box('on'); xlim([-0.05*limx,1.05*limx])
        
            subplot(P,2,2)
                hold all
                errorbar(y,meanI0,stdI0,'--.','MarkerSize',20)
                errorbar(0,PI(p,t),stdPI(p,t),'.b','MarkerSize',20);
                plot([0;y],fB([0;y]),'-k')
                xlabel('1/B^{D+1}'); ylabel('PI naive')
                limx=max(1./binNumber.^(2)');
                box('on');xlim([-0.05*limx,1.05*limx])
        else
                hold all
                errorbar(x,meanI,stdI,'--.','MarkerSize',20)
                plot([0;x],fMs{1,1}([0;x]),'-k')
                errorbar(zeros(B,1),PI(p,t),stdPI(p,t),'.r','MarkerSize',20)
                xlabel('1/M'); ylabel('PI naive');
                limx=max(x);
                xlim([-0.05*limx,1.05*limx])
                box('on')
        end
        
    end
end

end

