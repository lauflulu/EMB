function [PI,stdPI] = EMB_extrapolatePE(fun,g,K,m,plotON)
%EMB_EXTRAPOLATEPI Summary of this function goes here
%   Detailed explanation goes here
M=size(m,2);
[N,I,X,T]=size(g);

%all combinations of genes
perm= unique(nchoosek([zeros(1,I-1),1:I],I),'rows');
P=size(perm,1);
pi=zeros(K,M,P,X,T);


for k=1:K
    for n=1:M
            draw=randsample(N,m(n),false); %true/false=w/wo replacement
            h=g(draw,:,:,:);
            [pi(k,n,:,:,:),~]=fun(h);
    end
end

% extrapolate PI
PI=zeros(P,X,T);
stdPI=zeros(P,X,T);

y=1./m';
if plotON
    figure
end
for p=1:P
    for x=1:X
        for t=1:T
            meanI=squeeze(mean(pi(:,:,p,x,t),1));
            stdI=squeeze(std(pi(:,:,p,x,t),[],1));

            [fM,~] = fit(y,meanI','poly1');
            ci=confint(fM,0.68);
            PI(p,x,t)=fM(0);
            stdPI(p,x,t) = ci(2,2)-fM(0);    
        end
        if plotON
            subplot(1,5,x)
            hold all
            plot([0;y],fM([0;y]),'-k')
            errorbar(y,meanI,stdI,'--.','MarkerSize',20)
            errorbar(0,PI(p,x,t),stdPI(p,x,t),'.r','MarkerSize',20)
            xlabel('1/M'); ylabel('PE (dropelts)');
            limx=max(y);
            xlim([-0.05*limx,1.05*limx])
            box('on')
        end
    end
end

end

