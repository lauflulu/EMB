clear all
close all

%%
load('data\data_full_1mM_v2.mat')
fusionTime=91;
g=EMB_data2g(data,fusionTime);
[N,I,X,T]=size(g);

T=1;

% bin sizes
%bin=linspace(0.02,0.05,10);
bin=linspace(0.02,0.5,10);
B=length(bin);
% fraction of samples M<N
m=ceil([0.95,0.9,0.85,0.8,0.75,0.5]*N);
M=length(m);
% number of iterations for each (bin,m)
K=100;

%all combinations of genes
perm= unique(nchoosek([zeros(1,I-1),1:I],I),'rows');
P=size(perm,1);
pi=zeros(K,M,B,P,T);
PI=zeros(P,T);

for k=1:K
    for n=1:M
        for b=1:B
            %draw=randi(N,m(n),1); %with replacement
            draw=randsample(N,m(n),false); %true/false=w/wo replacement
            h=g(draw,:,:,91);
            
            % shuffle pdf as sanity check, pdf_x and pdf_g stay as is
            for j=1:m(n)
                shuffle=randsample(X,X,false); 
                h(j,:,:,:)=h(j,:,shuffle,:);
            end
            
            nPDFg=EMB_g2binPDF(h,bin(b));

            for p=1:P
                % marginal/joint conditional distributions for each permutation
                % p({g_i}|x,t)
                condensedDim=1:I;
                condensedDim=condensedDim(1:I~=perm(p,perm(p,:)~=0));
                keptDim=perm(p,perm(p,:)~=0);
                if size(condensedDim,2)~=0
                    pdf=sum(nPDFg,condensedDim);
                else
                    pdf=nPDFg;
                end
                
                
                pdf=pdf./sum(pdf,1:I); %conditional pdf p({g_i}|x,t)
                pdf_g=sum(pdf,I+1)./sum(pdf,1:(I+1)); %marginal pdf_g(g)
                pdf_x=sum(pdf,1:I)./sum(pdf,1:(I+1)); % marginal pdf_x(x)
                pdf_xg=pdf./sum(pdf,I+1);pdf_xg(isnan(pdf_xg))=0; % conditional pdf p(x|{g_i},t)
                
                % eq. (8) Tkacik (2015)
                pi(k,n,b,p,:)= EMB_pdf2entropy(pdf_g,keptDim)-sum(EMB_pdf2entropy(pdf,keptDim).*pdf_x,3);
                % eq. (10) Tkacik (2015)
                %pi(k,n,b,p,:)= EMB_pdf2entropy(pdf_x,I+1)-sum(EMB_pdf2entropy(pdf_xg,I+1).*pdf_g,keptDim);
                % eq. (7) Tkacik (2015)
%                 int=zeros(5,T);
%                 for x=1:5
%                     l=pdf(:,:,x,:)./(pdf_g(:,:,1,:));l(isnan(l))=0;
%                     l=log2(l); 
%                     l(l==-inf)=0;
%                     int(x,:)=sum(pdf(:,:,x,:).*l,keptDim);
%                 end
%                 pi(k,n,b,p,:)=sum(squeeze(pdf_x) .* int,1);
 
            end
        end
    end
    k
end

%%

for p=1%:3
    for t=1%:91
        x=[ones(length(m),1),1./m'];
        y=squeeze(mean(pi(:,:,:,p,t),1));
        b1=x\y;
        y0=[1,0]*b1;
        z=[ones(length(bin),1),bin.^2'];
        b2=z\y0';
        PI(p,t)=[1,0]*b2;
    end
end


%%
figure(1)
    subplot(1,2,1)
        hold all
        plot(x(:,2),y)
        plot([0;x(:,2)],[1,0;x]*b1,'-r')
    subplot(1,2,2)
        hold all
        plot(z(:,2),y0)
        plot([0;z(:,2)],[1,0;z]*b2,'-r')
        
%%
% figure(2)
%     hold all
%     plot((0:90)*5/60,PI')
%     plot([0,7.5],[1,1]*log(5)/log(2),'-k')        
    