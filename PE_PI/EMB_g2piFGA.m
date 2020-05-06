function [pi,perm] = EMB_g2piFGA(g,B)
% EMB_g2piDE computes naive first Gaussian assumption estimate for PI 
% g: M<=N subsamples of data
% B: number of bins in the interval 0-1, where 1 is normalized to max(g_i)
[N,I,X,T]=size(g);

%all combinations of genes
perm= unique(nchoosek([zeros(1,I-1),1:I],I),'rows');
P=size(perm,1);
pi=zeros(P,T);

nPDFg=EMB_g2binPDF(g,1/B);
[GaussPDFg,~]=EMB_g2gaussPDF(g,B);

for p=1:P
    % marginal/joint conditional distributions for each permutation
    % p({g_i}|x,t)
    condensedDim=1:I;
    condensedDim=condensedDim(1:I~=perm(p,perm(p,:)~=0));
    keptDim=perm(p,perm(p,:)~=0); D=size(keptDim,2);
    if size(condensedDim,2)~=0
        pdf=sum(nPDFg,condensedDim);
        GaussPDF=sum(GaussPDFg,condensedDim);
    else
        pdf=nPDFg;
        GaussPDF=GaussPDFg;
    end

    pdf_gx=pdf./sum(pdf,1:I); %conditional pdf p({g_i}|x,t)
    pdf_g=sum(pdf,I+1)./sum(pdf,1:(I+1)); %marginal pdf_g(g|t)
    pdf_x=sum(pdf,1:I)./sum(pdf,1:(I+1)); % marginal pdf_x(x|t)
    pdf_xg=pdf./sum(pdf,I+1);pdf_xg(isnan(pdf_xg))=0; % conditional pdf p(x|{g_i},t)
    
    % eq. (8) Tkacik (2015)
    totalS=EMB_pdf2entropy(pdf_g,keptDim);
    
%     noiseS=sum(EMB_pdf2entropy(GaussPDF./sum(GaussPDF,1:I),keptDim).*pdf_x,3);
    % eq. 25 something is wrong with the normalization
    C=EMB_g2covG(g);
    noiseS=zeros(size(totalS));
    for t=1:T
        for x=1:X
            noiseS(:,:,x,t)=0.5*log2((2*2*asin(1)*exp(1))^D*det(C(:,:,x,t)))-D*log2(B);
        end
    end
%     noiseS=sum(noiseS.*pdf_x,3);
    pi(p,:)= totalS-noiseS;
    % eq. (10) Tkacik (2015)
    %pi(p,:)= EMB_pdf2entropy(pdf_x,I+1)-sum(EMB_pdf2entropy(pdf_xg,I+1).*pdf_g,keptDim);
    % eq. (7) Tkacik (2015)
%                 int=zeros(5,T);
%                 for x=1:5
%                     l=pdf(:,:,x,:)./(pdf_g(:,:,1,:));l(isnan(l))=0;
%                     l=log2(l); 
%                     l(l==-inf)=0;
%                     int(x,:)=sum(pdf(:,:,x,:).*l,keptDim);
%                 end
%                 pi(p,:)=sum(squeeze(pdf_x) .* int,1);
end

end

