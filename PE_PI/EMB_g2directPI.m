function [PI,perm] = EMB_g2directPI(g)
%EMB_G2DIRECTPI Summary of this function goes here
[N,I,X,T]=size(g);

% bin sizes
bin=linspace(0.02,0.05,10);
B=length(bin);
% fraction of samples M<N
m=ceil([0.95,0.9,0.85,0.8,0.75,0.5]*N);
M=length(m);

%all combinations of genes
perm= unique(nchoosek([zeros(1,I-1),1:I],I),'rows');
P=size(perm,1);
pi=zeros(M,B);
PI=zeros(P,T);

for n=1:M
    for b=1:B
        draw=randi(N,m(n),1);
        h=g(draw,:,:,:);
        nPDFg=EMB_g2binPDF(h,bin(b));

        for p=1:P
            % marginal/joint conditional distributions for each permutation
            % p({g_i}|x,t)
            condensedDim=1:I;
            condensedDim=condensedDim(1:I~=perm(p,perm(p,:)~=0));
            if size(condensedDim,2)~=0
                pdf=sum(nPDFg,condensedDim);
            else
                pdf=nPDFg;
            end
            pdf=pdf./sum(pdf,1:I); %conditional pdf p({g_i}|x,t)
            pdf_g=sum(pdf,I+1)./sum(pdf,1:(I+1)); %marginal pdf_g(g)
            pdf_x=sum(pdf,1:I)./sum(pdf,1:(I+1)); % marginal pdf_x(x)
            for t=1:T
                int_YR=zeros(5,1);
                for x=1:5
                    i_YR=log2( pdf(:,:,x,t)./(pdf_g(:,:,1,t))); 
                    i_YR(i_YR==inf)=0;
                    int_YR(x,1)=nansum(nansum(pdf(:,:,x,t).*i_YR));
                end
                PI(p,t)=nansum(squeeze(pdf_x(1,1,:,t)) .* int_YR);
            end
        end
    end
end

end

