function [PI,perm] = EMB_g2cheatPI(g,B)
%EMB_g2CHEATPI calculates PI using inferred normal PDFs
[N,I,X,T]=size(g);
nPDFg=EMB_g2gaussPDF(g,B);

%all combinations of genes
perm= unique(nchoosek([zeros(1,I-1),1:I],I),'rows');
P=size(perm,1);
PI=zeros(P,T);

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
    pdf=pdf./sum(pdf,1:I); %conditional pdf
    pdf_g=sum(pdf,I+1)./sum(pdf,1:(I+1)); %marginal pdf g
    pdf_x=sum(pdf,1:I)./sum(pdf,1:(I+1)); % marginal pdf x
    for t=1:T
        int_YR=zeros(5,1);
        for x=1:5
            i_YR=log2( pdf(:,:,x,t)./(pdf_g(:,:,1,t))); 
            i_YR(i_YR==inf)=0;
            int_YR(x,1)=nansum(nansum(pdf(:,:,x,t)/pdf_x(1,1,x,t).*i_YR));
        end
        PI(p,t)=nansum(squeeze(pdf_x(1,1,:,t)) .* int_YR);
    end
end

end
    

