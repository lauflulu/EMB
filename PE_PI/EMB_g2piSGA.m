function [pi,perm] = EMB_g2piSGA(g,B)
%EMB_g2CHEATPI calculates PI using inferred normal PDFs
[N,I,X,T]=size(g);
nPDFg=EMB_g2gaussPDFadapt(g,B);

%all combinations of genes
perm= unique(nchoosek([zeros(1,I-1),1:I],I),'rows');
P=size(perm,1);
pi=zeros(P,T);

for p=1:P
    % marginal/joint conditional distributions for each permutation
    % p({g_i}|x,t)
    condensedDim=1:I;
    condensedDim=condensedDim(1:I~=perm(p,perm(p,:)~=0));
    keptDim=perm(p,perm(p,:)~=0); D=size(keptDim,2);
    if size(condensedDim,2)~=0
        pdf=EMB_sum(nPDFg,condensedDim);
    else
        pdf=nPDFg;
    end
    pdf_gx=pdf./EMB_sum(pdf,1:I); %conditional pdf p({g_i}|x,t)
    pdf_g=EMB_sum(pdf,I+1)./EMB_sum(pdf,1:(I+1)); %marginal pdf_g(g)
    pdf_x=EMB_sum(pdf,1:I)./EMB_sum(pdf,1:(I+1)); % marginal pdf_x(x)
    
    totalS=EMB_pdf2entropy(pdf_g,keptDim);
    noiseS=EMB_sum(EMB_pdf2entropy(pdf_gx,keptDim).*pdf_x,3);
    pi(p,:)= totalS-noiseS;
end

end
    

