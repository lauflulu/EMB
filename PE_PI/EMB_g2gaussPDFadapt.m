function [gaussPDF,maxG] = EMB_g2gaussPDFadapt(g,B)
%EMB_GAUSSPDF infers normal pdf from mean + cov of expression profiles
[N,I,X,T]=size(g);

meanG=EMB_g2meanG(g);
covG=EMB_g2covG(g);

maxG=zeros(I,1);
for i=1:I
    maxG(i,1)=max(max(meanG(:,i,:,:),[],3),[],4)+5*max(max(covG(i,i,:,:),[],3),[],4)^0.5;
end
B1=ceil(B*maxG(1,1)); B2=ceil(B*maxG(2,1));
[G1,G2]=ndgrid(linspace(0,maxG(1,1),B1),linspace(0,maxG(2,1),B2));
G = [G1(:) G2(:)];

gaussPDF=zeros(B1,B2,X,T);

for x=1:X
    for t=1:T
        gaussPDF(:,:,x,t)=reshape(mvnpdf(G,squeeze(meanG(:,:,x,t)),covG(:,:,x,t)),B1,B2);
    end
end

end


