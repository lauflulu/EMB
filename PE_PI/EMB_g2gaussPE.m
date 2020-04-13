function [PE,perm] = EMB_g2gaussPE(g)
%EMB_GAUSSPE computes PE from meanG and covG
% initialize
[N,I,X,T]=size(g);
meanG=EMB_g2meanG(g);
covG=EMB_g2covG(g);

%all combinations of genes
perm= unique(nchoosek([zeros(1,I-1),1:I],I),'rows');
P=size(perm,1);
J=zeros(P,X,T);
PE=zeros(P,X,T);

for p=1:P
    C=covG(perm(p,perm(p,:)~=0),perm(p,perm(p,:)~=0),:,:); %subset covariance matrix
    m=meanG(1,perm(p,perm(p,:)~=0),:,:);
    for t=1:T
        % derivative, element-wise
        dCdx=zeros(size(C));
        for i=1:sum(perm(p,:)~=0)
            for j=1:sum(perm(p,:)~=0)
                dCdx(i,j,:,t)=EMB_dfdx(C(i,j,:,t));
            end
        end
        dgdx=zeros(size(m));
        for i=1:sum(perm(p,:)~=0)
            dgdx(1,i,:,t)=EMB_dfdx(m(1,i,:,t));
        end
        for x=1:X
            J(p,x,t)=squeeze(dgdx(:,:,x,t))*inv(squeeze(C(:,:,x,t)))*squeeze(dgdx(:,:,x,t))' ...
                     + 0.5*trace(inv(C(:,:,x,t))*dCdx(:,:,x,t)*inv(C(:,:,x,t))*dCdx(:,:,x,t));
            PE(p,x,t)=(1./J(p,x,t)).^0.5;
        end
    end
end

end

