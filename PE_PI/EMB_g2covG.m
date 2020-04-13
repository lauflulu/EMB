function [covG] = EMB_g2covG(g)
%EMB_MEANG computes covariance for each x, t
[N,I,X,T]=size(g);
covG=zeros(I,I,X,T);

for x=1:X
    for t=1:T
        covG(:,:,x,t)=cov(g(:,:,x,t));
    end
end
end

