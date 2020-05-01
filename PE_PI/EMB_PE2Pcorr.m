function [Perr,lbPerr,ubPerr] = EMB_PE2Pcorr(meanPE,stdPE)
%EMB_PE2PERR Summary of this function goes here
%   Detailed explanation goes here
[P,X,T]=size(meanPE);

Perr=zeros(P,X,T);
lbPerr=zeros(P,X,T);
ubPerr=zeros(P,X,T);
d=0.001;
D=0:d:1;
for p=1:P
    for x=1:X
        for t=1:T
            Perr(p,x,t)=1-sum(normpdf(D,0.5,meanPE(p,x,t)))*d;
            lbPerr(p,x,t)=1-sum(normpdf(D,0.5,meanPE(p,x,t)-stdPE(p,x,t)))*d;
            ubPerr(p,x,t)=1-sum(normpdf(D,0.5,meanPE(p,x,t)+stdPE(p,x,t)))*d;
        end
    end
end

end
