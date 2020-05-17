function [Pcorr,stdPcorr] = EMB_PE2Pcorr(meanPE,stdPE)
%EMB_PE2PERR Summary of this function goes here
%   Detailed explanation goes here
[P,X,T]=size(meanPE);

Pcorr=zeros(P,X,T);
stdPcorr=zeros(P,X,T);

for p=1:P
    for x=1:X
        for t=1:T
            % analytical solutions (wolfram alpha)
            Pcorr(p,x,t)=erf(0.353553/meanPE(p,x,t));
            stdPcorr(p,x,t)= abs(0.398942*exp(-1/8./meanPE(p,x,t)^2)./meanPE(p,x,t)^2)*stdPE(p,x,t);
        end
    end
end

end
