function [meanG] = EMB_g2meanG(g)
%EMB_MEANG computes mean expression profiles
meanG=mean(g,1,'omitnan');
end

