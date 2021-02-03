function dx = EMB_degrade_v1(kdeg,Kdeg,x)
%EMB_DEGRADE_V1 Summary of this function goes here
%   Detailed explanation goes here
dx=-kdeg/(1+Kdeg/x);
end

