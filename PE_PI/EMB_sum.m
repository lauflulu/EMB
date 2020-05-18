function [output] = EMB_sum(input,dim)
%EMB_SUM Summary of this function goes here
%   Detailed explanation goes here
output=input;
D=length(dim);
for d=1:D
    output=sum(output,dim(d));
end

