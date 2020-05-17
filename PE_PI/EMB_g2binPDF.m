function [pdf, maxG] = EMB_g2binPDF(g,b)
% EMB_BINPDF computes the joint conditional pdf p({g_i}|x) by binning
% input
% g: raw fluorescence values ,
% b: normalized bin size (0-1)
[N,I,X,T]=size(g);

% to normalize bin edges from 0 to max(max(g)), normalized to max(mean(g))
max_mean=max(mean(g,1),[],[3,4]);
max_max=max(g,[],[1,3,4]);

% bin number B
B=zeros(1,I);
edges=cell(1,I);
maxG=zeros(1,I);

for i=1:I
    edges{1,i}=[0:b:ceil(max_max(1,i)/max_mean(1,i)/b)*b]*max_mean(1,i);
    maxG(1,i)=edges{1,i}(end);
    B(1,i)=length(edges{1,i})-1;
end

% compute joint conditional pdf, uses histcn
% Bruno Luong (2020). N-dimensional histogram 
% (https://www.mathworks.com/matlabcentral/fileexchange/23897-n-dimensional-histogram),
% MATLAB Central File Exchange. Retrieved April 10, 2020.
pdf=zeros([B,X,T]);
for t=1:T
    for x=1:X
    	pdf(:,:,x,t)=histcn(g(:,:,x,t),edges{1,1},edges{1,2})/N; % not dynamic
    end
end
end

