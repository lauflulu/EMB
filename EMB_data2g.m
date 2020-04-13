function [g] = EMB_data2g(data,fusionTime)
% EMB_data2g computes a data tensor g from data struct
% samples with fusion events before timepoint T are excluded
% normalization by reference dye is performed
% dimensions: [N,I,X,T]=[Samples,genes,time (frame),space (droplet)]

% exclude samples starting with nan between frame 1 and fusion time
a={data.r}; a=cat(3, a{:});
excluder=ones(size(a,3),1);
for n=1:length(data)
    if ~(sum(sum(isnan(a(1:fusionTime,:,n))))==0)
        excluder(n,1)=0;
    end
end
data=data(excluder==1);

% shuffle around data
a={data.YFP}; a=cat(3, a{:});
b={data.RFP}; b=cat(3, b{:});
c={data.Cy5}; c=cat(3, c{:});

% normalization by reference dye
FIy=a./c; FIy=FIy(1:fusionTime,2:6,:); FIy=permute(FIy,[3 2 1]);
FIr=b./c; FIr=FIr(1:fusionTime,2:6,:); FIr=permute(FIr,[3 2 1]);

% write in g
N=sum(excluder); I=2; T=fusionTime; X=5;
g=zeros(N,I,X,T);
g(:,1,:,:)=FIy; g(:,2,:,:)=FIr;
end