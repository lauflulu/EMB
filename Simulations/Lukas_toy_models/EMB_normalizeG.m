function normG = EMB_normalizeG(g)
%EMB_NORMALIZEG normalizes g
[N,I,X,T]=size(g);
normG=zeros(size(g));
for i=1:I
    normG(:,i,:,:)=g(:,i,:,:)/max(mean(g(:,i,:,:),1),[],'all');
end
end

