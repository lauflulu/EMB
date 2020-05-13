function [g] = EMB_sim2g(y_R,genes)
%EMB_SIM2G rearranges simulation data into variable g  ('NIXT' format)
N=size(y_R{1,1},2);
I=size(genes,2);
X=5;
T=size(y_R{1,1},1);

g=zeros(N,I,X,T);

G=size(y_R,1)/6;

for i=1:I
    for x=1:X
        g(:,i,x,:)=y_R{genes(i)+x*G,1}';
    end
    g(:,i,:,:)=g(:,i,:,:)/max(mean(g(:,i,:,:),1),[],'all');
end

end

