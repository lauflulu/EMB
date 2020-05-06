function [f,meanF,stdF] = EMB_bootstrap(fun,g,B)
%EMB_BOOTSTRAP generates distribution of output of fun, by bootstrapping g
[N,I,X,T]=size(g);
f=zeros([B,size(fun(g))]);

for b=1:B
    draw=randi(N,N,1);
    boot=g(draw,:,:,:);
	S.subs =  repmat({':'},1,ndims(f)); % the third row
    S.subs{1} = b;
	S.type = '()';
    [PE,~]=fun(boot);
	f = subsasgn(f,S,PE);
end
meanF=squeeze(mean(f,1));
stdF=squeeze(std(f,[],1));

end

