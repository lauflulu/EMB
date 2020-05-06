x=[0:1:10]';
y=x+2*randn(length(x),1);

[fM,gof,out] = fit(x,y,'poly1');



figure(1)
    hold all
    plot(x,y)
    plot(fM)
    
    
%%
conf=[0.51:0.01:0.95];
C= length(conf);
ci=zeros(1,C);

for c=1:C
    CI=confint(fM,conf(c));
    ci(1,c)=CI(2,2)-fM(0);
%(ci(2,2)-fM(0))/tinv(conf,length(x)-2)*(length(x)-2)^0.5;
end

figure(2)
    subplot(1,2,1)
    hold all
    plot(conf,ci)
    subplot(1,2,2)
    hold all
    plot(conf,ci./tinv((1-(1-conf)/2),length(x)-2))
    %%
    % generate data
x = 0:.1:10;
y = x.*x + randn(size(x));
w = linspace(0.01, 0.021,length(x));
x = x(:);
y = y(:);
w = w(:);
%plot data
plot(x,y,'.');
%fit
ft = fittype('poly2');
cf = fit(x,y,ft,'Weights',w);
% Plot fit
hold on
plot(cf,'fit',0.95);
ci=confint(cf);
ci(2,3)-ci(1,3)