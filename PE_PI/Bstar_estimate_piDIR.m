clear all
close all

%%
load('data\data_full_1mM_v2.mat')
fusionTime=91;
g=EMB_data2g(data,fusionTime);
g=g(:,:,:,91);
[N,I,X,T]=size(g);

%T=1;

% bin sizes
%bin=linspace(0.02,0.05,10);
binNumber=linspace(3,8,10);
% fraction of samples M<N
m=ceil([0.95,0.9,0.85,0.8,0.75,0.5]*N);
% number of iterations for each (bin,m)
K=100;
% shuffle?
shuffle=true;


%% find max bin number B*
%m=5:2:N;
tol=0.1;
binNumber=[1:20];
%binNumber=[50,100];
B=length(binNumber);
PIb=zeros(B-2,3);
stdPIb=zeros(B-2,3);

for i=3:B
    [PIb(i-2,:),stdPIb(i-2,:)] = EMB_extrapolatePI(@EMB_g2piDIR,g,K,m,binNumber(i-2:i),shuffle);
    i
end
figure(3)
    errorbar(binNumber(3:B)'*ones(1,3),PIb,stdPIb,'--o')
      
bb=binNumber(2:B)'*ones(1,3);
for p=1:2
    Bmax(p)=max(bb((PIb(:,p)+stdPIb(:,p))<tol,p));
end

%%
% figure(2)
%     hold all
%     plot((0:90)*5/60,PI')
%     plot([0,7.5],[1,1]*log(5)/log(2),'-k')        
    