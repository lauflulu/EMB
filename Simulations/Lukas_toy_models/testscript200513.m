clear all
close all

%% set parameters
N = 10; % number of samples
I = 8; % IPTG, free LacI, IPTG-LacI, RFP, RFPmat, TetR, YFP, YFPmat
X = 6; % 1 sender 5 receivers
T = 46;

% diffusion
%P = 1.67e-3*60; % um/min, for IPTG, determined from droplet fitting
P = 10e-3*60; % um/min, for IPTG, determined from droplet fitting

% geometry
cvV=0.35;
V = 1.2e6*(1+cvV*randn(N,X)); % um3 average volume of compartments, 1 nL
while min(V,[],'all')<0.5e6
    tmpV=1.2e6*(1+cvV*randn(N,X));
    V(V<0.5e6)=tmpV(V<0.5e6);
end
theta = 42*pi/180; % ?, interface angle
A = EMB_V2A(V, theta);

% gene expression
CValpha =0.25;
alphamax = 6e-3*(1+CValpha*randn(N,X)); % uM/min expression factor, around 100pM/s and up, from Schwarz-Schilling, Aufinger

aLacI = 3*1000/1077; % rNTPs length LacI, factor of expression normalized, and plasmid ratio
aRFP = 3*1000/693; % rNTPs length RFP
aTetR = 1000/618; % rNTPs length TetR
aYFP = 1000/717; % rNTPs length YFP
tau = 90; % min, cell-extract lifetime, decay
tmRNA = 15; % min, lifetime of RNA, from Karzbrun et al.
a0=0.1; % leak expression (fraction of max. expression)
alpha=[aLacI,aRFP,aTetR,aYFP,tau,tmRNA,a0];

% maturation
tmat1 = 33; % min, maturation time of mScarletI, from Balleza et al. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5765880/pdf/nihms916443.pdf
tmat2 = 12; % min, maturation time of YPet, from Balleza et al. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5765880/pdf/nihms916443.pdf
tmat=[tmat1,tmat2];

% Hill parameters
KdTetR = 0.4;%1.3e-1; % uM, determined from bulk fitting
nTetR = 1.3;%4.3; % Hill coefficient, TetR binding to pTetO, determined from bulk fitting
KdLacI = 0.1;%3.6e-2; % uM, LacI binding to pLacO, determined from bulk fitting
nLacI = 1.35; % Hill coefficient, LacI binding to pLacO, determined from bulk fitting
kHill=[KdTetR,nTetR,KdLacI,nLacI];

% IPTG induction
kon = 7.2; % /uM/min, on-rate of IPTG binding LacI, from Xu et al.
koff = 12.6; % /min, off-rate of IPTG binding LacI, from Xu et al.
k=[kon, koff];

% degradation
kdeg = 1e-3; % uM/min, rate of 0th-order protein degradation, from Karzbrun et al.
Kdeg = 1e-3; % uM, deducted from Karzbrun et al.
deg=[kdeg, Kdeg];

%% initialization
Y=zeros(N,I,X,T);

y0=zeros(I,X);
%initial values ~= 0
y0(1,1)=1*1e3; %mM to uM, IPTG sender
% create vector for odesolver
y0=reshape(y0,1,[]); y0=y0';

% duration of experiment
dt=10; tspan=0:dt:dt*(T-1);

options=odeset('AbsTol',1e-8,'RelTol',1e-4); % set tolerances

tic
for n=1:N
    [t,y]=ode23s(@(t,y)EMB_ODE1(t, y, I, X, A(n,:), V(n,:), P, alphamax(n,:), alpha, kHill, k, tmat, deg),tspan,y0,options);
    Y(n,:,:,:) = reshape(permute(y,[2,1]),I,X,T);
end
toc


%% PI
g = EMB_normalizeG(Y(:,[8,5],2:6,2:end));
time=t(2:end)/60;

meanY=squeeze(mean(Y(:,8,2:6,2:end),1));
meanR=squeeze(mean(Y(:,5,2:6,2:end),1));
stdY=squeeze(std(Y(:,8,2:6,2:end),[],1));
stdR=squeeze(std(Y(:,5,2:6,2:end),[],1));

[PIsga,perm] = EMB_g2piSGA(g,100); % without extrapolation since N=1000
[pdf,maxG]=EMB_g2binPDF(g,1/10);
%% plots

close all
figure(1)
    for x=1:5
    subplot(4,5,x)
        hold all
        plot(time, squeeze(Y(:,8,x+1,2:end)),'-k');
        plot(time, meanY(x,:), '-c','LineWidth',2);
        plot(time, meanY(x,:)+stdY(x,:), '--c','LineWidth',2);
        plot(time, meanY(x,:)-stdY(x,:), '--c','LineWidth',2);
        ylim([0,1.1*max(Y(:,8,:,2:end),[],'all')]); ylabel('YFPmat (au)');
        box('on'); 
        xlim([0,7.5]); xticks(0:2.5:7.5);
    subplot(4,5,5+x)
        hold all
        plot(time, squeeze(Y(:,5,x+1,2:end)),'-k')
        plot(time, meanR(x,:), '-r','LineWidth',2);
        plot(time, meanR(x,:)+stdR(x,:), '--r','LineWidth',2);
        plot(time, meanR(x,:)-stdR(x,:), '--r','LineWidth',2);
       ylim([0,1.1*max(Y(:,5,:,2:end),[],'all')]); ylabel('RFPmat (au)');
        box('on'); 
        xlim([0,7.5]); xticks(0:2.5:7.5);
    subplot(4,5,10+x)
        hold all
        
        plot([0,time(end)],[0.2,0.2],'--b','LineWidth',2)
        plot([0,time(end)],[2,2],'-b','LineWidth',2)
        plot([0,time(end)],[20,20],'--b','LineWidth',2)
        plot(time, squeeze(Y(:,1,x+1,2:end)),'-k')
        box('on'); ylabel('IPTG (mM)');
        xlim([0,7.5]); xticks(0:2.5:7.5);
    end
    subplot(4,1,4)
        hold all
        plot(time, PIsga', '-');set(gca,'ColorOrderIndex',1)
        %plot(time, PIsga+stdPIsga, '--');set(gca,'ColorOrderIndex',1)
        %plot(time, PIsga-stdPIsga, '--');
        box('on');% xlim([0,7.5]);ylim([-0.1,1.5]);
        xlabel('Time (h)');ylabel('I_{SGA}')
        legend('g1','g2','joint');%xticks([0:2.5:7.5]);
