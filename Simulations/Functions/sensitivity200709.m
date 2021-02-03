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

% maturation
tmat1 = 33; % min, maturation time of mScarletI, from Balleza et al. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5765880/pdf/nihms916443.pdf
tmat2 = 12; % min, maturation time of YPet, from Balleza et al. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5765880/pdf/nihms916443.pdf

% Hill parameters
KdTetR = 0.4;%1.3e-1; % uM, determined from bulk fitting
nTetR = 1.3;%4.3; % Hill coefficient, TetR binding to pTetO, determined from bulk fitting
KdLacI = 0.1;%3.6e-2; % uM, LacI binding to pLacO, determined from bulk fitting
nLacI = 1.35; % Hill coefficient, LacI binding to pLacO, determined from bulk fitting

% IPTG induction
kon = 7.2; % /uM/min, on-rate of IPTG binding LacI, from Xu et al.
koff = 12.6; % /min, off-rate of IPTG binding LacI, from Xu et al.

% degradation
kdeg = 1e-3; % uM/min, rate of 0th-order protein degradation, from Karzbrun et al.
Kdeg = 1e-3; % uM, deducted from Karzbrun et al.

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

%% simulation
tic

params=[P, ...
    aLacI, aRFP, aTetR, aYFP, tau, tmRNA, a0, ...
    KdTetR, nTetR, KdLacI, nLacI, ...
    tmat1, tmat2, ...
    kon, koff, ...
    kdeg, Kdeg];

paraNames={'P', ...
    'aLacI', 'aRFP', 'aTetR', 'aYFP', 'tau', 'tmRNA', 'a0', ...
    'KdTetR', 'nTetR', 'KdLacI', 'nLacI', ...
    'tmat1', 'tmat2', ...
    'kon', 'koff', ...
    'kdeg', 'Kdeg'};

K=size(params,2);
L=11;

PI=zeros(3,K,L);
varyparam=zeros(K,L);

for k=1
    newparams=params;
    varyparam(k,:)=logspace(log10(params(1,k)/10),log10(params(1,k)*10),L);
    for l=1:L
        newparams(1,k)=varyparam(k,l);
        for n=1:N
            [t,y] = ode23s(@(t,y)EMB_ODE_sensitivity(t, y, I, X, A(n,:), V(n,:), alphamax(n,:), newparams),tspan,y0,options);
            Y(n,:,:,:) = reshape(permute(y,[2,1]),I,X,T);
            
        end
        g = EMB_normalizeG(Y(:,[8,5],2:6,end));
        [PI(:,k,l),~] = EMB_g2piSGA(g,100);
    end
    k
    toc
end
toc

%% plots
figure(1)
    for k=1:K
        subplot(3,6,k)
            plot(varyparam(k,:),squeeze(PI(:,k,:)))
            box('on')
            ylabel('PI (bits)')
            xlabel(paraNames{1,k})
            xlim([0.95*varyparam(k,1),1.05*varyparam(k,L)]);
            set(gca,'XScale','log')
            xticks(varyparam(k,1:5:end))
    end
