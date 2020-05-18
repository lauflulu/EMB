%% 20-04-09 Figures

%% Figure 1: parameter screen for A
A=[1e-10,1e-9,1e-8,1e-7,1e-6];
for i=1:5;
[t, y_R{i,1}]=EMB_Diff_pLacO_A(A(i));
end

[t,outcome]= EMB_Diff_pLacO_A(a)

nspecies = 5; % IPTG, free LacI, IPTG-LacI, P, Pmat
ncompartment = 6; % 1 sender 5 receivers
D = 7.44e-10*60; % m2/min, free diffusion of IPTG
V = 1.2e-12; % m3 average volume of compartments, 1 nL
%l = 150.3e-6; % m, average length of compartments
l = 2*(4*V./(3*pi)).^(1./3); % m, average length of compartments

P1 = 1.67e-9*60; % m/min, for IPTG, determined from droplet fitting with A1
A1 = 2e-8; % m2, for unspecific diffusion, from previous paper with larger droplets
A2 = 8.6e-9; % m2, determined by Lukas
P2 = A1*P1*D./(A2*(D+P1*l)-l*A1*P1);
P = P2;

A = a; % m2

alphamax = 0.0036; % expression factor, determined from bulk and droplet fitting
tRNA = 15; % min, lifetime of RNA, from Karzbrun et al.
tmat = 7; % min, maturation time of GFP, from Iizuka et al.
KdLacI = 3.6e-8; % M, LacI binding to pLacO, determined from bulk fitting
nLacI = 1.35; % Hill coefficient, LacI binding to pLacO, determined from bulk fitting
kon = 7.2e6; % /M/min, on-rate of IPTG binding LacI, from Xu et al.
koff = 12.6; %/min, off-rate of IPTG binding LacI, from Xu et al.
R0 = 100e-9; % M, initial concentration of LacI in receivers

%initial values < > 0
x(1,1)=0.01; %M, IPTG sender

for i=1:1:ncompartment;
    x(2,i)=R0; % M, LacI in each receiver
end

% duration of experiment
dt=10;
duration=480; % min, duration of experiment: 8h

EMB_plot_screen(t, y_R);
% it seems to affect expression quite significantly, particularly
% differentiation across the array
% but the variation looked at here is drastically higher than in real
% experiments

%% Figure 2: the effect of the variation in A as characterized by Lukas
mu = 8.6e-9;
sigmarel = 0.19;
[t, y_R_1000, r_1000]=EMB_rand(mu, sigmarel, 1000);
EMB_plotting(t, y_R_1000);

[t,y_R, r]=EMB_rand(a,b,c)
mu = a; % mean of distribution
sigmarel = b; % CV of distribution
sigma = a*b; % stdev of distribution
runs = c; % number of runs with random valus
r = mu + sigma.*randn(runs,1); % random values within normal distribution
for i=1:runs;
    if r(i) <= 0;
        while r(i) <=0;
            r(i) = mu + sigma.*randn(1,1);
        end
    end
end
for i=1:runs;
    [t,outcome]=EMB_Diff_pLacO_A(r(i));
    [x1,y1]=size(outcome);
    for j=1:y1;
        y_R{j,1}(:,i)=outcome(:,j);
    end
end
[x2,y2]=size(y_R{1,1});

for j=1:y1;
    for i=1:x2;
        y_R{j,2}(i,1)=mean(y_R{j,1}(i,1:y2));
        y_R{j,3}(i,1)=std(y_R{j,1}(i,1:y2));
    end
end
    
t=t;
y_R=y_R;
r=r;

[t,outcome]= EMB_Diff_pLacO_A(a)
nspecies = 5; % IPTG, free LacI, IPTG-LacI, P, Pmat
ncompartment = 6; % 1 sender 5 receivers
D = 7.44e-10*60; % m2/min, free diffusion of IPTG
V = 1.2e-12; % m3 average volume of compartments, 1 nL
%l = 150.3e-6; % m, average length of compartments
l = 2*(4*V./(3*pi)).^(1./3); % m, average length of compartments

P1 = 1.67e-9*60; % m/min, for IPTG, determined from droplet fitting with A1
A1 = 2e-8; % m2, for unspecific diffusion, from previous paper with larger droplets
A2 = 8.6e-9; % m2, determined by Lukas
P2 = A1*P1*D./(A2*(D+P1*l)-l*A1*P1);
P = P2;

A = a; % m2

alphamax = 0.0036; % expression factor, determined from bulk and droplet fitting
tRNA = 15; % min, lifetime of RNA, from Karzbrun et al.
tmat = 7; % min, maturation time of GFP, from Iizuka et al.
KdLacI = 3.6e-8; % M, LacI binding to pLacO, determined from bulk fitting
nLacI = 1.35; % Hill coefficient, LacI binding to pLacO, determined from bulk fitting
kon = 7.2e6; % /M/min, on-rate of IPTG binding LacI, from Xu et al.
koff = 12.6; %/min, off-rate of IPTG binding LacI, from Xu et al.
R0 = 100e-9; % M, initial concentration of LacI in receivers

%initial values < > 0
x(1,1)=0.01; %M, IPTG sender

for i=1:1:ncompartment;
    x(2,i)=R0; % M, LacI in each receiver
end
% duration of experiment
dt=10;
duration=480; % min, duration of experiment: 8h

%% Figure 3: CV when varying A
for i=1:49;
for j=1:30;
y_R_1000{j,4}(i,1)=y_R_1000{j,3}(i,1)./y_R_1000{j,2}(i,1);
end
end
EMB_plotting_CV(t, y_R_1000, y_R_1000);

