%% 20-04-07 Figures

%% Figure 1-3: getting an idea of the effect of sender droplet volume on the
% dynamics of the system
[t, y_1pL]=EMB_Diff_pLacO_Vs(1e-15); 10 pL, 100 pL, 1 nL, 10 nL

[t,outcome]= EMB_Diff_pLacO_Vs(a)
nspecies = 5; % IPTG, free LacI, IPTG-LacI, P, Pmat
ncompartment = 6; % 1 sender 5 receivers
D = 7.44e-10*60; % m2/min, free diffusion of IPTG
P = 1.67e-9*60; % m/min, for IPTG, determined from droplet fitting
A = 2e-8; % m2, for unspecific diffusion
V = 1.2e-12; % m3 average volume of compartments, 1 nL
l = 2*(4*V./(3*pi)).^(1./3); % m, average length of compartments
Vs = a;
alphamax = 0.0036; % expression factor, determined from bulk and droplet fitting
tRNA = 15; % min, lifetime of RNA, from Karzbrun et al.
tmat = 7; % min, maturation time of GFP, from Iizuka et al.
KdLacI = 3.6e-8; % M, LacI binding to pLacO, determined from bulk fitting
nLacI = 1.35; % Hill coefficient, LacI binding to pLacO, determined from bulk fitting
kon = 7.2e6; % /M/min, on-rate of IPTG binding LacI, from Xu et al.
koff = 12.6; %/min, off-rate of IPTG binding LacI, from Xu et al.
R0 = 100e-9; % M, initial concentration of LacI in receivers

%initial values < > 0
x(1,1)=0.010; %M, IPTG sender
for i=1:1:ncompartment;
    x(2,i)=R0; % M, LacI in each receiver
end

% duration of experiment
dt=10;
duration=480; % min, duration of experiment: 8h

dx(1,1)=-(D*P*A)*(x(1,1)-x(1,2))./(Vs*(D+P*l));

% turns out unless it's much lower than 100 pL it doesnt change much.
% The concentration in the first droplets changes a lot, but in the others
% not really, and this has to do with the slow permeability. 

% Figure 1: IPTG
% Figure 2: IPTG-LacI which is basically the ratio of active pLacO
% Figure 3: protein expression

%% Figure 4: noise in droplet sender volume, how it affects the protein
% expression
% fixed the script to clean for negative values of V and avoid divergence
mu=1.2e-12;
sigmarel=0.35;
[t, y_R_200, r_200]=EMB_rand(mu, sigmarel, 200);
[t, y_R_100, r_100]=EMB_rand(mu, sigmarel, 100);
[t, y_R_20, r_20]=EMB_rand(mu, sigmarel, 20);
[t, y_R_20, r_20]=EMB_rand(mu, sigmarel, 1000);
[t, y_R_20, r_20]=EMB_rand(mu, sigmarel, 2000);
% it looks like 200 is still not enough to have reached a convergence, it's
% still different between 200 and 100. 
% Run with 1000: that's as good as 2000

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
    [t,outcome]=EMB_Diff_pLacO_Vs(r(i));
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

[t,outcome]= EMB_Diff_pLacO_Vs(a)
nspecies = 5; % IPTG, free LacI, IPTG-LacI, P, Pmat
ncompartment = 6; % 1 sender 5 receivers
D = 7.44e-10*60; % m2/min, free diffusion of IPTG
P = 1.67e-9*60; % m/min, for IPTG, determined from droplet fitting
A = 2e-8; % m2, for unspecific diffusion
V = 1.2e-12; % m3 average volume of compartments, 1 nL
l = 2*(4*V./(3*pi)).^(1./3); % m, average length of compartments
Vs = a;
alphamax = 0.0036; % expression factor, determined from bulk and droplet fitting
tRNA = 15; % min, lifetime of RNA, from Karzbrun et al.
tmat = 7; % min, maturation time of GFP, from Iizuka et al.
KdLacI = 3.6e-8; % M, LacI binding to pLacO, determined from bulk fitting
nLacI = 1.35; % Hill coefficient, LacI binding to pLacO, determined from bulk fitting
kon = 7.2e6; % /M/min, on-rate of IPTG binding LacI, from Xu et al.
koff = 12.6; %/min, off-rate of IPTG binding LacI, from Xu et al.
R0 = 100e-9; % M, initial concentration of LacI in receivers

%initial values < > 0
x(1,1)=0.010; %M, IPTG sender

for i=1:1:ncompartment;
    x(2,i)=R0; % M, LacI in each receiver
end

% duration of experiment
dt=10;
duration=480; % min, duration of experiment: 8h

ydot = EMB_Diff_pLacO_Vs_function(t, y, flags)
dx(1,1)=-(D*P*A)*(x(1,1)-x(1,2))./(Vs*(D+P*l));

EMB_plotting(t, y_R_1000); % that's figure 4

%% Figure 5: vary the volume of all droplets in an assembly
[t,y_R_20,r_20]=EMB_rand_V(mu, sigmarel, 20);
[t,y_R_1000,r_1000]=EMB_rand_V(mu, sigmarel, 1000);

[t,y_R, r]=EMB_rand_V(a,b,c)
mu = a; % mean of distribution
sigmarel = b; % CV of distribution
sigma = a*b; % stdev of distribution
ncompartment = 6; % number of compartments
runs = c.*ncompartment; % number of runs with random values
r = mu + sigma.*randn(runs,1); % random values within normal distribution
for i=1:runs;
    if r(i) <= 0;
        while r(i) <=0;
            r(i) = mu + sigma.*randn(1,1);
        end
    end
end
for i=1:c;
    [t,outcome]=EMB_Diff_pLacO_V(r(ncompartment*i-ncompartment+1:ncompartment*i));
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

[t,outcome]= EMB_Diff_pLacO_V(a)
nspecies = 5; % IPTG, free LacI, IPTG-LacI, P, Pmat
ncompartment = 6; % 1 sender 5 receivers
D = 7.44e-10*60; % m2/min, free diffusion of IPTG
P = 1.67e-9*60; % m/min, for IPTG, determined from droplet fitting
A = 2e-8; % m2, for unspecific diffusion
for i=1:6;
    V(i)=a(i);% m3 average volume of compartments, 1 nL
    l(i) = 2*(4*V(i)./(3*pi)).^(1./3); % m, average length of compartments
end
alphamax = 0.0036; % expression factor, determined from bulk and droplet fitting
tRNA = 15; % min, lifetime of RNA, from Karzbrun et al.
tmat = 7; % min, maturation time of GFP, from Iizuka et al.
KdLacI = 3.6e-8; % M, LacI binding to pLacO, determined from bulk fitting
nLacI = 1.35; % Hill coefficient, LacI binding to pLacO, determined from bulk fitting
kon = 7.2e6; % /M/min, on-rate of IPTG binding LacI, from Xu et al.
koff = 12.6; %/min, off-rate of IPTG binding LacI, from Xu et al.
R0 = 100e-9; % M, initial concentration of LacI in receivers
%initial values < > 0
x(1,1)=0.010; %M, IPTG sender
for i=1:1:ncompartment;
    x(2,i)=R0; % M, LacI in each receiver
end
% duration of experiment
dt=10;
duration=480; % min, duration of experiment: 8h

ydot = EMB_Diff_pLacO_V_function(t, y, flags)
dx(1,1)=-(D*P*A)*(x(1,1)-x(1,2))./(V(1)*(D+P*l(1)));
for icomp=2:1:ncompartment-1;
dx(1,icomp)=-(D*P*A)*(2*x(1,icomp)-x(1,icomp-1)-x(1,icomp+1))./(V(icomp)*(D+P*l(icomp)))-kon*x(1,icomp)*x(2,icomp)+koff*x(3,icomp);
end
dx(1,ncompartment)=-(D*P*A)*(x(1,ncompartment)-x(1,ncompartment-1))./(V(ncompartment)*(D+P*l(ncompartment)))-kon*x(1,ncompartment)*x(2,ncompartment)+koff*x(3,ncompartment);
