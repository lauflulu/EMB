%% 20-04-06 Figures

% Figure 1: trying randomization and figuring out what works
% created EMB_rand
% trying it for varying volume Vs of the sender droplet

mu=1.2e-12;
sigmarel=0.35; % this is CV. As reported by Lukas from an analysis of a 100 mM DP experiment

[t,y_R_10]=EMB_rand(mu, sigmarel, 10);, c = 100, c = 200

mu = a; % mean of distribution
sigmarel = b; % CV of distribution
sigma = a*b; % stdev of distribution
runs = c; % number of runs with random valus

r = mu + sigma.*randn(runs,1); % random values within normal distribution

for i=1:runs;
    [t,outcome]=EMB_Diff_rand_Vs(r(i));
    for j=1:6;
        y_R{j,1}(:,i)=outcome(:,j);
    end
end

[x,y]=size(y_R{1,1});

for j=1:6;
    for i=1:x;
        y_R{j,2}(i,1)=mean(y_R{j,1}(i,1:y));
        y_R{j,3}(i,1)=std(y_R{j,1}(i,1:y));
    end
end
    
t=t;
y_R=y_R;

end

function [t,outcome]= EMB_Diff_rand_Vs(a)
nspecies = 1; % IPTG
ncompartment = 6; % 1 sender 5 receivers
D = 7.44e-10*60; % m2/min, free diffusion of IPTG
P = 1.67e-9*60; % m/min, for IPTG, determined from droplet fitting
A = 2e-8; % m2, for unspecific diffusion
V = 1.2e-12; % m3 average volume of compartments, 1 nL
l = 2*(4*V./(3*pi)).^(1./3); % m, average length of compartments
Vs = a;

%initial values < > 0
x(1,1)=0.01; %M, IPTG sender

% duration of experiment
dt=10;
duration=480; % min, duration of experiment: 8h
%output
t=t;
outcome=y(:,:);

end


function ydot = EMB_Diff_rand_Vs_function(t, y, flags)

dx(1,1)=-(D*P*A)*(x(1,1)-x(1,2))./(Vs*(D+P*l));

for icomp=2:1:ncompartment-1;
dx(1,icomp)=-(D*P*A)*(2*x(1,icomp)-x(1,icomp-1)-x(1,icomp+1))./(V*(D+P*l));
end

dx(1,ncompartment)=-(D*P*A)*(x(1,ncompartment)-x(1,ncompartment-1))./(V*(D+P*l));

end

at c = 300 runs, the program gets buggy and the standard deviation of the simulated concentrations diverge. 
We compare 10, 100 and 200 runs, to confirm that 100 is roughly good enough. stdev is identical, mean is slightly higher at 200 runs. 
for i=1:6;
patch(t([1:end, end:-1:1]), [y_R_10{i,2}(1:end,1)+y_R_10{i,3}(1:end,1); y_R_10{i,2}(end:-1:1,1)-y_R_10{i,3}(end:-1:1,1)], [1, 0.2, 1]);
hold on; axis square;
plot(t, y_R_10{i,2}(:,1));
end

for i=1:6;
patch(t([1:end, end:-1:1]), [y_R_100{i,2}(1:end,1)+y_R_100{i,3}(1:end,1); y_R_100{i,2}(end:-1:1,1)-y_R_100{i,3}(end:-1:1,1)], [1, 0.2, 1]);
hold on; axis square;
plot(t, y_R_100{i,2}(:,1));
end

for i=1:6;
patch(t([1:end, end:-1:1]), [y_R_200{i,2}(1:end,1)+y_R_200{i,3}(1:end,1); y_R_200{i,2}(end:-1:1,1)-y_R_200{i,3}(end:-1:1,1)], [1, 0.2, 1]);
hold on; axis square;
plot(t, y_R_200{i,2}(:,1));
end

% Figure 2: gene expression in presence of noise in sender droplet volume




















