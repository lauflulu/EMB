%% 20-04-05

% Figure 1: showing the discrete diffusion
[t,outcome]=EMB_Diff(6), 5 receivers, 10 mM IPTG source

function [t,outcome]= EMB_Diff(a)
nspecies = 1; % IPTG
ncompartment = a; % 1 sender 5 receivers
D = 7.44e-10*60; % m2/min, free diffusion of IPTG
P = 1.67e-9*60; % m/min, for IPTG, determined from droplet fitting
A = 2e-8; % m2, for unspecific diffusion
l = 150.3e-6; % m, average length of compartments
V = 1e-12; % m3 average volume of compartments, 1 nL

%initial values < > 0
x(1,1)=0.01; %M, IPTG sender

% duration of experiment
dt=10;
duration=480; % min, duration of experiment: 8h


l = 150.3e-6;
x1=[0,l,l,2*l, 2*l, 3*l, 3*l, 4*l, 4*l, 5*l, 5*l, 6*l]
x1=x1*1e6;
j=[1,7,13,25,49];

for i=1:5;
plot(x2, [outcome(j(i),1), outcome(j(i),1), outcome(j(i),2), outcome(j(i),2), outcome(j(i),3), outcome(j(i),3), outcome(j(i),4), outcome(j(i),4), outcome(j(i),5), outcome(j(i),5), outcome(j(i),6), outcome(j(i),6)]);
end
plotting at 0h, 2h, 4h, 6h