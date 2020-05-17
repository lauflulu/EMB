%%20-03-29

%Figure 1: Diffusion of IPTG in a 5 receivers assembly for different source
%concentrations (identical to final experiments)
EMB_Diff(0.1), 0.01, 0.001, 0.0001

function [t,outcome]= EMB_Diff(a)
nspecies = 1; % IPTG
ncompartment = 6; % 1 sender 5 receivers
D = 7.44e-10*60; % m2/min, free diffusion of IPTG
P = 1.67e-9*60; % m/min, for IPTG, determined from droplet fitting
A = 2e-8; % m2, for unspecific diffusion
l = 150.3e-6; % m, average length of compartments
V = 1e-12; % m3 average volume of compartments, 1 nL

x(1,1)=a; %M, IPTG sender

% duration of experiment
dt=10;
duration=480; % min, duration of experiment: 8h

%Figure 2: Same as Figure 1 but source concentrations are compared for each
%receiver

%Figure 3: Diffusion of IPTG for different lengths of assemblies (5, 7 and
%9 receivers), with sender IPTG 10 mM (identical to 18-08-18 experiment)
EMB_Diff(6), 8, 10

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

% Figure 4: like Figure 1, but on a y-log scale, and a bar is added at
% 0.8e-6 M, the Kd of LacI to IPTG

% Figure 5: Introducing protein expression from a pLacO promoter with 100
% nM LacI in the receivers, 5 receivers, different IPTG concentrations
EMB_Diff_pLacO(0.0001), 0.001, 0.01, 0.1

function [t,outcome]= EMB_Diff_pLacO(a)
nspecies = 5; % IPTG, free LacI, IPTG-LacI, P, Pmat
ncompartment = 6; % 1 sender 5 receivers
D = 7.44e-10*60; % m2/min, free diffusion of IPTG
P = 1.67e-9*60; % m/min, for IPTG, determined from droplet fitting
A = 2e-8; % m2, for unspecific diffusion
l = 150.3e-6; % m, average length of compartments
V = 1e-12; % m3 average volume of compartments, 1 nL
alphamax = 0.0036; % expression factor, determined from bulk and droplet fitting
tRNA = 15; % min, lifetime of RNA, from Karzbrun et al.
tmat = 7; % min, maturation time of GFP, from Iizuka et al.
KdLacI = 3.6e-8; % M, LacI binding to pLacO, determined from bulk fitting
nLacI = 1.35; % Hill coefficient, LacI binding to pLacO, determined from bulk fitting
kon = 7.2e6; % /M/min, on-rate of IPTG binding LacI, from Xu et al.
koff = 12.6; %/min, off-rate of IPTG binding LacI, from Xu et al.
R0 = 100e-9; % M, initial concentration of LacI in receivers

%initial values < > 0
x(1,1)=a; %M, IPTG sender

for i=1:1:ncompartment;
    x(2,i)=R0; % M, LacI in each receiver
end

% duration of experiment
dt=10;
duration=480; % min, duration of experiment: 8h
