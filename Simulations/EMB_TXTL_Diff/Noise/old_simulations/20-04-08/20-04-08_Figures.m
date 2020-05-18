%% 20-04-08 Figures

%% Figure 1: comparing the CV in IPTG, active pLac O and reporter when Vs or all V are varied 
% importing the data y_RVs_1000, y_RV_1000

% calculating the CV
for i=1:49;
for j=1:30;
y_RV_1000{j,4}(i,1)=y_RV_1000{j,3}(i,1)./y_RV_1000{j,2}(i,1);
y_RVs_1000{j,4}(i,1)=y_RVs_1000{j,3}(i,1)./y_RVs_1000{j,2}(i,1);
end
end

EMB_plotting_CV(t, y_RVs_1000, y_RV_1000); (top is CV with Vs varied, bottom is CV with all V varied)

%% Figure 2: effect of variation in l
% calculated l with A1(Aurore), A2(Lukas), full sphere, truncated sphere
A1 = 2.0000e-08; % m2
A2 = 8.6000e-09; % m2
V1 = 1.0000e-12; % m3
V2 = 1.2000e-12; % m3

% EMB_Geo_param to calculate l1 and l2 sphere and tsphere

l1_sphere = 1.5030e-04; % m
l1_tsphere = 3.5808e-05; % m
l2_sphere = 1.5972e-04; % m
l2_tsphere = 6.8830e-05; % m

% the bilayer to bilayer length in microscope experiment (DP_10mM xy1_2) is
% definitely around 150 ?m and not 35 ?m

[t,y_R_l1s]=EMB_Diff_pLacO_l(l1_sphere);
[t,y_R_l1ts]=EMB_Diff_pLacO_l(l1_tsphere);
[t,y_R_l2ts]=EMB_Diff_pLacO_l(l2_tsphere);
[t,y_R_l2s]=EMB_Diff_pLacO_l(l2_sphere);

[t,outcome]= EMB_Diff_pLacO_l(a)
nspecies = 5; % IPTG, free LacI, IPTG-LacI, P, Pmat
ncompartment = 6; % 1 sender 5 receivers
D = 7.44e-10*60; % m2/min, free diffusion of IPTG
P = 1.67e-9*60; % m/min, for IPTG, determined from droplet fitting
A = 2e-8; % m2, for unspecific diffusion
V = 1.2e-12; % m3 average volume of compartments, 1 nL
l = a; % m, average length of compartments

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
dt=10;
duration=480; % min, duration of experiment: 8h

% Figure 1: plotted IPTG, active promoter (IPTG-LacI) and reporter
% first line is l1 sphere
% second line is l1 t sphere
% third line is l2 sphere
% four line is l2 t sphere

% it makes zero difference

%% Figure 3: effect of variation in l, systematically

l=[1e-6;10e-6;100e-6;1e-3;10e-3];
for i=1:5;
[t, y_R{i,1}]=EMB_Diff_pLacO_l(l(i)); % see above for script
end

% plotted IPTG, ratio of active pLacO (IPTG-LacI) and reporter

% basically l has zero influence because diffusion is permeation limited