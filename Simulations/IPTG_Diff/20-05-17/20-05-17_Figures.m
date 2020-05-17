%% 20-05-17 Figures

%% Figure 1: Fitting P IPTG with Kd from 19-01-05 only

% experiment from 19-08-15, the minus a-HL plus IPTG data (Rpm)

mRpm_nonan=mRpm; %% all NaN are replaced with copies of the previous data point
R=RawDataTrace(t.*60, mRpm_nonan);
dR=diff(R, 8);
scharFit = Fit.Scharfit(@TXTL_Diff_fit_repression);

D_IPTG_m2s = 7.4442e-10;
D_IPTG_um2min=D_IPTG_m2s*1e12*60;

% radii of all droplets assuming a full sphere are measured with ImageJ
% multiplied by 2 and listed in l_s_measured_um
l_s_um=nanmean(nanmean(l_s_measured_um));
% which gives l_s_um = 152.7687;
V_calc_um3=(4./3).*pi.*((l_s_measured_um./2)).^(3);
V_um3=nanmean(V_calc_um3);

% problems
beta = 0;
V_um3 = 1.9008e6; % as calculated above
kon_Mmin = 7.2e6;
koff_min = 12.6;
n = 1.5; % from fit of 19-01-05
Kdtot_M = 2.9e-8; % from fit of 19-01-05
alphamax_uMmin = 0.03;
tRNA_min = 15; % also named t2
tmat_min = 7; 
I0_M = 0.003;  
R0_M = 1e-7;
t0_min = 1;  

% fitting for alpha in range 0 - 1e6, not log
% fitting for K in range 1e-3 - 1e6, not log
% scharSize: 10, fit also derivative, factor: 130
% fitting all data traces with one fit. 

scharFit.fit(dR);

% results: 

alpha_auuM = 275.9517;
K_um3min = 11073;

% calculating P

theta = 42*pi/180; % ?, interface angle
A_um2 = Bilayer_area(l_s_um, l_s_um, theta);
% which gives A_um2 = 8.2069e3; consistent with Lukas's fit
l_ts_um = 2.*(((4*V_um3./(3*pi))+((sqrt(A_um2./pi)+sqrt(A_um2./pi))./2).^3).^(1/3)-(sqrt(A_um2./pi)+sqrt(A_um2./pi))./2); 
% which gives l_ts_um = 93.7119;

P_ummin=K_um3min*D_IPTG_um2min/(D_IPTG_um2min*A_um2-K_um3min*l_ts_um);
P_ms=(P_ummin./60)*1e-6;

% Figure 1 shows the data and the fit, as well as permeabilities 1 order of
% magnitude away each


%% Figure 2: Diffusion of IPTG in assemblies

IPTG=[0.1,0.01,0.001,0.0001];
tic
[t_100mM, y_R_100mM]=EMB_Diff(IPTG(1));
[t_10mM, y_R_10mM]=EMB_Diff(IPTG(2));
[t_1mM, y_R_1mM]=EMB_Diff(IPTG(3));
[t_100uM, y_R_100uM]=EMB_Diff(IPTG(4));
toc


    [t,outcome]= EMB_Diff(a)
nspecies = 1; % IPTG
ncompartment = 6; % 1 sender 5 receivers
D = 7.44e-10*60*1e12; % um2/min, free diffusion of IPTG
P = 1.3531; % um/min, for IPTG, determined from droplet fitting
A = 8600; % um2, for unspecific diffusion
l = 68.83; % um, average length of compartments
V = 1.2e6; % um3 average volume of compartments, 1.2 nL

%initialization
for ispecies=1:1:nspecies;
    for icomp=1:1:ncompartment;
        x(ispecies, icomp)=0;
    end
end

%initial values < > 0
x(1,1)=a*1e6; %uM, IPTG sender

% create vector for odesolver
icount=0;

for icomp=1:1:ncompartment
    for ispecies=1:1:nspecies
        icount=icount+1;
        x0(icount)=x(ispecies,icomp);
    end
end

% duration of experiment
dt=10;
duration=480; % min, duration of experiment: 8h

% simulation run
tspan=[0:dt:duration];

options=odeset('AbsTol',1e-9,'RelTol',1e-5); % set tolerances
[t,y]=ode23s(@EMB_Diff_function,tspan,x0,options);

%output
t=t;
outcome=y(:,:);


ydot = EMB_Diff_function(t, y, flags)
% translate back:
ncount=size(y,1);

for icount=1:1:ncount
    icomp=floor((icount-1)/nspecies)+1;
    ispecies=icount-(icomp-1)*nspecies;
    x(ispecies,icomp)=y(icount);
end

dx=zeros(nspecies,ncompartment);


dx(1,1)=-(D*P*A)*(x(1,1)-x(1,2))./(V*(D+P*l));

for icomp=2:1:ncompartment-1;
dx(1,icomp)=-(D*P*A)*(2*x(1,icomp)-x(1,icomp-1)-x(1,icomp+1))./(V*(D+P*l));
end

dx(1,ncompartment)=-(D*P*A)*(x(1,ncompartment)-x(1,ncompartment-1))./(V*(D+P*l));


icount=0;
% construct derivative vector
for icomp=1:1:ncompartment
    for ispecies=1:1:nspecies
        icount=icount+1;
        dy(icount)=dx(ispecies,icomp);
    end
end

ydot=dy';


plot(t_100uM, y_R_100uM);
plot(t_1mM, y_R_1mM);
plot(t_10mM, y_R_10mM);
plot(t_100mM, y_R_100mM);

%% Figure 3: plottig diffusion profile for 10 mM source
plot([1,2,3,4,5,6], [y_R_10mM(1,:)]); % 0h
% 7 - 1 h
% 13 - 2 h
% 19 - 3 h
% 25 - 4 h
% 49 - 8 h

%% Figure 4: noise in V with new P

mu=1.2e-12;
sigmarel=0.35;
IPTG=[0.1,0.01,0.001,0.0001];

[t_100mM, y_R_100mM, r_100mM]=EMB_rand_V(mu, sigmarel, 100, IPTG(1));
[t_10mM, y_R_10mM, r_10mM]=EMB_rand_V(mu, sigmarel, 100, IPTG(2));
[t_1mM, y_R_1mM, r_1mM]=EMB_rand_V(mu, sigmarel, 100, IPTG(3));
[t_100uM, y_R_100uM, r_100uM]=EMB_rand_V(mu, sigmarel, 100, IPTG(4));

EMB_plotting_circuit(t_100uM, y_R_100uM);
EMB_plotting_circuit(t_1mM, y_R_1mM);
EMB_plotting_circuit(t_10mM, y_R_10mM);
EMB_plotting_circuit(t_100mM, y_R_100mM);








