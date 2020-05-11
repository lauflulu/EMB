%% 20-04-15 Figures

%% Figure 1: Fitting 19-01-05 t1 and t2 with updated TXTL_decay

t=data(:,1)./60;
dataR=data(:,2:11);
R=RawDataTrace(t(1:201), dataR(1:201,:));
dR=diff(R,8);
model=Fit.Scharfit(@TXTL_repression_decay_fit);

beta = 3300; % (fitted by Igor)
Kd = 2.9e-8; % M
n = 1.5;
alphamax = 0.0036;
tRNA = 15; % min
tmat = 7;
R0 = 0, 1e-9, 3e-9, 1e-8 3e-8, 1e-7, 3e-7, 1e-6, 3e-6, 6e-6; % M
t0 = 1;

model.fit(dR);

alpha = 13000; 
t1 = 150; % min
t2 = 170; % min

%% Figure 2: Testing the whole circuit with decay

% Created EMB_Diff_circuit_decay and EMB_circuit_decay
[t,yR_0mM]=EMB_circuit_decay(0);
[t,yR_10uM]=EMB_circuit_decay(10e-6);
[t,yR_100uM]=EMB_circuit_decay(100e-6);
[t,yR_1mM]=EMB_circuit_decay(0.001);
[t,yR_10mM]=EMB_circuit_decay(0.01);
[t,yR_100mM]=EMB_circuit_decay(0.1);

[t,outcome]= EMB_circuit_decay(a)
nspecies = 8; % IPTG, free LacI, IPTG-LacI, RFP, RFPmat, TetR, YFP, YFPmat
alphamax = 0.0036; % expression factor, determined from bulk and droplet fitting
tRNA = 15; % min, lifetime of RNA, from Karzbrun et al.
tmat1 = 33; % min, maturation time of mScarletI, from Balleza et al. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5765880/pdf/nihms916443.pdf
tmat2 = 12; % min, maturation time of YPet, from Balleza et al. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5765880/pdf/nihms916443.pdf
t1 = 150; % min, delay before cell-extract degrades
t2 = 170; % min, characteristic time of cell-extract degradation

KdLacI = 3.6e-8; % M, LacI binding to pLacO, determined from bulk fitting
nLacI = 1.35; % Hill coefficient, LacI binding to pLacO, determined from bulk fitting
KdTetR = 1.3e-7; % M, determined from bulk fitting
nTetR = 4.3; % Hill coefficient, TetR binding to pTetO, determined from bulk fitting

kon = 7.2e6; % /M/min, on-rate of IPTG binding LacI, from Xu et al.
koff = 12.6; %/min, off-rate of IPTG binding LacI, from Xu et al.

%initial values < > 0
x0(1,1)=a; %M, IPTG sender

% duration of experiment
dt=10;
duration=480; % min, duration of experiment: 8h

plot(t./60, yR_0mM);

%% Figure 3: Compare different IPTG concentrations in bulk

subplot(2,2,1); hold on; axis square; plot(t./60, [yR_0uM(:,2),yR_0uM(:,5:6), yR_0uM(:,8)]);
subplot(2,2,2); hold on; axis square; plot(t./60, [yR_10uM(:,2),yR_10uM(:,5:6), yR_10uM(:,8)]);
subplot(2,2,3); hold on; axis square; plot(t./60, [yR_100uM(:,2),yR_100uM(:,5:6), yR_100uM(:,8)]);
subplot(2,2,4); hold on; axis square; plot(t./60, [yR_1mM(:,2),yR_1mM(:,5:6), yR_1mM(:,8)]);
