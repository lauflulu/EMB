%% 20-04-15 Figures

%% Figure 1: Fitting 19-01-05 t1 and t2 with updated TXTL_decay

t=data(:,1)./60;
dataR=data(:,2:11);
R=RawDataTrace(t(1:201), dataR(1:201,:));
dR=diff(R,8);
model=Fit.Scharfit(@TXTL_repression_decay_fit);

alpha = 210; % (fitted by Igor)
beta = 3300;
Kd = 2.9e-8; % M
n = 1.5;
alphamax = 0.0033;
tRNA = 15; % min
tmat = 7;
R0 = 0, 1e-9, 3e-9, 1e-8 3e-8, 1e-7, 3e-7, 1e-6, 3e-6, 6e-6; % M
t0 = 1;

model.fit(dR);

t1 = 1200; % min
t2 = 270; % min

%% Figure 2: Fitting 19-01-05 alpha and t2

t=data(:,1)./60;
dataR=data(:,2:11);
R=RawDataTrace(t(1:201), dataR(1:201,:));
dR=diff(R,8);
model=Fit.Scharfit(@TXTL_repression_decay_fit);

beta = 3300; % fitted by Igor
Kd = 2.9e-8; % M
n = 1.5;
alphamax = 0.0033;
tRNA = 15; % min
tmat = 7;
t1 = 0; % min
R0 = 0, 1e-9, 3e-9, 1e-8 3e-8, 1e-7, 3e-7, 1e-6, 3e-6, 6e-6; % M
t0 = 1;

model.fit(dR);

alpha = 1.9e4; 
t2 = 270; % min
