function [t,outcome]= EMB_circuit_decay(a)

global  nspecies alphamax tRNA tmat1 tmat2 t1 t2 KdLacI nLacI KdTetR nTetR kon koff kdeg Kdeg
nspecies = 8; % IPTG, free LacI, IPTG-LacI, RFP, RFPmat, TetR, YFP, YFPmat

alphamax = 100e-12*60; % expression factor, around 100pM/s and up, from Schwarz-Schilling, Aufinger et al. 
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
kdeg = 15e-9; % M/min, rate of 0th-order protein degradation, from Karzbrun et al.
Kdeg = 1e-9; % M, deducted from Karzbrun et al.



%initialization
for ispecies=1:1:nspecies;
    x0(ispecies, 1)=0;
end

%initial values < > 0
x0(1,1)=a; %M, IPTG sender

% duration of experiment
dt=10;
duration=600; % min, duration of experiment: 8h

% simulation run
tspan=[0:dt:duration];

options=odeset('AbsTol',1e-9,'RelTol',1e-5); % set tolerances
[t,y]=ode23s(@EMB_circuit_decay_function,tspan,x0,options);

%output
t=t;
outcome=y(:,:);

end


function ydot = EMB_circuit_decay_function(t, y, flags)

global  nspecies alphamax tRNA tmat1 tmat2 t1 t2 KdLacI nLacI KdTetR nTetR kon koff

dy(1,1)=-kon*y(1,1)*y(2,1)+koff*y(3,1); % IPTG

if t < t1; % LacI
    dy(2,1)=alphamax*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-kon*y(1,1)*y(2,1)+koff*y(3,1);
else;
    dy(2,1)=alphamax*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))*exp(-(t-t1)/t2)-kon*y(1,1)*y(2,1)+koff*y(3,1);
end

dy(3,1)=kon*y(1,1)*y(2,1)-koff*y(3,1); % IPTG-LacI

if t < t1; % RFP
    dy(4,1)=alphamax*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-y(4,1)/tmat1;
else;
    dy(4,1)=alphamax*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))*exp(-(t-t1)/t2)-y(4,1)/tmat1;
end

dy(5,1)=y(4,1)/tmat1; % RFPmat

if t < t1; % TetR
    dy(6,1)=alphamax*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI));
else;
    dy(6,1)=alphamax*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))*exp(-(t-t1)/t2);
end

if t < t1; % YFP
    dy(7,1)=alphamax*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))-y(7,1)/tmat2;
else;
    dy(7,1)=alphamax*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))*exp(-(t-t1)/t2)-y(7,1)/tmat2;
end

dy(8,1)=y(7,1)/tmat2; % YFPmat

ydot=dy;

end