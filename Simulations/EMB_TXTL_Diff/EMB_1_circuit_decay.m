function [t,outcome]= EMB_1_circuit_decay(a)

global  nspecies alphamax aLacI aRFP aTetR aYFP tRNA tmat1 tmat2 t1 t2 KdLacI nLacI KdTetR nTetR kon koff kdeg
nspecies = 8; % IPTG, free LacI, IPTG-LacI, RFP, RFPmat, TetR, YFP, YFPmat

% alphamax = 0.0036; % expression factor, determined from bulk and droplet fitting
alphamax = 100e-12*60; % expression factor, around 100pM/s and up, from Schwarz-Schilling, Aufinger et al.

aLacI = 3*1000/1077; % rNTPs length LacI, factor of expression normalized, and plasmid ratio
aRFP = 3*1000/693; % rNTPs length RFP
aTetR = 1000/618; % rNTPs length TetR
aYFP = 1000/717; % rNTPs length YFP

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
[t,y]=ode23s(@EMB_1_circuit_decay_function,tspan,x0,options);

%output
t=t;
outcome=y(:,:);

end


function ydot = EMB_1_circuit_decay_function(t, y, flags)

global  nspecies alphamax aLacI aRFP aTetR aYFP tRNA tmat1 tmat2 t1 t2 KdLacI nLacI KdTetR nTetR kon koff kdeg

dy(1,1)=-kon*y(1,1)*y(2,1)+koff*y(3,1); % IPTG

if t < t1; 
    dy(2,1)=alphamax*aLacI*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-kdeg*y(2,1)-kon*y(1,1)*y(2,1)+koff*y(3,1); % LacI
    dy(4,1)=alphamax*aRFP*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-y(4,1)/tmat1-kdeg*y(4,1); % RFP
    dy(6,1)=alphamax*aTetR*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI)); % TetR
    dy(7,1)=alphamax*aYFP*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))-y(7,1)/tmat2; % YFP
else;
    dy(2,1)=alphamax*aLacI*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))*exp(-(t-t1)/t2)-kdeg*y(2,1)-kon*y(1,1)*y(2,1)+koff*y(3,1);
    dy(4,1)=alphamax*aRFP*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))*exp(-(t-t1)/t2)-y(4,1)/tmat1-kdeg*y(4,1);
    dy(6,1)=alphamax*aTetR*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))*exp(-(t-t1)/t2);
    dy(7,1)=alphamax*aYFP*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))*exp(-(t-t1)/t2)-y(7,1)/tmat2;
end

dy(3,1)=kon*y(1,1)*y(2,1)-koff*y(3,1)-kdeg*y(3,1); % IPTG-LacI

dy(5,1)=y(4,1)/tmat1-kdeg*y(5,1); % RFPmat

dy(8,1)=y(7,1)/tmat2; % YFPmat

ydot=dy;

end