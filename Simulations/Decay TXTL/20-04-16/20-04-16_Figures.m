%% 20-04-16 Figures

%% Figure 1: EMB_2_circuit_decay taken plasmid ratio into consideration

IPTG=[0, 1e-6, 10e-6, 100e-6, 1e-3, 10e-3, 100e-3];

for i=1:7;
[t, yR{i,1}]=EMB_circuit_decay(IPTG(i));
end

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
duration=600; % min, duration of experiment: 8h

ydot = EMB_circuit_decay_function(t, y, flags)

dy(1,1)=-kon*y(1,1)*y(2,1)+koff*y(3,1); % IPTG

if t < t1; % LacI
    dy(2,1)=3*alphamax*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-kon*y(1,1)*y(2,1)+koff*y(3,1);
else;
    dy(2,1)=3*alphamax*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))*exp(-(t-t1)/t2)-kon*y(1,1)*y(2,1)+koff*y(3,1);
end

dy(3,1)=kon*y(1,1)*y(2,1)-koff*y(3,1); % IPTG-LacI

if t < t1; % RFP
    dy(4,1)=3*alphamax*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-y(4,1)/tmat1;
else;
    dy(4,1)=3*alphamax*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))*exp(-(t-t1)/t2)-y(4,1)/tmat1;
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

for i=1:7;
bar(i-0.2, yR{i,1}(61, 8)); % YFPmat
end
for i=1:7;
bar(i+0.2, yR{i,1}(61, 5)); % RFPmat
end

%% Figure 2: Comparing topologies with degradation 
% degradation building up kdeg/(1+Kdeg/P)
% alphamax replaced by DNA concentration with correct plasmid ratios

IPTG=[0, 1e-6, 10e-6, 100e-6, 1e-3, 10e-3, 100e-3];
for i=1:7;
[t, yR_1{i,1}]=EMB_1_circuit_decay(IPTG(i));
[t, yR_2{i,1}]=EMB_2_circuit_decay(IPTG(i));
[t, yR_3{i,1}]=EMB_3_circuit_decay(IPTG(i));
end

[t,outcome]= EMB_1_circuit_decay(a) % and other topologies accordingly
nspecies = 8; % IPTG, free LacI, IPTG-LacI, RFP, RFPmat, TetR, YFP, YFPmat

% alphamax = 0.0036; % expression factor, determined from bulk and droplet fitting
alphamax = 1e-9*100; % M, plasmid concentration*amplification from DNA to protein

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

Kdeg = 10e-9; % M, degradation thermodynamic constant, assumed from Karzbrun et al.
kdeg = 15e-9; % M/min, rate of 0th-order protein degradation, from Karzbrun et al.

%initial values < > 0
x0(1,1)=a; %M, IPTG sender

% duration of experiment
dt=10;
duration=600; % min, duration of experiment: 8h


ydot = EMB_1_circuit_decay_function(t, y, flags)

dy(1,1)=-kon*y(1,1)*y(2,1)+koff*y(3,1); % IPTG
if t < t1; % LacI
    dy(2,1)=3*alphamax*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-kdeg/(1+(Kdeg./y(2,1)))-kon*y(1,1)*y(2,1)+koff*y(3,1);
else;
    dy(2,1)=3*alphamax*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))*exp(-(t-t1)/t2)-kdeg/(1+(Kdeg./y(2,1)))-kon*y(1,1)*y(2,1)+koff*y(3,1);
end
dy(3,1)=kon*y(1,1)*y(2,1)-koff*y(3,1)-kdeg/(1+(Kdeg./y(3,1))); % IPTG-LacI
if t < t1; % RFP
    dy(4,1)=3*alphamax*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-y(4,1)/tmat1-kdeg/(1+(Kdeg./y(4,1)));
else;
    dy(4,1)=3*alphamax*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))*exp(-(t-t1)/t2)-y(4,1)/tmat1-kdeg/(1+(Kdeg./y(4,1)));
end
dy(5,1)=y(4,1)/tmat1-kdeg/(1+(Kdeg./y(5,1))); % RFPmat
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



subplot(1,3,1); hold on; axis square; plotyy(1,1,1,1);
subplot(1,3,2); hold on; axis square; plotyy(1,1,1,1);
subplot(1,3,3); hold on; axis square; plotyy(1,1,1,1);

for i=1:7;
bar(i-0.2, yR_1{i,1}(61, 8)); % YFPmat
end
for i=1:7;
bar(i+0.2, yR_1{i,1}(61, 5)); % RFPmat
end

for i=1:7;
bar(i-0.2, yR_2{i,1}(61, 8)); % YFPmat
end
for i=1:7;
bar(i+0.2, yR_2{i,1}(61, 5)); % RFPmat
end

for i=1:7;
bar(i-0.2, yR_3{i,1}(61, 7)); % YFPmat
end
for i=1:7;
bar(i+0.2, yR_3{i,1}(61, 5)); % RFPmat
end


%% Figure 3: Same as Figure 2 but with amplification proportional to aa length

for i=1:7;
[t, yR_1{i,1}]=EMB_1_circuit_decay(IPTG(i));
end

[t,outcome]= EMB_1_circuit_decay(a)
nspecies = 8; % IPTG, free LacI, IPTG-LacI, RFP, RFPmat, TetR, YFP, YFPmat

% alphamax = 0.0036; % expression factor, determined from bulk and droplet fitting
alphamax = 1e-9*100; % M, plasmid concentration*amplification from DNA to protein
lLacI = 1077; % rNTPs length LacI
lRFP = 693; % rNTPs length RFP
lTetR = 618; % rNTPs length TetR
lYFP = 717; % rNTPs length YFP
aLacI = (1/lLacI)*(4/(lLacI/3))*60; % max protein/DNA/min assuming 1 rNTP/s and 4 aa/s, from Karzbrun et al.
aRFP = (1/lRFP)*(4/(lRFP/3))*60;
aTetR = (1/lTetR)*(4/(lTetR/3))*60;
aYFP = (1/lYFP)*(4/(lYFP/3))*60;

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

Kdeg = 10e-9; % M, degradation thermodynamic constant, assumed from Karzbrun et al.
kdeg = 15e-9; % M/min, rate of 0th-order protein degradation, from Karzbrun et al.

%initial values < > 0
x0(1,1)=a; %M, IPTG sender

% duration of experiment
dt=10;
duration=600; % min, duration of experiment: 8h

ydot = EMB_1_circuit_decay_function(t, y, flags)

dy(1,1)=-kon*y(1,1)*y(2,1)+koff*y(3,1); % IPTG
if t < t1; % LacI
    dy(2,1)=3*alphamax*aLacI*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-kdeg/(1+(Kdeg./y(2,1)))-kon*y(1,1)*y(2,1)+koff*y(3,1);
else;
    dy(2,1)=3*alphamax*aLacI*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))*exp(-(t-t1)/t2)-kdeg/(1+(Kdeg./y(2,1)))-kon*y(1,1)*y(2,1)+koff*y(3,1);
end
dy(3,1)=kon*y(1,1)*y(2,1)-koff*y(3,1)-kdeg/(1+(Kdeg./y(3,1))); % IPTG-LacI
if t < t1; % RFP
    dy(4,1)=3*alphamax*aRFP*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-y(4,1)/tmat1-kdeg/(1+(Kdeg./y(4,1)));
else;
    dy(4,1)=3*alphamax*aRFP*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))*exp(-(t-t1)/t2)-y(4,1)/tmat1-kdeg/(1+(Kdeg./y(4,1)));
end
dy(5,1)=y(4,1)/tmat1-kdeg/(1+(Kdeg./y(5,1))); % RFPmat
if t < t1; % TetR
    dy(6,1)=alphamax*aTetR*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI));
else;
    dy(6,1)=alphamax*aTetR*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))*exp(-(t-t1)/t2);
end
if t < t1; % YFP
    dy(7,1)=alphamax*aYFP*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))-y(7,1)/tmat2;
else;
    dy(7,1)=alphamax*aYFP*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))*exp(-(t-t1)/t2)-y(7,1)/tmat2;
end
dy(8,1)=y(7,1)/tmat2; % YFPmat

% doesn't improve anything, it just changes the raw levels

