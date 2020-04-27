%% 20-04-17 Figures

%% Figure 1: solving the bilayer area for 1 droplet of fixed volume and
% one droplet of variable volume (variable sender scenario)
mu=1.2e-12;
sigmarel=0.35;
[x1, x2, r]=Bilayer_radius_Vs(a,b,c)

tic 

mu sigmarel sigma runs 
mu = a; % mean of distribution
sigmarel = b; % CV of distribution
sigma = a*b; % stdev of distribution
runs = c; % number of runs with random valus

r = mu + sigma.*randn(runs,1); % random values within normal distribution

for i=1:runs;
    if r(i) <= 0;
        while r(i) <=0;
            r(i) = mu + sigma.*randn(1,1);
        end
    end
end

for i=1:runs;
    [xc, xd]=Bilayer_radius_Vs_function(r(i));
    x1(i,1)=xc;
    x2(i,1)=xd;
end
    
toc

[x1, x2]=Bilayer_radius_Vs_function(a)

global A B C V l_s theta ls_s Vs
V = 1.2e-12; % m3 average volume of compartments, 1 nL
%l = 150.3e-6; % m, average length of compartments
l_s = 2*(4*V./(3*pi)).^(1./3); % m, average length of compartments assuming
%sphere
theta = 60*pi/180; % ?, interface angle

Vs = a;
ls_s = 2*(4*Vs./(3*pi)).^(1./3); % m, average length of sender compartment
A=ls_s./2;
B=l_s./2;
C=sin(2*theta);

x1 = sqrt((A^4*B^2*C^2 + A^2*B^4*C^2 - 2*sqrt(-A^6*B^6*C^4*(C^2-1)))/(A^4 + 4*A^2*B^2*C^2 - 2*A^2*B^2 + B^4));
x2 = sqrt((A^4*B^2*C^2 + A^2*B^4*C^2 + 2*sqrt(-A^6*B^6*C^4*(C^2-1)))/(A^4 + 4*A^2*B^2*C^2 - 2*A^2*B^2 + B^4));

A1 = pi*x1^2;
A2 = pi*x2^2;


[x1, x2, r]=Bilayer_radius_Vs(mu, sigmarel, 1000);
A1 = pi*x1.^2;
A2 = pi*x2.^2;
Amu=8.6e-9;
Asigma=1.6e-9;
rA= Amu + Asigma.*randn(1000,1);


hist(A1, 100);
hold on; axis square;
hist(A2, 100);
hist(rA, 100);

%% Figure 2: solving the bilayer area for 2 droplets of varying volumes

mu=1.2e-12;
sigmarel=0.35;
[x1, x2, r1, r2]=Bilayer_radius(mu, sigmarel, 1000);
Amu=8.6e-9;
Asigma=1.6e-9;
rA= Amu + Asigma.*randn(1000,1);
A1 = pi*x1.^2;
A2 = pi*x2.^2;
hist(A1, 100);
hist(A2, 100);

[x1, x2, r1, r2]=Bilayer_radius(a,b,c)

tic 

global mu sigmarel sigma runs 
mu = a; % mean of distribution
sigmarel = b; % CV of distribution
sigma = a*b; % stdev of distribution
runs = c; % number of runs with random valus

r1 = mu + sigma.*randn(runs,1); % random values within normal distribution
r2 = mu + sigma.*randn(runs,1); % random values within normal distribution


for i=1:runs;
    if r1(i) <= 0;
        while r1(i) <=0;
            r1(i) = mu + sigma.*randn(1,1);
        end
    end
    if r2(i) <= 0;
        while r2(i) <=0;
            r2(i) = mu + sigma.*randn(1,1);
        end
    end
end

for i=1:runs;
    [xc, xd]=Bilayer_radius_function(r1(i), r2(i));
    x1(i,1)=xc;
    x2(i,1)=xd;
end
    
toc


[x1, x2]=Bilayer_radius_function(a, b)

global A B C V l_s theta ls_s Vs
V = a; % m3 average volume of compartments, 1 nL
%l = 150.3e-6; % m, average length of compartments
l_s = 2*(4*V./(3*pi)).^(1./3); % m, average length of compartments assuming
%sphere
theta = 60*pi/180; % ?, interface angle

Vs = b;
ls_s = 2*(4*Vs./(3*pi)).^(1./3); % m, average length of sender compartment
A=ls_s./2;
B=l_s./2;
C=sin(2*theta);

x1 = sqrt((A^4*B^2*C^2 + A^2*B^4*C^2 - 2*sqrt(-A^6*B^6*C^4*(C^2-1)))/(A^4 + 4*A^2*B^2*C^2 - 2*A^2*B^2 + B^4));
x2 = sqrt((A^4*B^2*C^2 + A^2*B^4*C^2 + 2*sqrt(-A^6*B^6*C^4*(C^2-1)))/(A^4 + 4*A^2*B^2*C^2 - 2*A^2*B^2 + B^4));

A1 = pi*x1^2;
A2 = pi*x2^2;

%% Figure 3: Back to studying EMB circuit modelling: with the real alphamax

IPTG=[0, 1e-6, 10e-6, 100e-6, 1e-3, 10e-3, 100e-3];
for i=1:7;
[t, y_R{i,1}]=EMB_circuit_decay(IPTG(i));
end

[t,outcome]= EMB_circuit_decay(a)

global  nspecies alphamax tRNA tmat1 tmat2 t1 t2 KdLacI nLacI KdTetR nTetR kon koff
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

ydot = EMB_circuit_decay_function(t, y, flags)

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

%% Figure 4: Comparing EMB topologies with first order degradation (it solves all the crazy negative or complex values)
IPTG=[0, 1e-6, 10e-6, 100e-6, 1e-3, 10e-3, 100e-3];

for i=1:7;
[t, yR_1{i,1}]=EMB_1_circuit_decay(IPTG(i));
[t, yR_2{i,1}]=EMB_2_circuit_decay(IPTG(i));
[t, yR_3{i,1}]=EMB_3_circuit_decay(IPTG(i));
end

[t,outcome]= EMB_1_circuit_decay(a) % and other topologies accordingly

global  nspecies alphamax tRNA tmat1 tmat2 t1 t2 KdLacI nLacI KdTetR nTetR kon koff kdeg
nspecies = 8; % IPTG, free LacI, IPTG-LacI, RFP, RFPmat, TetR, YFP, YFPmat

% alphamax = 0.0036; % expression factor, determined from bulk and droplet fitting
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


ydot = EMB_1_circuit_decay_function(t, y, flags)

global  nspecies alphamax tRNA tmat1 tmat2 t1 t2 KdLacI nLacI KdTetR nTetR kon koff kdeg

dy(1,1)=-kon*y(1,1)*y(2,1)+koff*y(3,1); % IPTG

if t < t1; 
    dy(2,1)=alphamax*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-kdeg*y(2,1)-kon*y(1,1)*y(2,1)+koff*y(3,1); % LacI
    dy(4,1)=alphamax*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-y(4,1)/tmat1-kdeg*y(4,1); % RFP
    dy(6,1)=alphamax*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI)); % TetR
    dy(7,1)=alphamax*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))-y(7,1)/tmat2; % YFP
else;
    dy(2,1)=alphamax*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))*exp(-(t-t1)/t2)-kdeg*y(2,1)-kon*y(1,1)*y(2,1)+koff*y(3,1);
    dy(4,1)=alphamax*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))*exp(-(t-t1)/t2)-y(4,1)/tmat1-kdeg*y(4,1);
    dy(6,1)=alphamax*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))*exp(-(t-t1)/t2);
    dy(7,1)=alphamax*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))*exp(-(t-t1)/t2)-y(7,1)/tmat2;
end

dy(3,1)=kon*y(1,1)*y(2,1)-koff*y(3,1)-kdeg*y(3,1); % IPTG-LacI

dy(5,1)=y(4,1)/tmat1-kdeg*y(5,1); % RFPmat

dy(8,1)=y(7,1)/tmat2; % YFPmat

ydot=dy;

%% Figure 5: now accounting for difference in plasmid ratio and aa length

IPTG=[0, 1e-6, 10e-6, 100e-6, 1e-3, 10e-3, 100e-3];

for i=1:7;
[t, yR_1{i,1}]=EMB_1_circuit_decay(IPTG(i));
[t, yR_2{i,1}]=EMB_2_circuit_decay(IPTG(i));
[t, yR_3{i,1}]=EMB_3_circuit_decay(IPTG(i));
end

[t,outcome]= EMB_1_circuit_decay(a) % and accordingly for 2 and 3

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

ydot = EMB_1_circuit_decay_function(t, y, flags)

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

%% Figure 6: Back to bilayer area, using a new graphical solution we found with Lukas

mu=1.2e-12;
sigmarel=0.35;
[x, A, r1, r2]=Bilayer_radius_2(mu, sigmarel, 10000);
Amu=8.6e-9;
Asigma=1.6e-9;
rA= Amu + Asigma.*randn(10000,1);
hist(A, 100);
hist(rA, 100);

[x, A, r1, r2]=Bilayer_radius_2(a,b,c)

tic 

global mu sigmarel sigma runs 
mu = a; % mean of distribution
sigmarel = b; % CV of distribution
sigma = a*b; % stdev of distribution
runs = c; % number of runs with random valus

r1 = mu + sigma.*randn(runs,1); % random values within normal distribution
r2 = mu + sigma.*randn(runs,1); % random values within normal distribution


for i=1:runs;
    if r1(i) <= 0;
        while r1(i) <=0;
            r1(i) = mu + sigma.*randn(1,1);
        end
    end
    if r2(i) <= 0;
        while r2(i) <=0;
            r2(i) = mu + sigma.*randn(1,1);
        end
    end
end

for i=1:runs;
    [x1, A1]=Bilayer_radius_2_function(r1(i), r2(i));
    x(i,1)=x1;
    A(i,1)=A1;
end
    
toc



[x1, A1]=Bilayer_radius_2_function(a, b)

global A B C D V l_s theta ls_s Vs
V = a; % m3 average volume of compartments, 1 nL
%l = 150.3e-6; % m, average length of compartments
l_s = 2*(4*V./(3*pi)).^(1./3); % m, average length of compartments assuming
%sphere
theta = 60*pi/180; % ?, interface angle

Vs = b;
ls_s = 2*(4*Vs./(3*pi)).^(1./3); % m, average length of sender compartment

A=ls_s./2;
B=l_s./2;
C=cos(2*theta);

D = sqrt(A^2+B^2+2*A*B*C);

x1 = sqrt(-D^4+2*D^2*A^2+2*D^2*B^2-A^4+2*A^2*B^2-B^4)/(2*D);

A1 = pi*x1.^2;

%% Figure 8: Back to EMB_circuit, comparing simulations to experiments
% importing simulations from Figure 6 and data from 19-10-09

%consider values at 10h
for i=1:6;
YFP_th(i,1)=yR_1{i,1}(61,8);
YFP_th(i,2)=yR_2{i,1}(61,8);
YFP_th(i,3)=yR_3{i,1}(61,7);
RFP_th(i,1)=yR_1{i,1}(61,5);
RFP_th(i,2)=yR_2{i,1}(61,5);
RFP_th(i,3)=yR_3{i,1}(61,5);
end

YFP_exp=[YFP(201,1:6); YFP(201,15:20); YFP(201, 8:13)]';
RFP_exp=[RFP(201,1:6); RFP(201,15:20); RFP(201, 8:13)]';

%defining a factor of proportionality
YFP_factor=YFP_th(6,1)/YFP_exp(6,1);
RFP_factor=RFP_th(6,1)/RFP_exp(6,1);

subplot(2,3,1); hold on; axis square; plot(IPTG_exp, YFP_exp(:,1)*YFP_factor); plot(IPTG_exp, YFP_th(:,1));
subplot(2,3,2); hold on; axis square; plot(IPTG_exp, YFP_exp(:,2)*YFP_factor); plot(IPTG_exp, YFP_th(:,2));
subplot(2,3,3); hold on; axis square; plot(IPTG_exp, YFP_exp(:,3)*YFP_factor); plot(IPTG_exp, YFP_th(:,3));
subplot(2,3,4); hold on; axis square; plot(IPTG_exp, RFP_exp(:,1)*RFP_factor); plot(IPTG_exp, RFP_th(:,1));
subplot(2,3,5); hold on; axis square; plot(IPTG_exp, RFP_exp(:,2)*RFP_factor); plot(IPTG_exp, RFP_th(:,2));
subplot(2,3,6); hold on; axis square; plot(IPTG_exp, RFP_exp(:,3)*RFP_factor); plot(IPTG_exp, RFP_th(:,3));

% IPTG = 1e-8 graphically corresponds to IPTG = 0 experimentally

%% Figure 8: get a smoother curve for the simulations

IPTG_th=[0, logspace(-6,-2)];


for i=1:51;
[t, outcome1]=EMB_1_circuit_decay(IPTG_th(i)); % same code as in Figure 6
YFP_th(i,1)=outcome1(61,8);
RFP_th(i,1)=outcome1(61,5);
[t, outcome2]=EMB_2_circuit_decay(IPTG_th(i));
YFP_th(i,2)=outcome2(61,8);
RFP_th(i,2)=outcome2(61,5);
[t, outcome3]=EMB_3_circuit_decay(IPTG_th(i));
YFP_th(i,3)=outcome3(61,7);
RFP_th(i,3)=outcome3(61,5);
end

IPTG_th_plot=[1e-8, logspace(-6,-2)];

plot(IPTG_th_plot, YFP_th(:,1));
plot(IPTG_th_plot, YFP_th(:,2));
plot(IPTG_th_plot, YFP_th(:,3));
plot(IPTG_th_plot, RFP_th(:,1));
plot(IPTG_th_plot, RFP_th(:,2));

