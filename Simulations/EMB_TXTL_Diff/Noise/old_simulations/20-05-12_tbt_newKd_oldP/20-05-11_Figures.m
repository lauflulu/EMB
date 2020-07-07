%% 20-05-11 Figures
% new Kd, old P, topology 1
mu=1.2e-12;
sigmarel=0.35;
IPTG=[100e-3, 10e-3, 1e-3, 100e-6];

[t_100mM, y_R_100mM, r_100mM]=EMB_rand_V(mu, sigmarel, 1000, IPTG(1));
[t_10mM, y_R_10mM, r_10mM]=EMB_rand_V(mu, sigmarel, 1000, IPTG(2));
[t_1mM, y_R_1mM, r_1mM]=EMB_rand_V(mu, sigmarel, 1000, IPTG(3));
[t_100uM, y_R_100uM, r_100uM]=EMB_rand_V(mu, sigmarel, 1000, IPTG(4));


[t,y_R, r]=EMB_rand_V(a,b,c,d)

tic 

mu = a; % mean of distribution
sigmarel = b; % CV of distribution
sigma = a*b; % stdev of distribution
ncompartment = 6; % number of compartments
runs = c.*ncompartment; % number of runs with random values
I0 = d; % IPTG, M

r = mu + sigma.*randn(runs,1); % random values within normal distribution

for i=1:runs;
    if r(i) <= 0;
        while r(i) <=0;
            r(i) = mu + sigma.*randn(1,1);
        end
    end
end

for i=1:c;
    [t,outcome]=EMB_Diff_circuit_decay_V(r(ncompartment*i-ncompartment+1:ncompartment*i), I0);
    [x1,y1]=size(outcome);
    for j=1:y1;
        y_R{j,1}(:,i)=outcome(:,j);
    end
end

[x2,y2]=size(y_R{1,1});

for j=1:y1;
    for i=1:x2;
        y_R{j,2}(i,1)=mean(y_R{j,1}(i,1:y2));
        y_R{j,3}(i,1)=std(y_R{j,1}(i,1:y2));
    end
end
    
t=t;
y_R=y_R;
r=r;

toc


    [t,outcome]= EMB_Diff_circuit_decay_V(a, b)

global  nspecies ncompartment ... 
    D P V l_s l_ts A theta ...
    alphamax aLacI aRFP aTetR aYFP tRNA tmat1 tmat2 t1 t2 ...
    KdLacI nLacI KdTetR nTetR kon koff kdeg Kdeg

nspecies = 8; % IPTG, free LacI, IPTG-LacI, RFP, RFPmat, TetR, YFP, YFPmat
ncompartment = 6; % 1 sender 5 receivers

D = 7.44e-10*60*1e12; % um2/min, free diffusion of IPTG
P = 1.67e-9*60*1e6; % um/min, for IPTG, determined from droplet fitting

for i=1:ncompartment;
    V(i)=a(i)*1e18;% um3 average volume of compartments, 1 nL
    l_s(i) = 2*(4*V(i)./(3*pi)).^(1./3); % m, average length of compartments
end
theta = 42*pi/180; % ?, interface angle

for i=1:ncompartment-1;
    A(i) = Bilayer_area(l_s(i), l_s(i+1), theta);
end

l_ts(1) = l_s(1)+(((4*V(1)./(3*pi))+(sqrt(A(1)./pi)).^3).^(1/3)-(sqrt(A(1)./pi))); % um, 
for i=2:ncompartment-1;
    l_ts(i) = 2.*(((4*V(i)./(3*pi))+((sqrt(A(i-1)./pi)+sqrt(A(i)./pi))./2).^3).^(1/3)-(sqrt(A(i-1)./pi)+sqrt(A(i)./pi))./2); % um, 
end
l_ts(ncompartment) = l_s(ncompartment)+(((4*V(ncompartment)./(3*pi))+(sqrt(A(ncompartment-1)./pi)).^3).^(1/3)-(sqrt(A(ncompartment-1)./pi))); % um, 
% average length of compartments, assuming truncated sphere

alphamax = 100e-6*60; % uM/min expression factor, around 100pM/s and up, from Schwarz-Schilling, Aufinger

aLacI = 3*1000/1077; % rNTPs length LacI, factor of expression normalized, and plasmid ratio
aRFP = 3*1000/693; % rNTPs length RFP
aTetR = 1000/618; % rNTPs length TetR
aYFP = 1000/717; % rNTPs length YFP

tRNA = 15; % min, lifetime of RNA, from Karzbrun et al.
tmat1 = 33; % min, maturation time of mScarletI, from Balleza et al. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5765880/pdf/nihms916443.pdf
tmat2 = 12; % min, maturation time of YPet, from Balleza et al. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5765880/pdf/nihms916443.pdf
t1 = 150; % min, delay before cell-extract degrades
t2 = 170; % min, characteristic time of cell-extract degradation

KdLacI = 2.9e-2; % uM, LacI binding to pLacO, determined from bulk fitting
nLacI = 1.5; % Hill coefficient, LacI binding to pLacO, determined from bulk fitting
KdTetR = 1.3e-1; % uM, determined from bulk fitting
nTetR = 4.3; % Hill coefficient, TetR binding to pTetO, determined from bulk fitting

kon = 7.2; % /uM/min, on-rate of IPTG binding LacI, from Xu et al.
koff = 12.6; %/min, off-rate of IPTG binding LacI, from Xu et al.

kdeg = 1e-3; % uM/min, rate of 0th-order protein degradation, from Karzbrun et al.
Kdeg = 1e-3; % uM, deducted from Karzbrun et al.


%initialization
for ispecies=1:1:nspecies;
    for icomp=1:1:ncompartment;
        x(ispecies, icomp)=0;
    end
end

%initial values < > 0
x(1,1)=b*1e6; %uM, IPTG sender

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

options=odeset('AbsTol',1e-8,'RelTol',1e-4); % set tolerances
[t,y]=ode23s(@EMB_Diff_circuit_decay_V_function,tspan,x0,options);

%output
t=t;
outcome=y(:,:);



ydot = EMB_Diff_circuit_decay_V_function(t, y, flags)

% translate back:
ncount=size(y,1);

for icount=1:1:ncount
    icomp=floor((icount-1)/nspecies)+1;
    ispecies=icount-(icomp-1)*nspecies;
    x(ispecies,icomp)=y(icount);
end

dx=zeros(nspecies,ncompartment);


dx(1,1)=-(D*P*A(1))*(x(1,1)-x(1,2))./(V(1)*(D+P*l_ts(1))); % IPTG

for icomp=2:1:ncompartment-1; 
    dx(1,icomp)=-((D*P*A(icomp-1))*(x(1,icomp)-x(1,icomp-1))+(D*P*A(icomp))*(x(1,icomp)-x(1,icomp+1)))./(V(icomp)*(D+P*l_ts(icomp)))-kon*x(1,icomp)*x(2,icomp)+koff*x(3,icomp);
end

dx(1,ncompartment)=-(D*P*A(ncompartment-1))*(x(1,ncompartment)-x(1,ncompartment-1))./(V(ncompartment)*(D+P*l_ts(ncompartment)))-kon*x(1,ncompartment)*x(2,ncompartment)+koff*x(3,ncompartment);

for icomp=2:1:ncompartment;
    if t < t1;
        dx(2,icomp)=alphamax*aLacI*(1-exp(-t/tRNA))*(1./(1+(x(6,icomp)./KdTetR).^nTetR))-kdeg/(1+Kdeg/x(2,icomp))-kon*x(1,icomp)*x(2,icomp)+koff*x(3,icomp); % LacI
        dx(3,icomp)=kon*x(1,icomp)*x(2,icomp)-koff*x(3,icomp)-kdeg/(1+Kdeg/x(3,icomp)); % IPTG-LacI
        dx(4,icomp)=alphamax*aRFP*(1-exp(-t/tRNA))*(1./(1+(x(6,icomp)./KdTetR).^nTetR))-x(4,icomp)/tmat1-kdeg/(1+Kdeg/x(4,icomp)); % RFP
        dx(5,icomp)=x(4,icomp)/tmat1-kdeg/(1+Kdeg/x(5,icomp)); % RFP mat
        dx(6,icomp)=alphamax*aTetR*(1-exp(-t/tRNA))*(1./(1+(x(2,icomp)./KdLacI).^nLacI)); % TetR
        dx(7,icomp)=alphamax*aYFP*(1-exp(-t/tRNA))*(1./(1+(x(2,icomp)./KdLacI).^nLacI))-x(7,icomp)/tmat2; % YFP
    else;
        dx(2,icomp)=(alphamax*aLacI*(1-exp(-t/tRNA))*(1./(1+(x(6,icomp)./KdTetR).^nTetR))-kdeg/(1+Kdeg/x(2,icomp)))*exp(-(t-t1)/t2)-kon*x(1,icomp)*x(2,icomp)+koff*x(3,icomp);
        dx(3,icomp)=kon*x(1,icomp)*x(2,icomp)-koff*x(3,icomp)-kdeg/(1+Kdeg/x(3,icomp))*exp(-(t-t1)/t2);
        dx(4,icomp)=(alphamax*aRFP*(1-exp(-t/tRNA))*(1./(1+(x(6,icomp)./KdTetR).^nTetR))-kdeg/(1+Kdeg/x(4,icomp)))*exp(-(t-t1)/t2)-x(4,icomp)/tmat1;
        dx(5,icomp)=x(4,icomp)/tmat1-kdeg/(1+Kdeg/x(5,icomp))*exp(-(t-t1)/t2);
        dx(6,icomp)=alphamax*aTetR*(1-exp(-t/tRNA))*(1./(1+(x(2,icomp)./KdLacI).^nLacI))*exp(-(t-t1)/t2);
        dx(7,icomp)=alphamax*aYFP*(1-exp(-t/tRNA))*(1./(1+(x(2,icomp)./KdLacI).^nLacI))*exp(-(t-t1)/t2)-x(7,icomp)/tmat2;
    end
    dx(8,icomp)=x(7,icomp)/tmat2; % YFP mat
end

icount=0;
% construct derivative vector
for icomp=1:1:ncompartment
    for ispecies=1:1:nspecies
        icount=icount+1;
        dy(icount)=dx(ispecies,icomp);
    end
end

ydot=dy';

