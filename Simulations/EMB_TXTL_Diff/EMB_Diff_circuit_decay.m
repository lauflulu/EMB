function [t,outcome]= EMB_Diff_circuit_decay(a)

global  nspecies ncompartment ...
    D P A l_s l_ts V theta ...
    alphamax aLacI aRFP aTetR aYFP ...
    tRNA tmat1 tmat2 t1 t2 KdLacI nLacI KdTetR nTetR kon koff kdeg
nspecies = 8; % IPTG, free LacI, IPTG-LacI, RFP, RFPmat, TetR, YFP, YFPmat
ncompartment = 6; % 1 sender 5 receivers

D = 7.44e-10*60; % m2/min, free diffusion of IPTG
P = 1.67e-9*60; % m/min, for IPTG, determined from droplet fitting
V = 1.2e-12; % m3 average volume of compartments, 1 nL
%l = 150.3e-6; % m, average length of compartments
l_s = 2*(4*V./(3*pi)).^(1./3); % m, average length of compartments assuming
%sphere
theta = 42*pi/180; % interface angle

A = Bilayer_area(l_s, l_s, theta);
l_ts = 2.*(((4*V./(3*pi))+(sqrt(A./pi)).^3).^(1/3)-sqrt(A./pi)); % m, 
% average length of compartments, assuming truncated sphere

alphamax = 0.0036; % expression factor, determined from bulk and droplet fitting

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
    for icomp=1:1:ncompartment;
        x(ispecies, icomp)=0;
    end
end

%initial values < > 0
x(1,1)=a; %M, IPTG sender

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
[t,y]=ode23s(@EMB_Diff_pLacO_function,tspan,x0,options);

%output
t=t;
outcome=y(:,:);

end




function ydot = EMB_Diff_pLacO_function(t, y, flags)

global  nspecies ncompartment D P A l_ts V ...
    alphamax aLacI aRFP aTetR aYFP ...
    tRNA tmat1 tmat2 t1 t2 KdLacI nLacI KdTetR nTetR kon koff kdeg

% translate back:
ncount=size(y,1);

for icount=1:1:ncount
    icomp=floor((icount-1)/nspecies)+1;
    ispecies=icount-(icomp-1)*nspecies;
    x(ispecies,icomp)=y(icount);
end

dx=zeros(nspecies,ncompartment);


dx(1,1)=-(D*P*A)*(x(1,1)-x(1,2))./(V*(D+P*l_ts));

for icomp=2:1:ncompartment-1;
    dx(1,icomp)=-(D*P*A)*(2*x(1,icomp)-x(1,icomp-1)-x(1,icomp+1))./(V*(D+P*l_ts))-kon*x(1,icomp)*x(2,icomp)+koff*x(3,icomp);
end

dx(1,ncompartment)=-(D*P*A)*(x(1,ncompartment)-x(1,ncompartment-1))./(V*(D+P*l_ts))-kon*x(1,ncompartment)*x(2,ncompartment)+koff*x(3,ncompartment);

for icomp=2:1:ncompartment;
    if t < t1;
        dx(2,icomp)=alphamax*aLacI*(1-exp(-t/tRNA))*(1./(1+(x(6,icomp)./KdTetR).^nTetR))-kdeg*x(2,icomp)-kon*x(1,icomp)*x(2,icomp)+koff*x(3,icomp); % LacI
        dx(4,icomp)=alphamax*aRFP*(1-exp(-t/tRNA))*(1./(1+(x(6,icomp)./KdTetR).^nTetR))-x(4,icomp)/tmat1-kdeg*x(4,icomp); % RFP
        dx(6,icomp)=alphamax*aTetR*(1-exp(-t/tRNA))*(1./(1+(x(2,icomp)./KdLacI).^nLacI)); % TetR
        dx(7,icomp)=alphamax*aYFP*(1-exp(-t/tRNA))*(1./(1+(x(2,icomp)./KdLacI).^nLacI))-x(7,icomp)/tmat2;
    else;
        dx(2,icomp)=alphamax*aLacI*(1-exp(-t/tRNA))*(1./(1+(x(6,icomp)./KdTetR).^nTetR))*exp(-(t-t1)/t2)-kdeg*x(2,icomp)-kon*x(1,icomp)*x(2,icomp)+koff*x(3,icomp);
        dx(4,icomp)=alphamax*aRFP*(1-exp(-t/tRNA))*(1./(1+(x(6,icomp)./KdTetR).^nTetR))*exp(-(t-t1)/t2)-x(4,icomp)/tmat1-kdeg*x(4,icomp);
        dx(6,icomp)=alphamax*aTetR*(1-exp(-t/tRNA))*(1./(1+(x(2,icomp)./KdLacI).^nLacI))*exp(-(t-t1)/t2);
        dx(7,icomp)=alphamax*aYFP*(1-exp(-t/tRNA))*(1./(1+(x(2,icomp)./KdLacI).^nLacI))*exp(-(t-t1)/t2)-x(7,icomp)/tmat2;
    end
    dx(3,icomp)=kon*x(1,icomp)*x(2,icomp)-koff*x(3,icomp)-kdeg*x(3,icomp); % IPTG-LacI
    dx(5,icomp)=x(4,icomp)/tmat1-kdeg*x(5,icomp); % RFP mat
    dx(8,icomp)=x(7,icomp)/tmat2;
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

end