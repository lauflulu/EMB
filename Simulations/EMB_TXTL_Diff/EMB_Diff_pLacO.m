function [t,outcome]= EMB_Diff_pLacO(a)

global  nspecies ncompartment D P A l V alphamax tRNA tmat KdLacI nLacI kon koff R0
nspecies = 5; % IPTG, free LacI, IPTG-LacI, P, Pmat
ncompartment = 6; % 1 sender 5 receivers
D = 7.44e-10*60; % m2/min, free diffusion of IPTG
P = 1.67e-9*60; % m/min, for IPTG, determined from droplet fitting
A = 2e-8; % m2, for unspecific diffusion
V = 1.2e-12; % m3 average volume of compartments, 1 nL
%l = 150.3e-6; % m, average length of compartments
l = 2*(4*V./(3*pi)).^(1./3); % m, average length of compartments
alphamax = 0.0036; % expression factor, determined from bulk and droplet fitting
tRNA = 15; % min, lifetime of RNA, from Karzbrun et al.
tmat = 7; % min, maturation time of GFP, from Iizuka et al.
KdLacI = 3.6e-8; % M, LacI binding to pLacO, determined from bulk fitting
nLacI = 1.35; % Hill coefficient, LacI binding to pLacO, determined from bulk fitting
kon = 7.2e6; % /M/min, on-rate of IPTG binding LacI, from Xu et al.
koff = 12.6; %/min, off-rate of IPTG binding LacI, from Xu et al.
R0 = 100e-9; % M, initial concentration of LacI in receivers


%initialization
for ispecies=1:1:nspecies;
    for icomp=1:1:ncompartment;
        x(ispecies, icomp)=0;
    end
end

%initial values < > 0
x(1,1)=a; %M, IPTG sender

for i=1:1:ncompartment;
    x(2,i)=R0; % M, LacI in each receiver
end

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

global  nspecies ncompartment D P A l V alphamax tRNA tmat KdLacI nLacI kon koff

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
dx(1,icomp)=-(D*P*A)*(2*x(1,icomp)-x(1,icomp-1)-x(1,icomp+1))./(V*(D+P*l))-kon*x(1,icomp)*x(2,icomp)+koff*x(3,icomp);
end

dx(1,ncompartment)=-(D*P*A)*(x(1,ncompartment)-x(1,ncompartment-1))./(V*(D+P*l))-kon*x(1,ncompartment)*x(2,ncompartment)+koff*x(3,ncompartment);

for icomp=2:1:ncompartment;
    dx(2,icomp)=-kon*x(1,icomp)*x(2,icomp)+koff*x(3,icomp);
    dx(3,icomp)=kon*x(1,icomp)*x(2,icomp)-koff*x(3,icomp);
    dx(4,icomp)=alphamax*(1-exp(-t/tRNA))*(KdLacI.^nLacI./(KdLacI.^nLacI+x(2,icomp).^nLacI))-x(4,icomp)/tmat;
    dx(5,icomp)=x(4,icomp)/tmat;
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