function [t,outcome]= EMB_Diff(a)

global  nspecies ncompartment D P A l V 
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

end


function ydot = EMB_Diff_function(t, y, flags)

global  nspecies ncompartment D P A l V 

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

end