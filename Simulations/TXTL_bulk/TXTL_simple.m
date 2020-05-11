function [t,outcome]=TXTL_simple

global nspecies alphamax tRNA tmat 
nspecies = 2; % P, Pmat
alphamax = 0.0036; % expression factor, determined from bulk and droplet fitting
tRNA = 15; % min, lifetime of RNA, from Karzbrun et al.
tmat = 7; % min, maturation time of GFP, from Iizuka et al.

% initialization
for ispecies=1:1:nspecies;
    x0(ispecies, 1)=0;
end

% duration of experiment
dt=10;
duration=480; % min, duration of experiment: 8h

% simulation run
tspan=[0:dt:duration];

options=odeset('AbsTol',1e-9,'RelTol',1e-5); % set tolerances
[t,y]=ode23s(@TXTL_simple_function,tspan,x0,options);

%output
t=t;
outcome=y(:,:);

end

function ydot=TXTL_simple_function(t,y,flags)

global nspecies alphamax tRNA tmat 

dy(1,1)=alphamax*(1-exp(-t/tRNA))-y(1,1)/tmat;
dy(2,1)=y(1,1)/tmat;

ydot=dy;

end
